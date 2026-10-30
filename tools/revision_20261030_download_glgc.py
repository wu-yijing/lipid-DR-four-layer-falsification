# -*- coding: utf-8 -*-
"""Parallel ranged downloader for the three GLGC 2013 summary-statistics files.

EBI throttles each connection to roughly 10 KB/s, so a single stream takes hours.
The files are served by Apache with byte-range support; splitting into chunks and
running several in parallel raises aggregate throughput substantially.
"""
import os, sys, time, gzip, shutil, threading
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000'
FILES = [
    ('TG.txt.gz', f'{BASE}/GCST002216/jointGwasMc_TG.txt.gz'),
    ('LDL.txt.gz', f'{BASE}/GCST002222/jointGwasMc_LDL%20(1).txt.gz'),
    ('HDL.txt.gz', f'{BASE}/GCST002223/jointGwasMc_HDL%20(1).txt.gz'),
]
WORK = os.path.dirname(os.path.abspath(__file__))
CHUNK = 4 * 1024 * 1024
THREADS = 10


def head_size(url):
    req = urllib.request.Request(url, method='HEAD')
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers['Content-Length'])


def fetch(url, start, end, path, tries=6):
    for k in range(tries):
        try:
            have = os.path.getsize(path) if os.path.exists(path) else 0
            want = end - start + 1
            if have >= want:
                return path
            req = urllib.request.Request(url, headers={'Range': f'bytes={start + have}-{end}'})
            with urllib.request.urlopen(req, timeout=120) as r, open(path, 'ab') as f:
                while True:
                    buf = r.read(65536)
                    if not buf:
                        break
                    f.write(buf)
            if os.path.getsize(path) >= want:
                return path
        except Exception as e:
            time.sleep(2 + 2 * k)
    raise RuntimeError(f'failed chunk {start}-{end} of {url}')


def one(name, url):
    size = head_size(url)
    parts = []
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = []
        for i, start in enumerate(range(0, size, CHUNK)):
            end = min(start + CHUNK - 1, size - 1)
            p = os.path.join(WORK, f'.{name}.part{i:03d}')
            parts.append((p, start, end))
            futs.append(ex.submit(fetch, url, start, end, p))
        total = 0
        for f in as_completed(futs):
            f.result()
            total += 1
    out = os.path.join(WORK, name)
    with open(out, 'wb') as w:
        for p, _, _ in parts:
            with open(p, 'rb') as r:
                shutil.copyfileobj(r, w, 1024 * 1024)
            os.remove(p)
    ok = os.path.getsize(out) == size
    print(f'[done] {name}: {size/1e6:.1f} MB, size_ok={ok}', flush=True)
    return ok


if __name__ == '__main__':
    for name, url in FILES:
        t0 = time.time()
        print(f'[start] {name}', flush=True)
        one(name, url)
        print(f'[time] {name}: {time.time()-t0:.0f}s', flush=True)
