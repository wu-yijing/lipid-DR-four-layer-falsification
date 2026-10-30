#!/bin/bash
# Extract the 134 jointly clumped variants from the three GLGC files and from FinnGen R13.
set -euo pipefail
cd /e/workbuddy/_r1_review_20261029/mvmr

awk -F'\t' 'NR>1{print $3}' clump_union.clumps | sort -u > snp134.txt
wc -l < snp134.txt

# GLGC: SNP_hg18 SNP_hg19 rsid A1 A2 beta se N P Freq
for f in TG LDL HDL; do
  echo "== $f"
  zcat $f.txt.gz | awk -F'\t' 'NR==FNR{id[$1]=1;next} ($3 in id){print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10}' snp134.txt - > extract_$f.tsv
  wc -l < extract_$f.tsv
done

echo "== FinnGen"
# chrom pos ref alt rsids ... beta(9) sebeta(10) af_alt(11)
zcat /e/workbuddy/2026-09-04-18-19-16/phase1_out/finngen_R13_DM_RETINOPATHY.gz \
 | awk -F'\t' 'NR==FNR{id[$1]=1;next} {n=split($5,a,","); for(i=1;i<=n;i++){ if(a[i] in id){print a[i]"\t"$3"\t"$4"\t"$9"\t"$10"\t"$11; break} }}' snp134.txt - > extract_FG.tsv
wc -l < extract_FG.tsv
echo done
