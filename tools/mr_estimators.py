# -*- coding: utf-8 -*-
"""Two-sample MR estimators shared by the phase-11 analyses.

These reproduce exactly the estimators used in 1_mr/phase8_exmore_sensitivity.py
(IVW with multiplicative random-effects standard errors, MR-Egger, weighted
median), so a positive-control analysis is comparable with the primary one rather
than merely similar to it. Extracted here so the new scripts do not fork the
definitions.

Pure Python/NumPy; no TwoSampleMR dependency.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import stats

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement(a: str) -> str:
    return COMP.get(str(a).upper(), str(a).upper())


def align(ea, nea, ref, alt):
    """Return +1 / -1 if the exposure effect allele aligns with the outcome's
    reference/alternate allele, else None (ambiguous or mismatched)."""
    ea, nea, ref, alt = [str(x).upper() for x in (ea, nea, ref, alt)]
    if ea == alt and nea == ref:
        return +1
    if ea == ref and nea == alt:
        return -1
    if complement(ea) == alt and complement(nea) == ref:
        return +1
    if complement(ea) == ref and complement(nea) == alt:
        return -1
    return None


def ivw(b, se):
    """Inverse-variance weighted estimate with multiplicative random effects."""
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1 / se**2
    est = float(np.sum(w * b) / np.sum(w))
    se0 = float(np.sqrt(1 / np.sum(w)))
    Q = float(np.sum(w * (b - est) ** 2))
    df = len(b) - 1
    se_re = se0 * math.sqrt(max(1.0, Q / df)) if df > 0 else se0
    return est, se_re, Q, df


def egger(b, se, bx):
    """MR-Egger slope/intercept with their standard errors."""
    b, se, bx = np.asarray(b, float), np.asarray(se, float), np.asarray(bx, float)
    X = np.column_stack([np.ones(len(bx)), bx])
    W = np.diag(1 / se**2)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ b)
    resid = b - X @ beta
    cov = np.linalg.inv(X.T @ W @ X) * ((resid.T @ W @ resid) / (len(b) - 2))
    return beta[1], float(np.sqrt(cov[1, 1])), beta[0], float(np.sqrt(cov[0, 0]))


def weighted_median(b, se):
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1 / se**2
    o = np.argsort(b)
    b, w = b[o], w[o]
    p = np.cumsum(w) / np.sum(w) - 0.5 * w / np.sum(w)
    k = min(max(int(np.searchsorted(p, 0.5)), 1), len(b) - 1)
    return float(b[k - 1] + (b[k] - b[k - 1]) * (0.5 - p[k - 1]) / (p[k] - p[k - 1]))


def p_from_z(z):
    return float(2 * (1 - stats.norm.cdf(abs(z))))


def ci_from_est(est, se, z=1.959963984540054):
    return math.exp(est - z * se), math.exp(est + z * se)


def ci_p_consistent(lo, hi, p, tol=5e-3):
    """For a log-scale interval the point estimate is the geometric mean of the
    limits, so the p-value implied by the interval half-width must match the
    reported one. Returns (p_from_ci, consistent)."""
    if not (lo > 0 and hi > 0 and 0 < p < 1):
        return float("nan"), False
    se = (math.log(hi) - math.log(lo)) / (2 * 1.959963984540054)
    if se <= 0:
        return float("nan"), False
    z = abs(math.log(math.sqrt(lo * hi))) / se
    p_ci = p_from_z(z)
    return p_ci, abs(p_ci - p) <= tol


def heterogeneity(b, se):
    """Cochran's Q, its df, p-value and I-squared."""
    est, _, Q, df = ivw(b, se)
    p = float(1 - stats.chi2.cdf(Q, df)) if df > 0 else float("nan")
    i2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 and df > 0 else 0.0
    return Q, df, p, i2
