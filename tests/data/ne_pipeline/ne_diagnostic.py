#!/usr/bin/env python3
"""
Diagnostic: Show the grid-alignment bias of the discrete BB + Bin MMLE
at high read depth, and verify the diffusion (continuous Beta) likelihood
has no such bias.
"""
import math
import random
from math import lgamma as gammaln
def betaln(a, b): return gammaln(a) + gammaln(b) - gammaln(a + b)

def log_comb(n, k):
    if k < 0 or k > n: return -math.inf
    return gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)

def log_bb_pmf(n, k, alpha, beta):
    if k < 0 or k > n: return -math.inf
    return (log_comb(n,k)
            + gammaln(k+alpha) + gammaln(n-k+beta) - gammaln(n+alpha+beta)
            - gammaln(alpha) - gammaln(beta) + gammaln(alpha+beta))

def log_bin_pmf(n, k, p):
    if k < 0 or k > n: return -math.inf
    if p <= 0: return 0.0 if k == 0 else -math.inf
    if p >= 1: return 0.0 if k == n else -math.inf
    return log_comb(n,k) + k*math.log(p) + (n-k)*math.log1p(-p)

def lse(a, b):
    if a == -math.inf: return b
    if b == -math.inf: return a
    if a > b:
        return a + math.log1p(math.exp(b-a))
    return b + math.log1p(math.exp(a-b))

def log_lik_discrete(m_alt, m_ref, c_alt, c_dp, Ne):
    alpha = m_alt + 1; beta = m_ref + 1
    total = -math.inf
    for k in range(Ne+1):
        lbb = log_bb_pmf(Ne, k, alpha, beta)
        p1 = k / Ne
        ls = log_bin_pmf(c_dp, c_alt, p1)
        total = lse(total, lbb + ls)
    return total

def log_lik_continuous(m_alt, m_ref, c_alt, c_dp, Ne, N=200):
    if Ne <= 1:
        return -math.inf
    a_m = m_alt + 1; b_m = m_ref + 1
    total = -math.inf
    log_norm_const = -betaln(a_m, b_m)
    for i in range(1, N):
        p = i / N
        log_pm_post = log_norm_const + (a_m-1)*math.log(p) + (b_m-1)*math.log1p(-p)
        ap = p * (Ne - 1); bp = (1.0 - p) * (Ne - 1)
        if ap <= 0 or bp <= 0: continue
        log_bb = log_bb_pmf(c_dp, c_alt, ap, bp)
        total = lse(total, log_pm_post + log_bb - math.log(N))
    return total

def simulate(n_pairs, true_ne, m_dp, c_dp, vaf_low, vaf_high, seed):
    rng = random.Random(seed)
    pairs = []
    for _ in range(n_pairs):
        p_m = rng.uniform(vaf_low, vaf_high)
        k = sum(1 for _ in range(true_ne) if rng.random() < p_m)
        p_c_true = k / true_ne
        m_alt = sum(1 for _ in range(m_dp) if rng.random() < p_m)
        c_alt = sum(1 for _ in range(c_dp) if rng.random() < p_c_true)
        pairs.append((m_alt, m_dp - m_alt, c_alt, c_dp))
    return pairs

def fit_ne(pairs, ne_min, ne_max, log_lik_fn):
    lls = {}
    for ne in range(ne_min, ne_max+1):
        ll = sum(log_lik_fn(*p, ne) for p in pairs)
        lls[ne] = ll
    best = max(lls, key=lls.get)
    return best, lls

def kimura_ne(pairs):
    num = 0.0; den = 0.0
    for m_alt, m_ref, c_alt, c_dp in pairs:
        m_dp = m_alt + m_ref
        if m_dp <= 0 or c_dp <= 0: continue
        pm = m_alt / m_dp; pc = c_alt / c_dp
        w = pm*(1-pm)
        if w <= 0: continue
        s = pm*(1-pm)/m_dp + pc*(1-pc)/c_dp
        d = (pc - pm)**2
        num += d - s
        den += w
    if den <= 0: return float('inf'), float('nan')
    V = max(0.0, num/den)
    b = max(1e-9, min(1-1e-9, 1-V))
    return 1/(1-b), b

print(f"{'true_Ne':>8} {'depth':>8} {'pairs':>6} {'MMLE_disc':>9} {'MMLE_cont':>9} {'Kimura':>8}")
for true_ne in [5, 10, 20, 50]:
    for m_dp, c_dp in [(500, 500), (2000, 2000)]:
        pairs = simulate(200, true_ne, m_dp, c_dp, 0.10, 0.90, seed=42)
        ne_disc, _ = fit_ne(pairs, 1, 150, log_lik_discrete)
        ne_cont, _ = fit_ne(pairs, 2, 150, log_lik_continuous)
        ne_kim, b = kimura_ne(pairs)
        print(f"{true_ne:8d} {m_dp:8d} {200:6d} {ne_disc:9d} {ne_cont:9d} {ne_kim:8.2f}")

# ----- Inject outliers (NUMTs / sequencing errors)
def simulate_with_outliers(n_pairs, true_ne, m_dp, c_dp, vaf_low, vaf_high,
                            outlier_frac, seed):
    rng = random.Random(seed)
    pairs = []
    for _ in range(n_pairs):
        p_m = rng.uniform(vaf_low, vaf_high)
        if rng.random() < outlier_frac:
            # outlier: artificial drift to 0 or 1
            p_c_true = rng.choice([0.0, 1.0])
        else:
            k = sum(1 for _ in range(true_ne) if rng.random() < p_m)
            p_c_true = k / true_ne
        m_alt = sum(1 for _ in range(m_dp) if rng.random() < p_m)
        c_alt = sum(1 for _ in range(c_dp) if rng.random() < p_c_true)
        pairs.append((m_alt, m_dp - m_alt, c_alt, c_dp))
    return pairs

def kimura_trimmed(pairs, trim_frac=0.10):
    """Robust Kimura: trim top trim_frac of pairs by |pc - pm|^2 / w."""
    rows = []
    for m_alt, m_ref, c_alt, c_dp in pairs:
        m_dp = m_alt + m_ref
        if m_dp <= 0 or c_dp <= 0: continue
        pm = m_alt / m_dp; pc = c_alt / c_dp
        w = pm*(1-pm)
        if w <= 0: continue
        s = pm*(1-pm)/m_dp + pc*(1-pc)/c_dp
        d = (pc - pm)**2
        rows.append((d - s, w, (d - s) / max(w, 1e-9)))
    if not rows: return float('inf'), float('nan')
    rows.sort(key=lambda x: x[2])
    keep = rows[: int(len(rows) * (1.0 - trim_frac))]
    num = sum(r[0] for r in keep); den = sum(r[1] for r in keep)
    if den <= 0: return float('inf'), float('nan')
    V = max(0.0, num/den)
    b = max(1e-9, min(1-1e-9, 1-V))
    return 1/(1-b), b


print()
print('--- user-scale: 442 pairs, depth 2000, varying true Ne and outlier rate ---')
print(f"{'true_Ne':>8} {'outl%':>6} {'MMLE_disc':>9} {'Kimura':>8} {'Kim_t10':>8} {'Kim_t20':>8}")
for true_ne in [20, 30, 50]:
    for outl in [0.0, 0.10, 0.20]:
        pairs = simulate_with_outliers(442, true_ne, 2000, 2000,
                                       0.10, 0.90, outl, seed=2026)
        ne_disc, _ = fit_ne(pairs, 1, 150, log_lik_discrete)
        ne_kim, _  = kimura_ne(pairs)
        ne_t10, _  = kimura_trimmed(pairs, 0.10)
        ne_t20, _  = kimura_trimmed(pairs, 0.20)
        print(f"{true_ne:8d} {outl*100:5.0f}% {ne_disc:9d} {ne_kim:8.2f}"
              f" {ne_t10:8.2f} {ne_t20:8.2f}")

print()
print('--- user-scale with light contamination + heteroplasmy bias ---')
# Add asymmetry: some pairs have child VAF biased to 0 (lost) more than to 1 (fixed),
# mimicking transmission of nearly-lost variants.
def simulate_realistic(n_pairs, true_ne, m_dp, c_dp,
                       loss_rate, seed):
    rng = random.Random(seed)
    pairs = []
    for _ in range(n_pairs):
        p_m = rng.uniform(0.05, 0.95)
        if rng.random() < loss_rate:
            p_c_true = 0.0  # lost in child
        else:
            k = sum(1 for _ in range(true_ne) if rng.random() < p_m)
            p_c_true = k / true_ne
        m_alt = sum(1 for _ in range(m_dp) if rng.random() < p_m)
        c_alt = sum(1 for _ in range(c_dp) if rng.random() < p_c_true)
        # Apply user-style filter: mother VAF in [0.10, 0.90]
        m_vaf = m_alt / m_dp
        if m_vaf < 0.10 or m_vaf > 0.90: continue
        pairs.append((m_alt, m_dp - m_alt, c_alt, c_dp))
    return pairs

print(f"{'true_Ne':>8} {'loss%':>6} {'kept':>5} {'MMLE_disc':>9} {'Kimura':>8} {'Kim_t10':>8} {'Kim_t20':>8}")
for true_ne in [30]:
    for loss in [0.0, 0.10, 0.20]:
        pairs = simulate_realistic(800, true_ne, 2000, 2000, loss, seed=2026)
        ne_disc, _ = fit_ne(pairs, 1, 150, log_lik_discrete)
        ne_kim, _  = kimura_ne(pairs)
        ne_t10, _  = kimura_trimmed(pairs, 0.10)
        ne_t20, _  = kimura_trimmed(pairs, 0.20)
        print(f"{true_ne:8d} {loss*100:5.0f}% {len(pairs):5d} {ne_disc:9d} {ne_kim:8.2f}"
              f" {ne_t10:8.2f} {ne_t20:8.2f}")

