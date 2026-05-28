import numpy as np
from tools.kimura import KimuraDistribution


def test_normalization_and_cdf_bounds():
    dist = KimuraDistribution(p=0.3, b=1.0)

    # PDF continuous integral (approx) + point masses should be ~1
    eps = 1e-6

    def cont_integrand(t):
        # continuous pdf excluding atoms
        return np.exp(dist.b * (t - dist.p)) / (dist._compute_normalization_constant() * t * (1 - t))

    integral, _ = __import__("scipy.integrate").integrate.quad(cont_integrand, eps, 1 - eps)
    total_prob = integral + float(dist.pdf(0)) + float(dist.pdf(1))

    assert 0.999 < total_prob < 1.001

    # CDF bounds
    assert 0.0 <= dist.cdf(0) <= 0.0 + 1e-12
    assert 0.9999 < dist.cdf(1) <= 1.0


def test_inverse_cdf_and_rvs():
    dist = KimuraDistribution(p=0.5, b=0.5)

    # random variates should be in [0,1]
    samples = dist.rvs(size=100, random_state=42)
    assert np.all(samples >= 0) and np.all(samples <= 1)

    # inverse cdf for some quantiles
    qs = [0.0, 0.05, 0.5, 0.95, 1.0]
    for q in qs:
        x = dist._inverse_cdf(q)
        assert 0.0 <= x <= 1.0


def test_pkimura_matches_r_example():
    from tools.kimura import pkimura

    p = 0.6
    b = 0.95
    intervals = np.arange(0, 1.1, 0.1)
    # Expected values from the R package README for pkimura(seq(0,1,0.1), p=0.6, b=0.95)
    expected = np.array([
        1.989732e-05, 2.165073e-05, 1.345668e-04, 3.560559e-03, 3.773170e-02,
        1.840206e-01, 4.908528e-01, 8.136989e-01, 9.711726e-01, 9.989777e-01,
        1.000000e+00
    ])
    computed = pkimura(intervals, p, b)
    assert np.allclose(computed, expected, rtol=1e-6, atol=1e-9)


def test_dkimura_atoms_match_readme():
    from tools.kimura import dkimura
    p = 0.6
    b = 0.95
    # Values from README
    assert abs(dkimura(0, p, b) - 1.989732e-05) < 1e-10
    assert abs(dkimura(1, p, b) - 9.694172e-06) < 1e-10
