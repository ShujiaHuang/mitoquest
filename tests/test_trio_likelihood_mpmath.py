#!/usr/bin/env python3
"""Cross-check the C++ G-M-C trio likelihood against high-precision mpmath."""

import json
import pathlib
import random
import subprocess
import sys
import tempfile

try:
    import mpmath as mp
except ImportError:
    print("SKIP: mpmath is not installed for this Python interpreter")
    raise SystemExit(77)


mp.mp.dps = 80


def log_comb(n, k):
    return mp.loggamma(n + 1) - mp.loggamma(k + 1) - mp.loggamma(n - k + 1)


def log_beta(alpha, beta):
    return mp.loggamma(alpha) + mp.loggamma(beta) - mp.loggamma(alpha + beta)


def log_sum_exp(values):
    maximum = max(values)
    return maximum + mp.log(sum(mp.exp(value - maximum) for value in values))


def unsigned_stirling_first_kind(limit):
    coefficients = [[0] * (limit + 1) for _ in range(limit + 1)]
    coefficients[0][0] = 1
    for n in range(1, limit + 1):
        for r in range(1, n + 1):
            coefficients[n][r] = (
                coefficients[n - 1][r - 1] + (n - 1) * coefficients[n - 1][r]
            )
    return coefficients


def exact_log_likelihood(case):
    g_dp, g_alt, m_dp, m_alt, c_dp, c_alt, ne = case
    ne = mp.mpf(ne)
    pg = mp.mpf(g_alt) / g_dp

    all_reference = m_alt == 0 and c_alt == 0
    all_alternate = m_alt == m_dp and c_alt == c_dp
    if ne == 1 or pg == 0 or pg == 1:
        if all_reference and pg < 1:
            return mp.log(1 - pg)
        if all_alternate and pg > 0:
            return mp.log(pg)
        return -mp.inf

    s = ne - 1
    alpha_g = s * pg
    beta_g = s * (1 - pg)
    child_ref = c_dp - c_alt
    coefficients = unsigned_stirling_first_kind(c_dp)
    terms = []
    for r in range(c_alt + 1):
        for q in range(child_ref + 1):
            coefficient = coefficients[c_alt][r] * coefficients[child_ref][q]
            if coefficient:
                terms.append(
                    mp.log(coefficient)
                    + (r + q) * mp.log(s)
                    + log_beta(alpha_g + m_alt + r, beta_g + m_dp - m_alt + q)
                )

    return (
        log_comb(m_dp, m_alt)
        + log_comb(c_dp, c_alt)
        - log_beta(alpha_g, beta_g)
        - (mp.loggamma(s + c_dp) - mp.loggamma(s))
        + log_sum_exp(terms)
    )


def direct_integral_log_likelihood(case):
    g_dp, g_alt, m_dp, m_alt, c_dp, c_alt, ne = case
    ne = mp.mpf(ne)
    s = ne - 1
    pg = mp.mpf(g_alt) / g_dp
    alpha_g = s * pg
    beta_g = s * (1 - pg)

    def log_betabinom(p):
        return (
            log_comb(c_dp, c_alt)
            + log_beta(s * p + c_alt, s * (1 - p) + c_dp - c_alt)
            - log_beta(s * p, s * (1 - p))
        )

    def integrand(p):
        return mp.exp(
            log_comb(m_dp, m_alt)
            - log_beta(alpha_g, beta_g)
            + (alpha_g + m_alt - 1) * mp.log(p)
            + (beta_g + m_dp - m_alt - 1) * mp.log(1 - p)
            + log_betabinom(p)
        )

    return mp.log(mp.quad(integrand, [0, 1]))


def probe_log_likelihood(probe, case):
    command = [probe, *(str(value) for value in case)]
    output = subprocess.check_output(command, text=True).strip()
    return mp.mpf(output)


def cli_log_likelihood(mitoquest, case):
    g_dp, g_alt, m_dp, m_alt, c_dp, c_alt, ne = case
    with tempfile.TemporaryDirectory(prefix="mitoquest-trio-cli-") as directory:
        directory = pathlib.Path(directory)
        input_tsv = directory / "trio.tsv"
        output_json = directory / "result.json"
        input_tsv.write_text(
            "MOTHER_DP\tMOTHER_AD_ALT\tMOTHER_VAF\tCHILD_DP\tCHILD_AD_ALT\t"
            "QC\tHAS_G\tGRANDMOTHER_DP\tGRANDMOTHER_AD_ALT\n"
            f"{m_dp}\t{m_alt}\t{m_alt / m_dp:.17g}\t{c_dp}\t{c_alt}\tPASS\t"
            f"1\t{g_dp}\t{g_alt}\n"
        )
        subprocess.run(
            [
                mitoquest,
                "ne-estimate",
                "-i",
                str(input_tsv),
                "--model",
                "continuous",
                "--min-vaf",
                "0",
                "--max-vaf",
                "1",
                "--min-ne",
                str(ne),
                "--max-ne",
                str(ne),
                "-o",
                str(output_json),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        payload = "\n".join(
            line for line in output_json.read_text().splitlines()
            if not line.startswith("#")
        )
        return mp.mpf(str(json.loads(payload)["Max_Marginal_LogLik"]))


def assert_close(label, actual, expected, tolerance=mp.mpf("5e-10")):
    if mp.isinf(expected):
        if actual != expected:
            raise AssertionError(f"{label}: expected {expected}, received {actual}")
        return
    difference = abs(actual - expected)
    if difference > tolerance:
        raise AssertionError(
            f"{label}: C++={mp.nstr(actual, 20)}, mpmath={mp.nstr(expected, 20)}, "
            f"absolute difference={mp.nstr(difference, 8)}"
        )


def deterministic_random_cases():
    generator = random.Random(20260901)
    ne_values = ("1.001", "1.01", "1.1", "2", "5", "20", "100")
    cases = {}
    for index in range(24):
        g_dp = generator.randrange(20, 120)
        g_alt = generator.choice((1, g_dp - 1, generator.randrange(1, g_dp)))
        m_dp = generator.randrange(1, 35)
        m_alt = generator.choice((0, m_dp, generator.randrange(0, m_dp + 1)))
        c_dp = generator.randrange(1, 61)
        c_alt = generator.randrange(0, c_dp + 1)
        cases[f"random_{index + 1:02d}"] = (
            g_dp,
            g_alt,
            m_dp,
            m_alt,
            c_dp,
            c_alt,
            generator.choice(ne_values),
        )
    return cases


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: test_trio_likelihood_mpmath.py "
            "<mitoquest_trio_probe> <mitoquest>"
        )

    probe = sys.argv[1]
    mitoquest = sys.argv[2]
    cases = {
        "moderate": (100, 30, 50, 15, 80, 20, "10"),
        "near_lower_boundary": (100, 1, 8, 0, 40, 0, "1.01"),
        "near_upper_boundary": (100, 99, 8, 8, 40, 40, "1.01"),
        "balanced_depth": (250, 83, 40, 17, 100, 41, "2.5"),
        "higher_depth": (500, 173, 60, 21, 200, 77, "50"),
        "adaptive_high_depth": (500, 173, 60, 21, 300, 117, "50"),
        "founder_reference": (100, 0, 50, 0, 80, 0, "10"),
        "founder_alternate": (100, 100, 50, 50, 80, 80, "10"),
        "impossible_founder_reference": (100, 0, 50, 1, 80, 0, "10"),
        "ne_one_reference": (100, 30, 50, 0, 80, 0, "1"),
        "ne_one_impossible": (100, 30, 50, 15, 80, 20, "1"),
    }
    cases.update(deterministic_random_cases())

    for label, case in cases.items():
        expected = exact_log_likelihood(case)
        actual = probe_log_likelihood(probe, case)
        assert_close(label, actual, expected)
        print(f"{label}: OK")

    integral = direct_integral_log_likelihood(cases["moderate"])
    finite_sum = exact_log_likelihood(cases["moderate"])
    if abs(integral - finite_sum) > mp.mpf("1e-50"):
        raise AssertionError(
            "direct integration disagrees with finite-sum oracle: "
            f"{mp.nstr(abs(integral - finite_sum), 8)}"
        )
    print("moderate direct integral: OK")

    cli_actual = cli_log_likelihood(mitoquest, cases["moderate"])
    assert_close("production CLI moderate", cli_actual, finite_sum, mp.mpf("5e-8"))
    print("production CLI moderate: OK")


if __name__ == "__main__":
    main()