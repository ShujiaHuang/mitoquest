"""
kimura: Two-parameter Kimura distributions for mtDNA heteroplasmy analysis.

This package provides tools for fitting Kimura distributions to heteroplasmy 
data and testing for evidence of selection pressure in mitochondrial DNA.

References:
    Wonnapinij, P., Chinnery, P. F., & Samuels, D. C. (2008). 
    The distribution of mitochondrial DNA heteroplasmy due to random genetic drift. 
    The American Journal of Human Genetics, 83(5), 582-593.
    
    Kimura, M. (1955). Solution of a process of random genetic drift with a continuous 
    model. Proceedings of the National Academy of Sciences, 41(3), 144.

Author: Shujia Huang
Date: 2025-12-30

1. 面向对象设计:

- KimuraDistribution 类：核心分布类
- KimuraFitter 类：参数估计
- KimuraTest 类：假设检验

2. 主要功能:

- PDF、CDF 计算
- 随机数生成（逆变换采样）
- 最大似然估计参数拟合
- Monte Carlo Kolmogorov-Smirnov 检验

3. 向后兼容:

提供了类似 R 包的便捷函数 (dkimura, pkimura, rkimura, test_kimura)

4. 示例代码:

在 __main__ 中可以直接运行这个脚本查看示例输出。

"""
import numpy as np
from scipy import integrate, optimize, stats
from scipy.special import comb as _comb
from typing import Union, Tuple, Optional
import warnings


class KimuraDistribution:
    """
    Two-parameter Kimura distribution for modeling mtDNA heteroplasmy.
    
    The Kimura distribution is characterized by two parameters:
    - p: the initial frequency of the allele
    - b: the selection coefficient (b = 4*N*s where N is population size and s is selection)
    
    Attributes:
        p (float): Initial allele frequency, must be in (0, 1)
        b (float): Selection coefficient, must be positive
    """
    def __init__(self, p: float, b: float):
        """
        Initialize a Kimura distribution with given parameters.
        
        Args:
            p: Initial allele frequency, must be in (0, 1)
            b: Selection coefficient, must be positive
            
        Raises:
            ValueError: If parameters are out of valid range
        """
        if not (0 < p < 1):
            raise ValueError("Parameter p must be in the interval (0, 1)")
        if b <= 0:
            raise ValueError("Parameter b must be positive")
            
        self.p = p
        self.b = b
        self._normalization_constant = None
        
    def _compute_normalization_constant(self) -> float:
        """
        Compute the normalization constant for the Kimura distribution.
        
        Returns:
            The normalization constant
        """
        if self._normalization_constant is not None:
            return self._normalization_constant

        # Avoid singularities at the endpoints 0 and 1 by integrating on a trimmed interval
        eps = 1e-10
        def integrand(x):
            # integrand for the continuous part of the distribution
            return np.exp(self.b * (x - self.p)) / (x * (1 - x))

        # Integrate over (eps, 1 - eps); quad will not evaluate exactly at endpoints
        result, _ = integrate.quad(integrand, eps, 1 - eps, limit=200, epsabs=1e-9, epsrel=1e-9)
        self._normalization_constant = result
        return result
    
    def pdf(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Delegate to the R-equivalent `dkimura` implementation for consistency."""
        return dkimura(x, self.p, self.b)
    
    def cdf(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Delegate to the R-equivalent `pkimura` implementation for consistency."""
        return pkimura(x, self.p, self.b)
    
    def rvs(self, size: Optional[int] = None, random_state: Optional[int] = None) -> Union[float, np.ndarray]:
        """Delegate to R-equivalent `rkimura` for sampling and use a local RNG."""
        rng = np.random.default_rng(random_state)
        if size is None:
            return rkimura(1, self.p, self.b, random_state=rng)[0]
            
        return rkimura(size, self.p, self.b, random_state=rng)
    
    def _inverse_cdf(self, u: float) -> float:
        """
        Inverse CDF (quantile function) using root finding.
        
        Args:
            u: Probability value in [0, 1]
            
        Returns:
            Quantile corresponding to probability u
        """
        if u <= 0:
            return 0
        if u >= 1:
            return 1
        
        # Check if u falls in the point mass at 0
        prob_0 = self.pdf(0)
        if u <= prob_0:
            return 0
        
        # Check if u falls in the point mass at 1
        prob_1_complement = self.pdf(1)
        if u >= 1 - prob_1_complement:
            return 1
        
        # Use root finding for interior points
        def objective(x):
            return self.cdf(x) - u
        
        try:
            root = optimize.brentq(objective, 1e-10, 1 - 1e-10)
            # brentq may sometimes return a (root, result) tuple if full_output=True,
            # or a result-like object; normalize to a numeric root value.
            if isinstance(root, tuple):
                root_val = root[0]
            elif hasattr(root, 'root'):
                root_val = root.root
            else:
                root_val = root
            return float(root_val)
        except ValueError:
            # Fallback to bisection if brentq fails; use root_scalar and return the numeric root
            res = optimize.root_scalar(objective, bracket=[1e-10, 1 - 1e-10], method='bisect')
            if getattr(res, 'converged', False):
                return float(res.root)
            else:
                raise ValueError("Inverse CDF root finding failed")
    
    def __repr__(self) -> str:
        return f"KimuraDistribution(p={self.p:.4f}, b={self.b:.4f})"


class KimuraFitter:
    """
    Fit Kimura distribution parameters to observed heteroplasmy data.
    """
    
    @staticmethod
    def fit(data: np.ndarray, method: str = 'mle') -> Tuple[float, float]:
        """
        Fit Kimura distribution parameters to data.
        
        Args:
            data: Array of heteroplasmy values in [0, 1]
            method: Fitting method ('mle' for maximum likelihood estimation)
            
        Returns:
            Tuple of (p, b) estimated parameters
            
        Raises:
            ValueError: If data is invalid
        """
        data = np.asarray(data)
        
        if np.any((data < 0) | (data > 1)):
            raise ValueError("All data values must be in [0, 1]")
        
        if len(data) == 0:
            raise ValueError("Data array cannot be empty")
        
        if method == 'mle':
            return KimuraFitter._fit_mle(data)
        else:
            raise ValueError(f"Unknown fitting method: {method}")
    
    @staticmethod
    def _fit_mle(data: np.ndarray) -> Tuple[float, float]:
        """
        Maximum likelihood estimation of Kimura parameters.
        
        Args:
            data: Array of heteroplasmy values
            
        Returns:
            Tuple of (p, b) estimated parameters
        """
        # Initial guess: p = mean of data, b = 1
        p0 = np.mean(data)
        p0 = np.clip(p0, 0.01, 0.99)  # Ensure valid range
        b0 = 1.0
        
        def neg_log_likelihood(params):
            p, b = params
            if p <= 0 or p >= 1 or b <= 0:
                return np.inf
            
            try:
                dist = KimuraDistribution(p, b)
                pdf_values = dist.pdf(data)
                
                # Avoid log(0)
                pdf_values = np.maximum(pdf_values, 1e-300)
                
                return -np.sum(np.log(pdf_values))
            except Exception:
                return np.inf
        
        # Use optimization to find MLE
        result = optimize.minimize(
            neg_log_likelihood,
            x0=[p0, b0],
            method='Nelder-Mead',
            bounds=[(0.01, 0.99), (0.01, 100)]
        )
        
        if result.success:
            return result.x[0], result.x[1]
        else:
            warnings.warn("MLE optimization did not converge, using moment estimates")
            return float(p0), float(b0)


class KimuraTest:
    """
    Hypothesis test for selection pressure in mtDNA heteroplasmy distributions.
    
    Performs a Kolmogorov-Smirnov test to assess whether observed heteroplasmy
    data is consistent with a Kimura distribution (neutral evolution with drift)
    or shows evidence of selection.
    """
    
    @staticmethod
    def test(data: np.ndarray, n_simulations: int = 10, 
             random_state: Optional[int] = None) -> dict:
        """
        Perform Monte Carlo Kolmogorov-Smirnov test for Kimura distribution fit.
        
        Args:
            data: Array of heteroplasmy values in [0, 1]
            n_simulations: Number of Monte Carlo simulations for p-value estimation
            random_state: Random seed for reproducibility
            
        Returns:
            Dictionary containing test results:
                - 'statistic': KS test statistic
                - 'p_value': Monte Carlo p-value
                - 'p_hat': Estimated p parameter
                - 'b_hat': Estimated b parameter
                - 'alternative': Test alternative hypothesis
        """
        data = np.asarray(data)
        
        if np.any((data < 0) | (data > 1)):
            raise ValueError("All data values must be in [0, 1]")
        
        # Fit Kimura distribution to data
        p_hat, b_hat = KimuraFitter.fit(data)
        
        # Compute KS statistic
        fitted_dist = KimuraDistribution(p_hat, b_hat)
        ks_statistic = KimuraTest._compute_ks_statistic(data, fitted_dist)
        
        # Monte Carlo simulation for p-value
        if random_state is not None:
            np.random.seed(random_state)
        
        simulated_statistics = []
        n = len(data)
        
        for _ in range(n_simulations):
            # Generate sample from fitted distribution
            sim_data = fitted_dist.rvs(size=n)

            # Ensure sim_data is an ndarray before fitting (rvs may return scalar)
            sim_data = np.atleast_1d(sim_data)

            # Fit distribution to simulated data
            try:
                p_sim, b_sim = KimuraFitter.fit(sim_data)
                sim_dist = KimuraDistribution(p_sim, b_sim)
                sim_ks = KimuraTest._compute_ks_statistic(sim_data, sim_dist)
                simulated_statistics.append(sim_ks)
            except Exception:
                # If fitting fails, skip this simulation
                continue
        
        # Compute p-value
        simulated_statistics = np.array(simulated_statistics)
        p_value = np.mean(simulated_statistics >= ks_statistic)
        
        return {
            'statistic': ks_statistic,
            'p_value': p_value,
            'p_hat': p_hat,
            'b_hat': b_hat,
            'alternative': 'one-sided',
            'method': 'Monte Carlo Kolmogorov-Smirnov'
        }
    
    @staticmethod
    def _compute_ks_statistic(data: np.ndarray, distribution: KimuraDistribution) -> float:
        """
        Compute Kolmogorov-Smirnov test statistic.
        
        Args:
            data: Observed data
            distribution: Fitted Kimura distribution
            
        Returns:
            KS test statistic (maximum absolute difference between ECDFs)
        """
        n = len(data)
        sorted_data = np.sort(data)
        
        # Empirical CDF
        ecdf = np.arange(1, n + 1) / n
        
        # Theoretical CDF
        tcdf = distribution.cdf(sorted_data)
        
        # KS statistic
        ks_stat = np.max(np.abs(ecdf - tcdf))
        
        return ks_stat


def print_test_results(results: dict):
    """
    Pretty print the results of a Kimura test.
    
    Args:
        results: Dictionary returned by KimuraTest.test()
    """
    print(f"\n{results['method']}\n")
    print(f"data:  Kimura({results['p_hat']:.3f}, {results['b_hat']:.4f})")
    print(f"D = {results['statistic']:.5f}, b = {results['b_hat']:.5f}, p-value = {results['p_value']:.3f}")
    print(f"alternative hypothesis: {results['alternative']}")


# Ported functions from the R package to match behavior and numerical results

def _hypgeo(i: int, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Recursive computation of F(i-1, i+2, 2, x) as in the R implementation."""
    # Support scalar or array x
    scalar = np.isscalar(x)
    xs = np.atleast_1d(x)
    out = np.empty_like(xs, dtype=float)
    for idx, xv in enumerate(xs):
        a = 1 - i
        b = i + 2
        c = 2
        f0 = 1.0
        f1 = 1.0 - (b * xv) / c
        if i == 1:
            out[idx] = f0
            continue
        if i == 2:
            out[idx] = f1
            continue
        f_prev = f0
        f_curr = f1
        f = 0.0
        for j in range(2, i):
            aa = 1 - j
            f = (aa * (1 - xv) * (f_curr - f_prev) + (aa + b * xv - c) * f_curr) / (aa - c)
            f_prev = f_curr
            f_curr = f
        out[idx] = f
    return float(out[0]) if scalar else out


def _f0(p: float, b: float) -> float:
    q = 1 - p
    N = 250
    sum_terms = np.zeros(N, dtype=float)
    signs = np.tile(np.array([-1.0, 1.0]), int(np.ceil(N / 2)))[:N]
    i_arr = np.arange(1, N + 1)
    i_coef = p * q * (2 * i_arr + 1) * signs
    b_coef = b ** _comb(np.arange(2, N + 2), 2)
    last_idx = 0
    for i in range(1, N + 1):
        q_hg = _hypgeo(i, q)
        sum_terms[i - 1] = q_hg * i_coef[i - 1] * b_coef[i - 1]
        last_idx = i
        if i > 1 and abs(sum_terms[i - 1] - sum_terms[i - 2]) < 1e-4:
            break

    if abs(sum_terms[N - 2] - sum_terms[N - 1]) > 1e-4:
        warnings.warn("Series not yet converged at N = 250.")

    # Use last computed index for safe slicing; if none computed, fall back to N
    if last_idx == 0:
        last_idx = N
    return max(q + np.sum(sum_terms[:last_idx]), 0.0)


def _f1(p: float, b: float) -> float:
    q = 1 - p
    N = 250
    sum_terms = np.zeros(N, dtype=float)
    i_arr = np.arange(1, N + 1)
    i_coef = p * q * (2 * i_arr + 1) * np.tile(np.array([-1.0, 1.0]), int(np.ceil(N / 2)))[:N]
    b_coef = b ** _comb(np.arange(2, N + 2), 2)
    last_idx = 0
    for i in range(1, N + 1):
        p_hg = _hypgeo(i, p)
        sum_terms[i - 1] = p_hg * i_coef[i - 1] * b_coef[i - 1]
        last_idx = i
        if i > 1 and abs(sum_terms[i - 1] - sum_terms[i - 2]) < 1e-4:
            break
    if abs(sum_terms[N - 2] - sum_terms[N - 1]) > 1e-4:
        warnings.warn("Series not yet converged at N = 250.")

    # Use last computed index for safe slicing; if none computed, fall back to N
    if last_idx == 0:
        last_idx = N
    return max(p + np.sum(sum_terms[:last_idx]), 0.0)


def _phi(x: Union[float, np.ndarray], p: float, b: float) -> Union[float, np.ndarray]:
    q = 1 - p
    N = 250
    x_arr = np.atleast_1d(x)
    lz = len(x_arr)
    sum_terms = np.zeros((N, lz), dtype=float)
    i_arr = np.arange(1, N + 1)
    i_coef = p * q * (i_arr * (i_arr + 1) * (2 * i_arr + 1))
    b_coef = b ** _comb(np.arange(2, N + 2), 2)

    # first term
    p_hg = _hypgeo(1, p)
    z_hg = np.array([_hypgeo(1, xv) for xv in x_arr])
    sum_terms[0, :] = p_hg * z_hg * i_coef[0] * b_coef[0]

    j_idx = np.arange(lz)
    last_idx = 1
    for i in range(2, N + 1):
        z_hg = np.array([_hypgeo(i, x_arr[j]) for j in j_idx])
        p_hg = _hypgeo(i, p)
        sum_terms[i - 1, j_idx] = p_hg * z_hg * i_coef[i - 1] * b_coef[i - 1]
        last_idx = i
        done_idx = np.abs(sum_terms[i - 2, j_idx] - sum_terms[i - 1, j_idx]) < 1e-4
        j_idx = j_idx[~done_idx]
        if len(j_idx) == 0:
            break
    if np.any(np.abs(sum_terms[N - 2, :] - sum_terms[N - 1, :]) > 1e-4):
        warnings.warn("Series not yet converged at N = 250.")
    # Use last computed row index for safe slicing
    res = np.sum(sum_terms[:last_idx, :], axis=0)
    res[res < 0] = 0.0
    return res if len(res) > 1 else float(res)


def dkimura(x: Union[float, np.ndarray], p: float, b: float) -> Union[float, np.ndarray]:
    """Density matching the R implementation."""
    x_arr = np.atleast_1d(x).astype(float)
    # Special cases
    if b == 0:
        # Bernoulli at p
        y = np.zeros_like(x_arr)
        y[x_arr == 0] = (1 - p)
        y[x_arr == 1] = p
        return y[0] if np.isscalar(x) else y
    
    if b == 1:
        y = np.zeros_like(x_arr)
        y[np.isclose(x_arr, p)] = 1.0
        return y[0] if np.isscalar(x) else y

    y = np.zeros_like(x_arr)
    y[np.isclose(x_arr, 0.0)] = _f0(p, b)
    y[np.isclose(x_arr, 1.0)] = _f1(p, b)

    mask = (x_arr > 0.0) & (x_arr < 1.0)
    z = x_arr[mask]
    lz = len(z)
    if lz == 0:
        return y[0] if np.isscalar(x) else y

    q = 1 - p
    N = 250
    sum_terms = np.zeros((N, lz), dtype=float)
    i_arr = np.arange(1, N + 1)
    i_coef = p * q * (i_arr * (i_arr + 1) * (2 * i_arr + 1))
    b_coef = b ** _comb(np.arange(2, N + 2), 2)

    p_hg = _hypgeo(1, p)
    z_hg = np.array([_hypgeo(1, zv) for zv in z])
    sum_terms[0, :] = p_hg * z_hg * i_coef[0] * b_coef[0]

    j_idx = np.arange(lz)
    # Track last computed term index to avoid referencing 'i' after the loop
    last_idx = 1
    for i in range(2, N + 1):
        z_hg = np.array([_hypgeo(i, z[j]) for j in j_idx])
        p_hg = _hypgeo(i, p)
        sum_terms[i - 1, j_idx] = p_hg * z_hg * i_coef[i - 1] * b_coef[i - 1]
        last_idx = i
        done_idx = np.abs(sum_terms[i - 2, j_idx] - sum_terms[i - 1, j_idx]) < 1e-4
        j_idx = j_idx[~done_idx]
        if len(j_idx) == 0:
            break
        
    if np.any(np.abs(sum_terms[N - 2, :] - sum_terms[N - 1, :]) > 1e-4):
        warnings.warn("Series not yet converged at N = 250.")

    # Use last computed index for slicing; if for some reason last_idx was not updated,
    # fall back to N
    if last_idx == 0:
        last_idx = N

    y[mask] = np.sum(sum_terms[:last_idx, :], axis=0)
    y[y < 0] = 0.0
    
    return y[0] if np.isscalar(x) else y


def _pkimura_full(p: float, b: float) -> np.ndarray:
    """Compute CDF on grid x=seq(0,1,1e-4) using trapezoidal rule (R's approach)."""
    x = np.linspace(0.0, 1.0, 10001)
    
    # Ensure pdf_x is always an ndarray (dkimura may return a scalar for scalar inputs)
    pdf_x = np.atleast_1d(dkimura(x, p, b))
    cdf_x = np.zeros_like(pdf_x)
    
    dx = 1e-4
    cdf_x[0] = pdf_x[0]
    cdf_x[1] = pdf_x[0] + 0.5 * dx * pdf_x[1]
    for i in range(2, len(x) - 1):
        cdf_x[i] = cdf_x[i - 1] + 0.5 * dx * (pdf_x[i - 1] + pdf_x[i + 1])
        
    cdf_x[-1] = 1.0
    return cdf_x


def pkimura(x: Union[float, np.ndarray], p: float, b: float) -> np.ndarray:
    """CDF matching the R implementation using interpolation from grid."""
    x_arr = np.atleast_1d(x).astype(float)
    # Special cases
    if b == 0:
        # CDF of Bernoulli
        return np.asarray((x_arr >= p).astype(float))
    if b == 1:
        return np.asarray((x_arr >= p).astype(float))

    y = np.zeros_like(x_arr)
    y[x_arr < 0] = 0.0
    y[np.isclose(x_arr, 0.0)] = _f0(p, b)
    y[x_arr >= 1.0] = 1.0

    mask = (x_arr > 0.0) & (x_arr < 1.0)
    z = x_arr[mask]
    if len(z) == 0:
        return y[0] if np.isscalar(x) else y

    # Ensure density arrays are ndarray so indexing (e.g. [ix]) works
    density_z = np.atleast_1d(dkimura(z, p, b))
    k = np.arange(1, int(np.max(z) * 1e4) + 2) * 1e-4
    # If k has length 0 (shouldn't happen for z>0) guard against empty arrays
    if k.size > 0:
        density_k = np.atleast_1d(dkimura(k, p, b))
        cdf_k = _pkimura_full(p, b)[:-1]  # drop last
    else:
        density_k = np.array([], dtype=float)
        cdf_k = np.array([], dtype=float)

    def get_cdf_for_idx(ix):
        xv = z[ix]
        ik = int(np.floor(xv * 1e4))
        if ik == 0:
            return _f0(p, b) + 0.5 * xv * density_z[ix]
        
        if (ik + 1) == int(1e4):
            d = 1 - xv
            return 1 - _f1(p, b) - 0.5 * d * density_z[ix]
        
        d = xv - k[ik - 1]
        return cdf_k[ik - 1] + 0.5 * d * (density_k[ik - 1] + density_z[ix])

    y[mask] = np.array([get_cdf_for_idx(i) for i in range(len(z))])
    return y[0] if np.isscalar(x) else y


def rkimura(n: int, p: float, b: float, random_state: Optional[Union[int, np.random.Generator]] = None) -> np.ndarray:
    """Sampling matching R's algorithm: select atoms or invert CDF by interpolation.

    random_state may be:
        - None: use a new Generator (np.random.default_rng())
        - int: interpreted as a seed for np.random.default_rng(seed)
        - np.random.Generator: used directly as the RNG
    """
    # Normalize RNG input to a Generator
    if isinstance(random_state, int):
        rng = np.random.default_rng(random_state)
    elif isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng()

    if b == 0:
        # Bernoulli draw
        return rng.binomial(1, p, size=n).astype(float)
    
    if b == 1:
        return np.full(n, p, dtype=float)

    x = np.linspace(0.0, 1.0, 10001)
    cdf_x = _pkimura_full(p, b)

    p0 = _f0(p, b)
    p1 = _f1(p, b)
    probs = [p0, 1 - p0 - p1, p1]
    choices = rng.choice([0.0, 0.5, 1.0], size=n, p=probs)
    
    idx = np.where(choices == 0.5)[0]
    m = len(idx)
    if m > 0:
        u = rng.uniform(low=cdf_x[0], high=1.0, size=m)
        # invert CDF by linear interpolation
        choices[idx] = np.interp(u, cdf_x, x)
        
    return choices


def test_kimura(data: np.ndarray, 
                n_simulations: int = 10, 
                random_state: Optional[int] = None) -> dict:
    """
    Test for selection pressure in heteroplasmy data.
    
    Args:
        data: Array of heteroplasmy values in [0, 1]
        n_simulations: Number of Monte Carlo simulations
        random_state: Random seed
        
    Returns:
        Dictionary of test results
    """
    results = KimuraTest.test(data, n_simulations, random_state)
    print_test_results(results)
    
    return results


if __name__ == "__main__":
    # Example usage matching the README
    print("=" * 60)
    print("Kimura Distribution Package - Example")
    print("=" * 60)
    
    # Load some heteroplasmy data
    h = np.array([0.06, 0.08, 0.27, 0.37, 0.40, 
                  0.45, 0.56, 0.61, 0.75, 0.79])
    
    # Carry out test for selection
    print("\nTest for selection pressure:")
    test_kimura(h, n_simulations=10, random_state=42)
    
    # Initialize Kimura parameters
    p = 0.6
    b = 0.95
    
    print(f"\n\nKimura distribution with p={p}, b={b}:")
    print("-" * 60)
    
    # Probability of allele loss
    print(f"\nProbability of allele loss (at x=0): {dkimura(0, p, b):.10f}")
    
    # Probability of fixing an allele
    print(f"Probability of fixation (at x=1): {dkimura(1, p, b):.10f}")
    
    # Kimura CDF at 0.1 intervals
    print(f"\nKimura({p}, {b}) CDF at 0.1 intervals:")
    
    intervals = np.arange(0, 1.1, 0.1)
    cdf_values = pkimura(intervals, p, b)
    for x, cdf_val in zip(intervals, cdf_values):
        print(f"  CDF({x:.1f}) = {cdf_val:.10f}")
    
    # Random number generation
    print(f"\n10 random variates from Kimura({p}, {b}):")
    # random_samples = rkimura(10, p, b, random_state=42)
    random_samples = rkimura(10, p, b)
    for i, sample in enumerate(random_samples, 1):
        print(f"  {i}: {sample:.7f}")
    
    print("\n" + "=" * 60)
