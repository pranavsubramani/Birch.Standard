/**
 * Multivariate linear-normal-inverse-gamma-Gaussian random variable with
 * delayed sampling.
 */
class DelayMultivariateLinearNormalInverseGammaGaussian(x:Random<Real[_]>,
    A:Real[_,_], μ:DelayMultivariateNormalInverseGamma, c:Real[_]) <
    DelayValue<Real[_]>(x) {
  /**
   * Scale.
   */
  A:Real[_,_] <- A;

  /**
   * Mean.
   */
  μ:DelayMultivariateNormalInverseGamma <- μ;

  /**
   * Offset.
   */
  c:Real[_] <- c;

  function doSimulate() -> Real[_] {
    return simulate_multivariate_linear_normal_inverse_gamma_gaussian(A, μ.μ,
        c, μ.Λ, μ.σ2.α, μ.σ2.β);
  }
  
  function doObserve(x:Real[_]) -> Real {
    return observe_multivariate_linear_normal_inverse_gamma_gaussian(x, A,
        μ.μ, c, μ.Λ, μ.σ2.α, μ.σ2.β);
  }

  function doCondition(x:Real[_]) {
    (μ.μ, μ.Λ, μ.σ2.α, μ.σ2.β) <- update_multivariate_linear_normal_inverse_gamma_gaussian(
        x, A, μ.μ, c, μ.Λ, μ.σ2.α, μ.σ2.β);
  }
}

function DelayMultivariateLinearNormalInverseGammaGaussian(x:Random<Real[_]>,
    A:Real[_,_], μ:DelayMultivariateNormalInverseGamma, c:Real[_]) ->
    DelayMultivariateLinearNormalInverseGammaGaussian {
  m:DelayMultivariateLinearNormalInverseGammaGaussian(x, A, μ, c);
  return m;
}