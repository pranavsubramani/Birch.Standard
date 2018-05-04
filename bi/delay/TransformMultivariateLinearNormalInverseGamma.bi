/*
 * Linear transformation of a multivariate normal-inverse-gamma variate.
 */
class TransformMultivariateLinearNormalInverseGamma(A:Real[_,_],
    x:DelayMultivariateNormalInverseGamma, c:Real[_]) <
    TransformMultivariateLinear<Real>(A, c) {
  /**
   * Random variate.
   */
  x:DelayMultivariateNormalInverseGamma <- x;
}

function TransformMultivariateLinearNormalInverseGamma(A:Real[_,_],
    x:DelayMultivariateNormalInverseGamma, c:Real[_]) ->
    TransformMultivariateLinearNormalInverseGamma {
  m:TransformMultivariateLinearNormalInverseGamma(A, x, c);
  return m;    
}