/*
 * Test multivariate normal-inverse-gamma Gaussian pdf evaluations.
 */
program test_pdf_multivariate_normal_inverse_gamma_multivariate_gaussian(D:Integer <- 4,
    N:Integer <- 20000, B:Integer <- 100, S:Integer <- 20) {
  m:TestMultivariateNormalInverseGammaMultivariateGaussian;
  m.play();
  test_pdf(m.marginal(), 5, N, B, S);
}
