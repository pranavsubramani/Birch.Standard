/*
 * Test linear matrix normal-inverse-gamma Gaussian pdf evaluations.
 */
program test_pdf_linear_matrix_normal_inverse_gamma_matrix_gaussian(
    N:Integer <- 20000, B:Integer <- 100, S:Integer <- 20) {
  m:TestLinearMatrixNormalInverseGammaMatrixGaussian;
  m.play();
  test_pdf(m.marginal(), 5, 2, N, B, S);
}