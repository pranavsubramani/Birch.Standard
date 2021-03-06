/*
 * Test multivariate Gaussian pdf evaluations.
 */
program test_pdf_multivariate_gaussian(D:Integer <- 4, N:Integer <- 10000,
    B:Integer <- 500, S:Integer <- 20) {
  μ:Real[D];
  Σ:Real[D,D];

  for i in 1..D {
    μ[i] <- simulate_uniform(-10.0, 10.0);
    for j in 1..D {
      Σ[i,j] <- simulate_uniform(-2.0, 2.0);
    }
  }
  Σ <- Σ*transpose(Σ);

  auto π <- Gaussian(μ, Σ);
  test_pdf(π, D, N, B, S);
}
