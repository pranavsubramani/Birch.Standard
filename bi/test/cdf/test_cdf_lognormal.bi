/*
 * Test lognormal cdf evaluations.
 */
program test_cdf_lognormal(N:Integer <- 10000) {
  auto μ <- simulate_uniform(-10.0, 10.0);
  auto σ2 <- simulate_uniform(0.0, 10.0);
  auto q <- LogNormal(μ, σ2);
  test_cdf(q, N);
}
