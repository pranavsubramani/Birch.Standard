/**
 * The gamma function.
 */
function gamma(x:Real64) -> Real64 {
  cpp {{
  return ::tgamma(x);
  }}
}

/**
 * The gamma function.
 */
function gamma(x:Real32) -> Real32 {
  cpp {{
  return ::tgammaf(x);
  }}
}

/**
 * Logarithm of the gamma function.
 */
function lgamma(x:Real64) -> Real64 {
  cpp {{
  return ::lgamma(x);
  }}
}

/**
 * Logarithm of the gamma function.
 */
function lgamma(x:Real32) -> Real32 {
  cpp {{
  return ::lgammaf(x);
  }}
}

/**
 * The multivariate gamma function.
 */
function gamma(x:Real64, p:Integer) -> Real64 {
  assert p > 0;
  auto y <- 0.25*(p*(p - 1))*log(π);
  for i in 1..p {
    y <- y*gamma(x + 0.5*(1 - i));
  }
  return y;
}

/**
 * The multivariate gamma function.
 */
function gamma(x:Real32, p:Integer) -> Real32 {
  assert p > 0;
  auto y <- Real32(0.25)*Real32(p*(p - 1))*log(Real32(π));
  for i in 1..p {
    y <- y*gamma(x + Real32(0.5)*Real32(1 - i));
  }
  return y;
}

/**
 * Logarithm of the multivariate gamma function.
 */
function lgamma(x:Real64, p:Integer) -> Real64 {
  assert p > 0;
  auto y <- 0.25*(p*(p - 1))*log(π);
  for i in 1..p {
    y <- y + lgamma(x + 0.5*(1 - i));
  }
  return y;
}

/**
 * Logarithm of the multivariate gamma function.
 */
function lgamma(x:Real32, p:Integer) -> Real32 {
  assert p > 0;
  auto y <- Real32(0.25)*Real32(p*(p - 1))*log(Real32(π));
  for i in 1..p {
    y <- y + lgamma(x + Real32(0.5)*Real32(1 - i));
  }
  return y;
}

/**
 * The beta function.
 */
function beta(x:Real64, y:Real64) -> Real64 {
  cpp {{
  return boost::math::beta(x, y);
  }}
}

/**
 * The beta function.
 */
function beta(x:Real32, y:Real32) -> Real32 {
  cpp {{
  return boost::math::beta(x, y);
  }}
}

/**
 * The incomplete beta function.
 */
function ibeta(a:Real64, b:Real64, x:Real64) -> Real64 {
  if x < 0.0 || x > 1.0 {
    return inf;
  }
  // whenever x < (a + 1) / (a + b + 2) -- CF converges
  if x > (a + 1.0)/(a + b + 2.0) {
    return 1.0 - ibeta(b, a, 1.0 - x);
  }

  auto lbeta_ab <- lgamma(a) + lgamma(b) - lgamma(a + b);
  auto front <- exp(log(x)*a + log(1.0 - x)*b - lbeta_ab) / a;

  // setup for Lentz's algorithm
  auto f <- 1.0;
  auto c <- 1.0;
  auto d <- 0.0;

  for auto i in 0..200 {
    auto m = i/2;

    auto numerator <- 1.0;

    if i == 0 {
      numerator <- 1.0;
    } else if mod(i, 2.0) {

    }

  }

}

/**
 * The incomplete beta function.
 */
function ibeta(a:Real32, b:Real32, x:Real32) -> Real32 {
  cpp {{
    return boost::math::ibeta(a, b, x);
  }}
}

/**
 * Logarithm of the beta function.
 */
function lbeta(x:Real64, y:Real64) -> Real64 {
  return lgamma(x) + lgamma(y) - lgamma(x + y);
}

/**
 * Logarithm of the beta function.
 */
function lbeta(x:Real32, y:Real32) -> Real32 {
  return lgamma(x) + lgamma(y) - lgamma(x + y);
}

/**
 * The binomial coefficient.
 */
function choose(x:Real64, y:Real64) -> Real64 {
  assert 0.0 <= x;
  assert 0.0 <= y;
  assert x >= y;

  if (y == 0.0) {
    return 1.0;
  } else {
    // see Boost binomial_coefficient function for this implementation
    return 1.0/(y*beta(y, x - y + 1.0));
  }
}

/**
 * The binomial coefficient.
 */
function choose(x:Real32, y:Real32) -> Real32 {
  assert Real32(0.0) <= x;
  assert Real32(0.0) <= y;
  assert x >= y;

  if (y == Real32(0.0)) {
    return Real32(1.0);
  } else {
    // see Boost binomial_coefficient function for this implementation
    return Real32(1.0)/(y*beta(y, x - y + Real32(1.0)));
  }
}

/**
 * Logarithm of the binomial coefficient.
 */
function lchoose(x:Real64, y:Real64) -> Real64 {
  assert 0.0 <= x;
  assert 0.0 <= y;
  assert x >= y;

  if (y == 0.0) {
    return 0.0;
  } else {
    // see Boost binomial_coefficient function for this implementation
    return -log(y) - lbeta(y, x - y + 1.0);
  }
}

/**
 * Logarithm of the binomial coefficient.
 */
function lchoose(x:Real32, y:Real32) -> Real32 {
  assert Real32(0.0) <= x;
  assert Real32(0.0) <= y;
  assert x >= y;

  if (y == Real32(0.0)) {
    return log(Real32(1.0));
  } else {
    // see Boost binomial_coefficient function for this implementation
    return -log(y) - lbeta(y, x - y + Real32(1.0));
  }
}
