/**
 * Log normal distribution
 */

final class LogNormal(μ:Expression<Real>, σ2:Expression<Real>) < Distribution <Real> {
    /**
     * Mean of lognormal
     */
     μ:Expression<Real> <- μ;

     /**
      * Variance of lognormal
      */

      σ2:Expression<Real> <- σ2;

      function valueForward() -> Real {
        assert !delay?;
        return simulate_lognormal(μ, σ2);
      }

      function observeForward(x:Real) -> Real {
        assert !delay?;
        return exp(logpdf_gaussian(x, μ, σ2));
      }

      function graft(force:Boolean) {
      if delay? {
        delay!.prune();
      } else if force {
        delay <- DelayLogNormal(future, futureUpdate, μ, σ2);
      }
    }
}


/**
 * Creates a LogNormal Distribution.
 */

 function LogNormal(μ:Expression<Real>, σ2:Real) -> LogNormal {
    m:LogNormal(μ, σ2);
    return m;
}

 function LogNormal(μ:Expression<Real>, σ2:Real) -> LogNormal {
    return LogNormal(μ, Boxed(σ2));
}

 function LogNormal(μ:Real, σ2:Expression<Real>) -> LogNormal {
    return LogNormal(Boxed(μ), σ2);
}

 function LogNormal(μ:Real, σ2:Real) -> LogNormal {
    return LogNormal(Boxed(μ), Boxed(σ2));
}
