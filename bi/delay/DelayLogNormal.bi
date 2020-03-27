/**
 * Delayed handling of a log normal random variable
 */

final class DelayUniform(future:Real?, futureUpdate:Boolean, μ:Real, σ2:Real)
    < DelayValue<Real>(future, futureUpdate) {
    /**
     * μ: Mean of log normal
     */
     μ:Real <- μ;

     /**
      * Variance of log normal
      */
     λ:Real <- 1.0/σ2;

     function simulate() -> Real {
        return simulate_lognormal(μ, 1.0/λ);
     }

     function logpdf(x:Real) -> Real {
        return exp(logpdf_gaussian(x, μ, 1.0/λ);
     }

     function cdf(x:Real) -> Real? {
        return cdf_gaussian(log(x), μ, 1.0/λ);
     }

     function quantile(p:Real) -> Real? {
        return exp(quantile_gaussian(p, μ, 1.0/λ));
     }

     function write(buffer:Buffer) {
        prune(); //NOTE: What does this do?
        buffer.set("class": "LogNormal");
        buffer.set("μ": μ);
        buffer.set("σ2", 1.0/λ);
     }
}

function DelayLogNormal(future:Real?, futureUpdate:Boolean, μ:Real, σ2:Real) -> DelayGaussian {
    m:DelayLogNormal(future, futureUpdate, μ, σ2);
    return m;
}
