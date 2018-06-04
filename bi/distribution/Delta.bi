/**
 * Delta distribution, representing a distribution on a discrete space with
 * all probability mass at one location.
 */
class Delta(μ:Expression<Integer>) < Distribution<Integer> {
  /**
   * Location.
   */
  μ:Expression<Integer> <- μ;

  function graft() {
    if delay? {
      delay!.prune();
    } else {
      m:DelayDiscrete?;
      if (m <- μ.graftDiscrete())? {
        delay <- DelayDiscreteDelta(x, m!);
      } else {
        delay <- DelayDelta(x, μ);
      }
    }
  }
}

/**
 * Create delta distribution.
 *
 * - μ: Location.
 */
function Delta(μ:Expression<Integer>) -> Delta {
  m:Delta(μ);
  return m;
}

/**
 * Create delta distribution.
 *
 * - μ: Location.
 */
function Delta(μ:Integer) -> Delta {
  return Delta(Boxed(μ));
}