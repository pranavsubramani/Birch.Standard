/**
 * Particle filter.
 */
class ParticleFilter < ForwardSampler {  
  /**
   * Particles.
   */
  x:ForwardModel[_];
  
  /**
   * Log-weights.
   */
  w:Real[_];

  /**
   * Ancestor indices.
   */
  a:Integer[_];
  
  /**
   * Chosen path at the end of the filter.
   */
  x':ForwardModel?;
  
  /**
   * Number of particles.
   */
  N:Integer <- 1;
  
  /**
   * Number of steps.
   */
  T:Integer <- 1;
  
  /**
   * Threshold for resampling. Resampling is performed whenever the
   * effective sample size, as a proportion of `N`, drops below this
   * threshold.
   */
  trigger:Real <- 0.7;
  
  /**
   * For each checkpoint, the logarithm of the normalizing constant estimate
   * so far.
   */
  Z:List<Real>;
  
  /**
   * For each checkpoint, the effective sample size (ESS).
   */
  ess:List<Real>;
  
  /**
   * At each checkpoint, how much memory is in use?
   */
  memory:List<Integer>;
  
  /**
   * At each checkpoint, what is the elapsed wallclock time?
   */
  elapsed:List<Real>; 

  function sample() -> (Model, Real) {
    initialize();
    if verbose {
      stderr.print("steps:");
    }
    start();
    reduce();
    for auto t in 1..T {
      if verbose {
        stderr.print(" " + t);
      }
      if ess.back() < trigger*N {
        resample();
        copy();
      }
      step();
      reduce();
    }
    finish();
    if verbose {
      stderr.print(", log weight: " + sum(Z.walk()) + "\n");
    }
    finalize();
    return (x'!, sum(Z.walk()));
  }

  /**
   * Initialize.
   */
  function initialize() {
    Z.clear();
    ess.clear();
    memory.clear();
    elapsed.clear();

    w <- vector(0.0, N);
    a <- iota(1, N);
    x1:Vector<ForwardModel>;
    x1.enlarge(N, archetype!);
    x <- x1.toArray();
    parallel for auto n in 1..N {
      x[n] <- clone<ForwardModel>(x[n]);
    }
    tic();
  }
   
  /**
   * Start particles.
   */
  function start() {
    parallel for auto n in 1..N {
      w[n] <- w[n] + x[n].start();
    }
  }
    
  /**
   * Step particles.
   */
  function step() {
    parallel for auto n in 1..N {
      w[n] <- w[n] + x[n].step();
    }
  }

  /**
   * Compute summary statistics.
   */
  function reduce() {
    m:Real <- max(w);
    W:Real <- 0.0;
    W2:Real <- 0.0;
    
    for auto n in 1..length(w) {
      auto v <- exp(w[n] - m);
      W <- W + v;
      W2 <- W2 + v*v;
    }
    auto V <- log(W) + m - log(N);
    w <- w - V;  // normalize weights to sum to N
   
    /* effective sample size */
    ess.pushBack(W*W/W2);
    if !(ess.back() > 0.0) {  // > 0.0 as may be nan
      error("particle filter degenerated.");
    }
  
    /* normalizing constant estimate */
    Z.pushBack(V);
    memory.pushBack(memoryUse());
    elapsed.pushBack(toc());
  }

  /**
   * Resample particles.
   */
  function resample() {
    a <- permute_ancestors(ancestors(w));
  }
  
  /**
   * Copy around particles after resampling.
   */
  function copy() {
    auto x0 <- x;
    parallel for auto n in 1..N {
      //if a[n] != n {
        x[n] <- clone<ForwardModel>(x0[a[n]]);
      //}
      w[n] <- 0.0;
    }
  }

  /**
   * Finish.
   */
  function finish() {
    auto b <- ancestor(w);
    if b > 0 {
      x' <- x[b];
    } else {
      error("particle filter degenerated.");
    }
  }
  
  /**
   * Finalize.
   */
  function finalize() {
    //
  }

  function setArchetype(a:Model) {
    super.setArchetype(a);
    T <- archetype!.size();
  }

  function read(buffer:Buffer) {
    super.read(buffer);
    N <-? buffer.get("nparticles", N);
    trigger <-? buffer.get("trigger", trigger);
  }

  function write(buffer:Buffer) {
    super.write(buffer);
    buffer.set("nparticles", N);
    buffer.set("trigger", trigger);
    buffer.set("levidence", Z);
    buffer.set("ess", ess);
    buffer.set("elapsed", elapsed);
    buffer.set("memory", memory);
  }
}