/*
 * Lazy multivariate multiplication.
 */
final class MultivariateMultiply<Left,Right,Value>(left:Expression<Left>,
    right:Expression<Right>) < Expression<Value> {  
  /**
   * Left operand.
   */
  left:Expression<Left> <- left;
  
  /**
   * Right operand.
   */
  right:Expression<Right> <- right;
  
  function value() -> Value {
    return left.value()*right.value();
  }

  function graftMultivariateLinearGaussian() ->
      TransformMultivariateLinearGaussian? {
    y:TransformMultivariateLinearGaussian?;
    z:DelayMultivariateGaussian?;
    
    if (y <- right.graftMultivariateLinearGaussian())? {
      y!.leftMultiply(left.value());
    } else if (z <- right.graftMultivariateGaussian())? {
      y <- TransformMultivariateLinearGaussian(left.value(), z!);
    }
    return y;
  }
  
  function graftMultivariateLinearNormalInverseGamma(σ2:Expression<Real>) ->
      TransformMultivariateLinearNormalInverseGamma? {
    y:TransformMultivariateLinearNormalInverseGamma?;
    z:DelayMultivariateNormalInverseGamma?;

    if (y <- right.graftMultivariateLinearNormalInverseGamma(σ2))? {
      y!.leftMultiply(left.value());
    } else if (z <- right.graftMultivariateNormalInverseGamma(σ2))? {
      y <- TransformMultivariateLinearNormalInverseGamma(left.value(), z!);
    }
    return y;
  }
  
  function graftMultivariateScaledInverseGamma(σ2:Expression<Real>) ->
      TransformMultivariateScaledInverseGamma? {
    y:TransformMultivariateScaledInverseGamma?;
    z:DelayInverseGamma?;
    
    if (y <- right.graftMultivariateScaledInverseGamma(σ2))? {
      y!.leftMultiply(left.value());
    } else if Object(right) == σ2 && (z <- right.graftInverseGamma())? {
      y <- TransformMultivariateScaledInverseGamma(left.value(), z!);
    }
    return y;
  }
}

operator (left:Expression<Real[_,_]>*right:Expression<Real[_]>) ->
    MultivariateMultiply<Real[_,_],Real[_],Real[_]> {
  m:MultivariateMultiply<Real[_,_],Real[_],Real[_]>(left, right);
  return m;
}

operator (left:Real[_,_]*right:Expression<Real[_]>) ->
    MultivariateMultiply<Real[_,_],Real[_],Real[_]> {
  return Boxed(left)*right;
}

operator (left:Expression<Real[_,_]>*right:Real[_]) ->
    MultivariateMultiply<Real[_,_],Real[_],Real[_]> {
  return left*Boxed(right);
}

operator (left:Expression<Real[_,_]>*right:Expression<Real>) ->
    MultivariateMultiply<Real[_,_],Real,Real[_,_]> {
  m:MultivariateMultiply<Real[_,_],Real,Real[_,_]>(left, right);
  return m;
}

operator (left:Real[_,_]*right:Expression<Real>) ->
    MultivariateMultiply<Real[_,_],Real,Real[_,_]> {
  return Boxed(left)*right;
}

operator (left:Expression<Real[_,_]>*right:Real) ->
    MultivariateMultiply<Real[_,_],Real,Real[_,_]> {
  return left*Boxed(right);
}
