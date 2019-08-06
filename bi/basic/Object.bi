/**
 * Root class of all other classes.
 */
class Object {
  /**
   * Read.
   */
  function read(buffer:Buffer?) {
    if (buffer?) {
      read(buffer!);
    }
  }

  /**
   * Read.
   */
  function read(buffer:Buffer) {
    //
  }
  
  /**
   * Write.
   */
  function write(buffer:Buffer?) {
    if (buffer?) {
      write(buffer!);
    }
  }

  /**
   * Write.
   */
  function write(buffer:Buffer) {
    //
  }
  
  /**
   * Touch the object. This is a null operation, but forces an update of the
   * pointer, which can be useful for debugging lazy deep clone issues.
   */
  function touch() {
    //
  }
}

/**
 * Make an object.
 *
 *   - name: Name of the class.
 *
 * Return: an optional with a value of the object if the make is successful,
 * or with no value otherwise.
 *
 * This is used to construct an object where the class is given in a string,
 * e.g. from user input. The class must not have constructor parameters.
 * Internally, it uses `dlsym()` to search the current process for a symbol
 * `make_name_` with C linkage, as generated by the Birch compiler for all
 * compatible types.
 */
function make(name:String) -> Object? {
  symbol:String <- "make_" + name + "_";
  cpp{{
  using make_t = bi::type::Object*();
  void* addr = dlsym(RTLD_DEFAULT, symbol.c_str());
  if (addr) {
    return reinterpret_cast<make_t*>(addr)();
  } else {
    return libbirch::Nil();
  }
  }}
}

/**
 * Identity comparison.
 */
operator (x:Object == y:Object) -> Boolean;

/**
 * Identity comparison.
 */
operator (x:Object != y:Object) -> Boolean;

/**
 * Identity conversion.
 */
function Object(o:Object) -> Object {
  return o;
}
