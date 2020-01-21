/**
 * Singly-linked tape. Maintains a current position along the tape, and
 * provides $O(1)$ operations at that position, but in the worst case an
 * $O(N)$ seek is required before operations at another position. Beyond its
 * typical uses, because Tape is a recursive and single-linked data
 * structure, it provides excellent sharing under Birch's lazy deep copy
 * mechanism.
 *
 * !!! attention
 *     See note under List for possible stack overflow issues on the
 *     destruction of large Tape objects. 
 */
final class Tape<Type> {
  forward:StackNode<Type>?;
  backward:StackNode<Type>?;
  forwardCount:Integer <- 0;
  backwardCount:Integer <- 0;

  /**
   * Number of elements.
   */
  function size() -> Integer {
    return forwardCount + backwardCount;
  }

  /**
   * Is this empty?
   */
  function empty() -> Boolean {
    return forwardCount + backwardCount == 0;
  }

  /**
   * Clear all elements.
   */
  function clear() {
    forward <- nil;
    backward <- nil;
    forwardCount <- 0;
    backwardCount <- 0;
  }

  /**
   * Get the first element.
   */
  function front() -> Type {
    rewind();
    return here();
  }

  /**
   * Get the last element.
   */
  function back() -> Type {
    fastForward();
    return before();
  }

  /**
   * Is there an element at the current position?
   */
  function hasHere() -> Boolean {
    return forward?;
  }

  /**
   * Is there an element before the current position?
   */
  function hasBefore() -> Boolean {
    return backward?;
  }

  /**
   * Get the element at the current position.
   */
  function here() -> Type {
    assert forward?;
    return forward!.x;
  }

  /**
   * Get the element before the current position.
   */
  function before() -> Type {
    assert backward?;
    return backward!.x;
  }

  /**
   * Insert a new element at the start.
   *
   * - x: Value.
   */
  function pushFront(x:Type) {
    rewind();
    pushHere(x);
  }

  /**
   * Insert a new element at the end.
   *
   * - x: Value.
   */
  function pushBack(x:Type) {
    fastForward();
    pushBefore(x);
  }

  /**
   * Remove the first element and return it.
   */
  function popFront() -> Type {
    rewind();
    return popHere();
  }

  /**
   * Remove the last element and return it.
   */
  function popBack() -> Type {
    fastForward();
    return popBefore();
  }

  /**
   * Insert a new element at the current position.
   *
   * - x: Value.
   */
  function pushHere(x:Type) {
    node:StackNode<Type>(x);
    node.next <- forward;
    forward <- node;
    forwardCount <- forwardCount + 1;
  }

  /**
   * Insert a new element before the current position.
   *
   * - x: Value.
   */
  function pushBefore(x:Type) {
    node:StackNode<Type>(x);
    node.next <- backward;
    backward <- node;
    backwardCount <- backwardCount + 1;
  }

  /**
   * Remove the element at the current position and return it.
   */
  function popHere() -> Type {
    auto node <- forward!;
    forward <- forward!.next;
    forwardCount <- forwardCount - 1;
    return node.x;
  }

  /**
   * Remove the element before the current position and return it.
   */
  function popBefore() -> Type {
    auto node <- backward!;
    backward <- backward!.next;
    backwardCount <- backwardCount - 1;
    return node.x;
  }
  
  /**
   * Move the current position forward one.
   */
  function next() {
    auto node <- forward!;
    forward <- node.next;
    node.next <- backward;
    backward <- node;
    forwardCount <- forwardCount - 1;
    backwardCount <- backwardCount + 1;
  }
  
  /**
   * Move the current position back one.
   */
  function previous() {
    auto node <- backward!;
    backward <- node.next;
    node.next <- forward;
    forward <- node;
    forwardCount <- forwardCount + 1;
    backwardCount <- backwardCount - 1;
  }
  
  /**
   * Rewind to the head to the first element.
   */
  function rewind() {
    while hasBefore() {
      previous();
    }
  }

  /**
   * Fast-forward the head to one past the last element.
   */
  function fastForward() {
    while hasHere() {
      next();
    }
  }

  /**
   * Forward iteration.
   *
   * Return: a fiber object that yields each element in order.
   */
  fiber walk() -> Type {
    rewind();
    while hasHere() {
      yield here();
      next();
    }
  }

  function read(buffer:Buffer) {
    auto f <- buffer.walk();
    while f? {
      /* tricky, but works for both value and class types */
      auto x <- make<Type>();
      auto y <- f!.get(x);
      if y? {
        x <- Type?(y);  // cast needed for y:Object?
        pushBack(x!);
      }
    }
    rewind();
  }

  function write(buffer:Buffer) {
    buffer.setArray();
    auto f <- walk();
    while f? {
      buffer.push().set(f!);
    }
  }
}