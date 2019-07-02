/**
 * @file
 */
#pragma once

namespace libbirch {
/**
 * Atomic value.
 *
 * @tparam Value type.
 *
 * The implementation uses OpenMP atomics as opposed to std::atomic. The
 * advantage of this is ensured memory model consistency and the organic
 * disabling of atomics when OpenMP, and thus multithreading, is
 * disabled. The disadvantage is that OpenMP atomics do not support
 * compare-and-swap/compare-and-exchange, only swap/exchange, which requires
 * some clunkier client code, especially for read-write locks.
 *
 * Atomic provides the default constructor, copy and move constructors, copy
 * and move assignment operators, in order to be trivially copyable and so
 * a mappable type for the purposes of OpenMP. These constructors and
 * operators *do not* behave atomically, however.
 */
template<class T>
class Atomic {
public:
  /**
   * Default constructor.
   */
  Atomic() = default;

  /**
   * Constructor.
   *
   * @param value Initial value.
   */
  explicit Atomic(const T& value) {
    #pragma omp atomic write
    this->value = value;
  }

  /**
   * Load the value, atomically.
   */
  T load() const {
    T value;
    #pragma omp atomic read
    value = this->value;

    return value;
  }

  /**
   * Store the value, atomically.
   */
  void store(const T& value) {
    #pragma omp atomic write
    this->value = value;
  }

  /**
   * Exchange the value with another, atomically.
   *
   * @param value New value.
   *
   * @return Old value.
   */
  T exchange(const T& value) {
    T old;
    #pragma omp atomic capture
    {
      old = this->value;
      this->value = value;
    }
    return old;
  }

  template<class U>
  T operator+=(const U& value) {
    T result;
    #pragma omp atomic capture
    result = this->value += value;

    return result;
  }

  template<class U>
  T operator-=(const U& value) {
    T result;
    #pragma omp atomic capture
    result = this->value -= value;

    return result;
  }

  T operator++() {
    T value;
    #pragma omp atomic capture
    value = ++this->value;

    return value;
  }

  T operator++(int) {
    T value;
    #pragma omp atomic capture
    value = this->value++;

    return value;
  }

  T operator--() {
    T value;
    #pragma omp atomic capture
    value = --this->value;

    return value;
  }

  T operator--(int) {
    T value;
    #pragma omp atomic capture
    value = this->value--;

    return value;
  }

private:
  /**
   * Value.
   */
  T value;
};
}