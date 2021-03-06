/*
 * Test deep clone of an object, where an object is cloned, the original
 * modified, and then that original accessed via an alias pointer.
 */
program test_deep_clone_alias() {
  /* create a simple list */
  x:List<Integer>;
  x.pushBack(1);
  x.pushBack(2);
  
  /* alias it */
  auto z <- x;
  
  /* clone it */
  auto y <- clone<List<Integer>>(x);

  /* modify the original */
  x.set(1, 3);
  x.set(2, 4);
  
  /* check that the alias has redirected */
  if (z.get(1) != 3 || z.get(2) != 4) {
    exit(1);
  }
}
