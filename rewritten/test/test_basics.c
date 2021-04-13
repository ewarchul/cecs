#define UNIT_TESTING 1

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>

#include <cmocka.h>

static void null_test_success(void **state) { (void)state; }

static void first_cmocka_test(void **state) {
  (void) state;
  assert_int_equal(5, 5);
}

int main() {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(null_test_success),
      cmocka_unit_test(first_cmocka_test),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
