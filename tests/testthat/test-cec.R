test_that("all benchmark functions from CEC-2013 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 1:28,
    dim = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2013(func, rnorm(dim)), "double")
  })
})

test_that("all benchmark functions from CEC-2014 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 1:30,
    dim = c(10, 20, 30, 50, 100)
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2014(func, rnorm(dim)), "double")
  })
})

test_that("all benchmark functions from CEC-2015 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 1:15,
    dim = c(10, 30, 50, 100)
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2015(func, rnorm(dim)), "double")
  })
})

test_that("all benchmark functions from CEC-2017 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 1:30,
    dim = c(10, 30, 50, 100)
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2017(func, rnorm(dim)), "double")
  })
})

test_that("F1-F3 benchmark functions from CEC-2019 can be executed", {
  problem_dim_grid <- list(
    func = 1:3,
    dim = c(9, 16, 18)
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2019(func, rnorm(dim)), "double")
  })
})


test_that("F4-F10 benchmark functions from CEC-2019 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 4:10,
    dim = 10
  )
  purrr::pmap(problem_dim_grid, function(func, dim) {
    expect_type(cec2019(func, rnorm(dim)), "double")
  })
})

test_that("all benchmark functions from CEC-2021 can be executed", {
  problem_dim_grid <- expand.grid(
    func = 1:10,
    dim = c(10, 20),
    suite =
    c(
      "basic",
      "shift",
      "rot",
      "bias",
      "shift_rot",
      "bias_rot",
      "bias_shift",
      "bias_shift_rot"
    )
  )
  purrr::pmap(problem_dim_grid, function(func, dim, suite) {
    expect_type(cec2021(func, rnorm(dim), suite), "double")
  })
})


