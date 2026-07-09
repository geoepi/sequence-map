package_root <- normalizePath(
  testthat::test_path("..", ".."),
  winslash = "/",
  mustWork = TRUE
)

r_dir <- file.path(package_root, "R")
for (file in list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)) {
  sys.source(file, envir = globalenv())
}
