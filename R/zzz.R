.onAttach <- function(libname, pkgname) {
  packageStartupMessage("For guided examples, go to 'https://jenfb.github.io/bkmr/overview.html'")
}

release_questions <- function() {
  c(
    "Have you updated the vignette and posted to GitHub?"
  )
}