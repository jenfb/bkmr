#' Options for printing summary of model fit to the console
#'
#' Set options for what will be printed to the console when verbose = TRUE in the main kmbayes function
#'
#' @param verbose_freq After this percentage of iterations has been completed the summary of the model fit so far will be printed to the console 
#'
#' @export
#'
set_verbose_opts <- function(verbose_freq) {
  opts <- list(verbose_freq = verbose_freq)
  opts
}

print_verbose <- function(opts, curr_iter, tot_iter, chain, varsel, hier_varsel, ztest) {
  verbose_freq <- opts$verbose_freq
  s <- curr_iter
  nsamp <- tot_iter
  perc_iter_completed <- round(100*curr_iter/tot_iter, 1)
  
  all_iter <- 100*(1:nsamp)/nsamp
  sel_iter <- seq(verbose_freq, 100, by = verbose_freq)
  print_iter <- sapply(sel_iter, function(x) min(which(all_iter >= x)))
  
  if (s %in% print_iter) {
    message("iter: ", s)
    cat(round(colMeans(chain$acc.lambda[1:s, ,drop = FALSE]), 4), "   lam accept rate\n")
    cat(round(colMeans(chain$acc.r[2:s, ]),4), "   r nosel accept rate\n")
    if (varsel) {
      cat(round(mean(chain$acc.rdelta[2:s]),4), "   rdelt accept rate\n")
      cat(round(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 1]),4), "   rdelt[move 1] accept rate\n")
      cat(round(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 2]),4), "   rdelt[move 2] accept rate\n")
      if (hier_varsel) cat(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 3]), "   rdelt[move 3] accept rate\n")
      cat(round(colMeans(chain$delta[1:s,ztest ,drop = FALSE]), 4), "   post incl probs\n")
      cat(round(colMeans(chain$r[2:s,], na.rm = TRUE),4), "   post mean of r\n")
    }
    print(difftime(Sys.time(), chain$time1))
  }
}