#' Options for printing summary of model fit to the console
#'
#' Set options for what will be printed to the console when verbose = TRUE in the main kmbayes function
#'
#' @param verbose_freq After this percentage of iterations has been completed the summary of the model fit so far will be printed to the console 
#' @param verbose_show_ests TRUE or FALSE: flag indicating whether to print out summary statistics of all posterior samples obtained up until this point, for select parameters
#' @param verbose_digits Number of digits to be printed to the console 
#' 
#' @noRd
set_verbose_opts <- function(verbose_freq = NULL, verbose_show_ests = NULL, verbose_digits = NULL) {
  if (is.null(verbose_freq)) verbose_freq <- 10
  if (is.null(verbose_digits)) verbose_digits <- 5
  if (is.null(verbose_show_ests)) verbose_show_ests <- FALSE
  opts <- list(
    verbose_freq = verbose_freq,
    verbose_digits = verbose_digits,
    verbose_show_ests = verbose_show_ests
    )
  opts
}

print_diagnostics <- function(verbose, opts, curr_iter, tot_iter, chain, varsel, hier_varsel, ztest, Z, groups) {
  verbose_freq <- opts$verbose_freq
  verbose_digits <- opts$verbose_digits
  verbose_show_ests <- opts$verbose_show_ests
  s <- curr_iter
  nsamp <- tot_iter
  perc_iter_completed <- round(100*curr_iter/tot_iter, 1)
  
  all_iter <- 100*(1:nsamp)/nsamp
  sel_iter <- seq(verbose_freq, 100, by = verbose_freq)
  print_iter <- sapply(sel_iter, function(x) min(which(all_iter >= x)))
  
  elapsed_time <- difftime(Sys.time(), chain$time1)
  
  if (s %in% print_iter) {
    #if (verbose) message("------------------------------------------")
    if (verbose) cat("\n")
    message("Iteration: ", s, " (", perc_iter_completed, "% completed; ", round(elapsed_time, verbose_digits), " ", attr(elapsed_time, "units"), " elapsed)")
    
    if (verbose) {
      cat("Acceptance rates for Metropolis-Hastings algorithm:\n")
      accep_rates <- data.frame()
      ## lambda
      nm <- "lambda"
      rate <- colMeans(chain$acc.lambda[2:s, ,drop = FALSE])
      if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ## r_m
      if (!varsel) {
        nm <- "r"
        rate <- colMeans(chain$acc.r[2:s, , drop = FALSE])
        if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      } else {
        nm <- "r/delta (overall)"
        rate <- mean(chain$acc.rdelta[2:s])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        ##
        nm <- "r/delta  (move 1)"
        rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 1])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        ##
        nm <- "r/delta  (move 2)"
        rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 2])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        if (hier_varsel) {
          nm <- "r/delta  (move 3)"
          rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 3])
          accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        }
      }
      print(accep_rates)
      
      ## extra information
      if (verbose_show_ests) {
        cat("\nCurrent parameter estimates:\n")
        chain$varsel <- varsel
        class(chain) <- c("bkmrfit", class(chain))
        chain$Z <- Z
        if (hier_varsel) chain$groups <- groups
        ests <- ExtractEsts(chain, q = c(0.025, 0.975), sel = 2:s)
        #ests$h <- ests$h[c(1,2,nrow(ests$h)), ]
        summ <- with(ests, rbind(beta, sigsq.eps, r, lambda))
        summ <- data.frame(param = rownames(summ), round(summ, verbose_digits))
        rownames(summ) <- NULL
        print(summ)
        if (varsel) {
          cat("\nCurrent posterior inclusion probabilities:\n")
          pips <- ExtractPIPs(chain, sel = 2:s)
          pips[, -1] <- round(pips[, -1], verbose_digits)
          print(pips)
        }
      }
    }
  }
  
}