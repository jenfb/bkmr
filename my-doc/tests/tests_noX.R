devtools::load_all()
library(ggplot2)

family <- "gaussian"
family <- "binomial"

## example where there is no X matrix ####

n <- 100
M <- 5
sigsq.true <- ifelse(family == "gaussian", 0.05, 1)
Z <- matrix(rnorm(n * M), n, M)
X <- cbind(3*cos(Z[, 1]) + 2*rnorm(n))
eps <- rnorm(n, sd = sqrt(sigsq.true))
h <- apply(Z, 1, function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3))
eps <- rnorm(n)
y <- drop(h + eps)
if (family == "binomial") {
  ystar <- y
  y <- ifelse(ystar > 0, 1, 0)
}

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, iter = 5000, varsel = TRUE, family = family)

fit0

summary(fit0)

TracePlot(fit = fit0, par = "beta")
ExtractPIPs(fit0)

pred.resp.univar <- PredictorResponseUnivar(fit = fit0)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

pred.resp.bivar <- PredictorResponseBivar(fit = fit0, 
                                          min.plot.dist = 1)

pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.df = pred.resp.bivar, 
                                                       Z = Z, qs = c(0.25, 0.5, 0.75))

ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")

risks.overall.approx <- OverallRiskSummaries(fit = fit0, 
                                             qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5)
risks.overall.approx

risks.overall.exact <- OverallRiskSummaries(fit = fit0, 
                                             qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")
risks.overall.exact

risks.singvar <- SingVarRiskSummaries(fit = fit0,
                                      qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75))
risks.singvar

risks.int <- SingVarIntSummaries(fit = fit0, 
                                 qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75))
risks.int
