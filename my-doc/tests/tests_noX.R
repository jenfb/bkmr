devtools::load_all()
library(ggplot2)

## example where there is no X matrix ####

load("H:/Research/2. BKMR/BKMR R package/workspace_vignette.RData")

n <- nrow(X)
h <- apply(Z, 1, function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3))
eps <- rnorm(n)
y <- drop(h + eps)

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, iter = 10000, verbose = TRUE, varsel = TRUE)

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
  facet_grid(variable1 ~ variable2) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")

risks.overall.approx <- OverallRiskSummaries(fit = fit0, 
                                             qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5)
risks.overall.approx

risks.singvar <- SingVarRiskSummaries(fit = fit0,
                                      qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75))
risks.singvar

risks.int <- SingVarIntSummaries(fit = fit0, 
                                 qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75))
risks.int
