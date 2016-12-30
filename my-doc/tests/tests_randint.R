devtools::load_all()
library(ggplot2)

## example where there is no X matrix ####

load("H:/Research/2. BKMR/BKMR R package/workspace_vignette.RData")

n <- nrow(X)
id <- rep(1:(n/2), each = 2)
h <- apply(Z, 1, function(z, ind = 1) 4*plogis(z[ind[1]], 0, 0.3))
u <- rep(rnorm(n/2), each = 2)
eps <- rnorm(n, 0, 0.2)
y <- drop(h + u + eps)

set.seed(111)
fit0 <- kmbayes(y = y, Z = Z, iter = 5000, id = id, varsel = TRUE, control.params = list(verbose_show_ests = TRUE))

fit0

## tox application from BKMR paper ####

## load & prep data
## from Brent's HEI report

DIR <- "H:/Research/Completed Projects/2014 Bayesian kernel machine regression (Biostatistics)/Code/tox application/"

meanvals.dat <- read.table(paste0(DIR, "data/bp-xrf.csv"),sep=",",header=T)
exposure.dat <- read.table(paste0(DIR, "data/BC.csv"),sep=",",header=T)

caps.dat <- meanvals.dat[meanvals.dat$Exposure=="CAPs",]
sham.dat <- meanvals.dat[meanvals.dat$Exposure=="Sham",]

bccaps.dat <- merge(exposure.dat,caps.dat,by="DATE") 
bcsham.dat <- sham.dat
bcsham.dat$BC <- 0

bcall.dat <- rbind(bccaps.dat,bcsham.dat)
bcall.dat <- bcall.dat[order(bcall.dat$seq),]
bcall.dat <- bcall.dat[bcall.dat$seq!=112,]
bcall.dat <- bcall.dat[bcall.dat$seq!=11,]

bcall.dat$Dog <- as.numeric(bcall.dat$Dog) 
bcall.dat$exp <- rep(0,length(bcall.dat$Exposure))
bcall.dat$exp[bcall.dat$Exposure=="CAPs"] <- 1
bcall.dat$stat <- rep(0,length(bcall.dat$Status2))
bcall.dat$stat[bcall.dat$Status2=="Post-Occlusio"] <- 1
bcall.dat$stat2 <- rep(0,length(bcall.dat$Status2))
bcall.dat$stat2[bcall.dat$Status2=="Prazosin"] <- 1


n.times <- tapply(rep(1,length(bcall.dat$Dog)),bcall.dat$Dog,sum)
time.var <- NULL
for (i in 1:length(n.times))
{
  time.var <- c(time.var,1:n.times[i])	
}
bcall.dat$time <- time.var

#bcall.dat <- bcall.dat[bcall.dat$stat2==0,]

dat<-NULL
#corresponds to new simulation study (except doesn't include bc because variable wasn't in bcall.dat and added Mn because that was the observed effect in the original study)
varnames <- c("Al","Si","Ti","Ca","K","Cu","Mn", "Ni","V","Zn", "S", "Cl", "BC")
groups <- c(1,1,1,1,1,1,1, 2,2,2, 3, 4, 5)
dat$Z <- as.matrix(bcall.dat[,varnames])
mean.z <- apply(dat$Z,2,mean)
dat$Z <- sweep(dat$Z,2,mean.z)
sd.z <- apply(dat$Z,2,sd)
dat$Z <- sweep(dat$Z,2,sd.z,FUN="/")

dat$X <- cbind(bcall.dat$stat,bcall.dat$stat2,bcall.dat$exp)
dat$y <- matrix(bcall.dat$mmrate,length(bcall.dat$mmrate),1)
dat$id <- bcall.dat$Dog

## take out outliers
inds <- seq(1,dim(dat$X)[1])
noout.inds <- inds[dat$Z[,1] < 6 & dat$Z[,2] < 6 & dat$Z[,3] < 6 & dat$Z[,4] < 6 & dat$Z[,5] < 6 & dat$Z[,6] < 6 & dat$Z[,7] < 6 & dat$Z[,8] < 6 & dat$Z[,9] < 6 & dat$Z[,10] < 6 & dat$Z[,11] < 6 & dat$Z[,12] < 6 & dat$Z[,13] < 6]

dat$y <- as.matrix(dat$y[noout.inds])
dat$X <- as.matrix(dat$X[noout.inds,])
dat$Z <- dat$Z[noout.inds,]
# dat$Z <- scale(dat$Z)
dat$id <- dat$id[noout.inds]

## set up and fit

set.seed(111)
fitkm <- kmbayes(y = dat$y, X = dat$X, Z = dat$Z, id = dat$id, iter = 1000, verbose = TRUE, varsel = TRUE)