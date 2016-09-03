set.seed(12121)
library(sparcl)
load("../data/sspcl_funcs.RData")

#################################################
 KMeansSparseCluster.permute2 <- function (x, y, K = 2, nperms = 25,
                                          wbounds = NULL, silent = FALSE,
                                          nvals = 10,
                                          ws0=rep(1/sqrt(ncol(x)), ncol(x)),
                                          ...)
{
    if (is.null(wbounds))
        wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9),
            len = nvals))
    permx <- list()
    nnonzerows <- NULL
    nws <- sum(ws0!=0)
    for (i in 1:nperms) {
        permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
        for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
    }
    tots <- NULL
    out <- KMeansSparseCluster2(x, K, wbounds = wbounds, silent = silent,
    	   		        ws0=ws0, ...)
    for (i in 1:length(out)) {
        nnonzerows <- c(nnonzerows, sum(out[[i]]$ws != 0))
        bcss <- sparcl:::GetWCSS(x, out[[i]]$Cs)$bcss.perfeature
        tots <- c(tots, sum(out[[i]]$ws * bcss))
    }
    permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
    for (k in 1:nperms) {
        if (!silent)
            cat("Permutation ", k, "of ", nperms, fill = TRUE)
        perm.out <- KMeansSparseCluster2(permx[[k]], K, wbounds = wbounds,
            silent = silent, ws0=ws0, ...)
        for (i in 1:length(perm.out)) {
            perm.bcss <- sparcl:::GetWCSS(permx[[k]], perm.out[[i]]$Cs)$bcss.perfeature
            permtots[i, k] <- sum(perm.out[[i]]$ws * perm.bcss)
        }
    }
    gaps <- (log(tots) - apply(log(permtots), 1, mean))
    out <- list(tots = tots, permtots = permtots, nnonzerows = nnonzerows,
        gaps = gaps, sdgaps = apply(log(permtots), 1, sd), wbounds = wbounds,
        bestw = wbounds[which.max(gaps)])
    if (!silent)
        cat(fill = TRUE)
    class(out) <- "kmeanssparseperm"
    return(out)
}
#################################################



##################
#Create the data
##################
##################
qst <- read.csv("../data/QST.csv")
psych <- read.csv("../data/Psych.csv")
auto <- read.csv("../data/Autonomic.csv")
oppera.raw <- merge(qst, psych)
oppera.raw <- merge(oppera.raw, auto)
OPPERA <- data.frame(scale(oppera.raw)[,12:127])



##################
#Cluster the data
##################
##################
permout <- KMeansSparseCluster.permute(OPPERA, K = 2, nperm=5)
sparse <- KMeansSparseCluster(OPPERA, K = 2, wbounds=permout$bestw)[[1]]
weightsPRIMARY <- sparse$ws





##################
#Determine weights
##################
##################
sig.features <- function(x, Cs, perm=FALSE, B=1000)
{ library(sparcl)

  if (!perm) {
    n.Cs <- length(unique(Cs))
    junk <- sparcl:::GetWCSS(x, Cs)
    pvals <- 1-pf((junk$bcss.perfeature/(n.Cs-1))/
                  (junk$wcss.perfeature/(nrow(x)-n.Cs)), n.Cs-1,
                  nrow(x)-n.Cs)}
  else {
    x <- scale(x)
    x.bcss <- sparcl:::GetWCSS(x, Cs)$bcss.perfeature

   if (true) {
      nf <- ncol(x)
      xstar.bcss <- matrix(nrow=B, ncol=nf)
      for (i in 1:B) {
        x.star <- matrix(nrow=nrow(x), ncol=ncol(x))
        for (j in 1:ncol(x)) {
          x.star[,j] <- sample(x[,j])
        }
        xstar.bcss[i,] <- sparcl:::GetWCSS(x.star, Cs)$bcss.perfeature
      }
      pvals <- rep(NA, nf)
      for (i in 1:nf) {
        pvals[i] <- sum(x.bcss[i]<=xstar.bcss)/(B*nf)
      }    }  }
  return(pvals)
}






##Determine input weights for second round
signif <- sig.features(OPPERA,sparse$Cs)
threshold1 = 1*(signif > (0.05/116))
weights = threshold1 / sqrt(sum(threshold1))





##################
#Cluster
##################
##################
perm <- KMeansSparseCluster.permute2(OPPERA, yval, K=2, ws0=weights, nperms=5)
index <- which(perm$gaps == max(perm$gaps))
tuning <- perm$wbound[index]
sspcl.out <- KMeansSparseCluster2(OPPERA, K=2, ws0=weights, wbounds=tuning)
weightsSECONDARY <- sspcl.out$ws




##Determine input weights for third round
signif <- sig.features(OPPERA,sspcl.out$Cs)
threshold2 = 1*(signif > (0.05/116))

threshold = ((unname(threshold2) + unname(threshold1)) ==2) + 0

weights = threshold / sqrt(sum(threshold))




##################
#Cluster
##################
##################
perm <- KMeansSparseCluster.permute2(OPPERA, yval, K=2, ws0=weights, nperms=5)
index <- which(perm$gaps == max(perm$gaps))
tuning <- perm$wbound[index]
sspcl.out <- KMeansSparseCluster2(OPPERA, K=2, ws0=weights, wbounds=tuning)
weightsTERTIARY <- sspcl.out[[1]]$ws


setEPS()
postscript("Figure5.eps")
par(mfrow=c(3,1))
plot(weightsPRIMARY, ylab="Weight", xlab="Feature", main="Primary Cluster Weights")
plot(weightsSECONDARY, ylab="Weight", xlab="Feature", main="Secondary Cluster Weights")
plot(weightsTERTIARY, ylab="Weight", xlab="Feature", main="Tertiary Cluster Weights")
dev.off()
