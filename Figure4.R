set.seed(1234)
library(sparcl)

load("../data/sspcl_funcs.RData")

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

genesa=matrix(0, nrow=50, ncol=100)
genesa[1:20, 1:50] <- genesa[1:20, 1:50]+6
genesa[1:20, 51:100] <- genesa[1:20, 51:100]-6
genesa[11:30, 2*(1:50)-1] <- genesa[11:30, 2*(1:50)-1]+3
genesa[11:30, 2*(1:50)] <- genesa[11:30, 2*(1:50)]-3

x <- genesa + matrix(rnorm(50*100),ncol=100)
tran=t(x)
trans=scale(tran)

permout <- KMeansSparseCluster.permute(trans, K = 2, nperm=5,
                                       silent=TRUE)
sparse <- KMeansSparseCluster(trans, K = 2, wbounds=permout$bestw,
                              silent=TRUE)

signif <- sig.features(trans,sparse[[1]]$Cs)
threshold = 1*(signif > (0.05/50))
weights = threshold / sqrt(sum(threshold))

perm <- KMeansSparseCluster.permute2(trans, yval, K=2, ws0=weights,
                                     nperms=5, silent=TRUE)

sspcl.out <- KMeansSparseCluster2(trans, K=2, ws0=weights,
                                  wbounds=perm$bestw, silent=TRUE)

setEPS()
postscript("Figure4.eps")
par(mfrow=c(3,1))
plot(sparse[[1]]$ws, ylab="Weight", xlab="Feature",
     main="Feature Weights (primary)")
plot(weights, ylab="Weight", xlab="Feature",
     main="Initial Feature Weights (Secondary)")
plot(sspcl.out$ws, ylab="Weight", xlab="Feature",
     main="Final Feature Weights (Secondary)")
dev.off()
