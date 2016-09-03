set.seed(12221)
library(sparcl)
library(compHclust)

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

sim2.test <- function(a, nsim=1100, nkeep=1000) {
  library(doMC)
  registerDoMC()

  genesa=matrix(0, nrow=50, ncol=12)
  genesa[1:20, 1:5] <- 6
  genesa[1:20, 6] <- a
  genesa[1:20, 7] <- -1*a
  genesa[1:20, 8:12] <- -6
  genesa[31:50, c(1,3,5,7,9,11)] <- 3
  genesa[31:50, c(2,4,6,8,10,12)] <- -3

  prim1 <- c(1,1,1,1,1,1,2,2,2,2,2,2)
  prim2<- c(2,2,2,2,2,2,1,1,1,1,1,1)
  sec1 <- c(1,2,1,2,1,2,1,2,1,2,1,2)
  sec2 <-  c(2,1,2,1,2,1,2,1,2,1,2,1)
  results <- foreach(i=1:nsim, .combine=rbind) %dopar% {
     x <- genesa + matrix(rnorm(50*12),ncol=12)

     tran=t(x)
     trans=scale(tran)

     out.psc <- rep(NA, 7)
     test <- try(permout <- KMeansSparseCluster.permute(trans, K = 2, nperm=5,
                                                        silent=TRUE),
                 silent=TRUE)
     if (!class(test)=="try-error") {
       sparse <- KMeansSparseCluster(trans, K = 2, wbounds=permout$bestw,
                                     silent=TRUE)
       signif <- sig.features(trans,sparse[[1]]$Cs)
       threshold = 1*(signif > (0.05/50))
       weights = threshold / sqrt(sum(threshold))
       test2 <- try(perm <- KMeansSparseCluster.permute2(trans, yval, K=2,
                                                         ws0=weights,
                                                         nperms=5,
                                                         silent=TRUE),
                    silent=TRUE)
       if (!class(test2)=="try-error") {
           sspcl.out <- KMeansSparseCluster2(trans, K=2, ws0=weights,
                                             wbounds=perm$bestw, silent=TRUE)
           out.psc[1] <-  (((sum(sparse[[1]]$Cs == prim1) == 12) |
                            (sum(sparse[[1]]$Cs == prim2) == 12)) &
                           ((sum(sspcl.out$Cs == sec1) == 12) |
                            (sum(sspcl.out$Cs == sec2) == 12)))
           out.psc[2] <-  (((sum(sparse[[1]]$Cs == sec1) == 12) |
                            (sum(sparse[[1]]$Cs == sec2) == 12)) &
                           ((sum(sspcl.out$Cs == prim1) == 12) |
                            (sum(sspcl.out$Cs == prim2) == 12)))
           out.psc[3] <-  (((sum(sparse[[1]]$Cs == prim1) == 12) |
                            (sum(sparse[[1]]$Cs == prim2) == 12)) &
                           ((sum(sspcl.out$Cs == sec1) != 12) &
                            (sum(sspcl.out$Cs == sec2) != 12)))
           out.psc[4] <-  (((sum(sparse[[1]]$Cs == sec1) == 12) |
                            (sum(sparse[[1]]$Cs == sec2) == 12)) &
                           ((sum(sspcl.out$Cs == prim1) != 12) &
                            (sum(sspcl.out$Cs == prim2) != 12)))
           out.psc[5] <- (((sum(sparse[[1]]$Cs == prim1) != 12) &
                          (sum(sparse[[1]]$Cs == prim2) != 12) &
                          (sum(sparse[[1]]$Cs == sec1) != 12) &
                          (sum(sparse[[1]]$Cs == sec2) != 12))&
                          ((sum(sspcl.out$Cs==prim1) == 12) |
                           (sum(sspcl.out$Cs==prim2) == 12)))
           out.psc[6] <- (((sum(sparse[[1]]$Cs == prim1) != 12) &
                          (sum(sparse[[1]]$Cs == prim2) != 12) &
                           (sum(sparse[[1]]$Cs == sec1) != 12) &
                          (sum(sparse[[1]]$Cs == sec2) != 12)) &
                          ((sum(sspcl.out$Cs==sec1) == 12) |
                           (sum(sspcl.out$Cs==sec2) == 12)))
           out.psc[7] <- (((sum(sparse[[1]]$Cs == prim1) != 12) &
                          (sum(sparse[[1]]$Cs == prim2) != 12) &
                           (sum(sparse[[1]]$Cs == sec1) != 12) &
                          (sum(sparse[[1]]$Cs == sec2) != 12)) &
                          ((sum(sspcl.out$Cs==prim1) != 12) &
                           (sum(sspcl.out$Cs==prim2) != 12) &
                           (sum(sspcl.out$Cs==sec1) != 12) &
                           (sum(sspcl.out$Cs==sec2) != 12)))
       }
     }

     x.hc <- hclust(as.dist(1-cor(x)))
     x.chc <- compHclust(x,x.hc)
     xp <- x.chc$x.prime
     xp.hc <- hclust(as.dist(1-cor(xp)))
     primaryN = cutree(x.hc, k=2)
     secondaryN = cutree(xp.hc, k=2)
     out.chc1 <- rep(NA, 7)
     out.chc1[1] <-  (((sum(primaryN == prim1) == 12) |
                       (sum(primaryN == prim2) == 12)) &
                      ((sum(secondaryN == sec1) == 12) |
                       (sum(secondaryN == sec2) == 12)))
     out.chc1[2] <-  (((sum(primaryN == sec1) == 12) |
                       (sum(primaryN == sec2) == 12)) &
                      ((sum(secondaryN == prim1) == 12) |
                       (sum(secondaryN == prim2) == 12)))
     out.chc1[3] <-  (((sum(primaryN == prim1) == 12) |
                       (sum(primaryN == prim2) == 12)) &
                      ((sum(secondaryN == sec1) != 12) &
                       (sum(secondaryN == sec2) != 12)))
     out.chc1[4] <-  (((sum(primaryN == sec1) == 12) |
                       (sum(primaryN == sec2) == 12)) &
                      ((sum(secondaryN == prim1) != 12) &
                       (sum(secondaryN == prim2) != 12)))
     out.chc1[5] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12))&
                     ((sum(secondaryN==prim1) == 12) |
                      (sum(secondaryN==prim2) == 12)))
     out.chc1[6] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12)) &
                     ((sum(secondaryN==sec1) == 12) |
                      (sum(secondaryN==sec2) == 12)))
     out.chc1[7] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12)) &
                     ((sum(secondaryN==prim1) != 12) &
                      (sum(secondaryN==prim2) != 12) &
                      (sum(secondaryN==sec1) != 12) &
                      (sum(secondaryN==sec2) != 12)))

     x.hc <- hclust(dist(t(x)))
     x.chc <- compHclust(x,x.hc)
     xp <- x.chc$x.prime
     xp.hc <- hclust(dist(t(xp)))
     primaryN = cutree(x.hc, k=2)
     secondaryN = cutree(xp.hc, k=2)
     out.chc2 <- rep(NA, 7)
     out.chc2[1] <-  (((sum(primaryN == prim1) == 12) |
                      (sum(primaryN == prim2) == 12)) &
                     ((sum(secondaryN == sec1) == 12) |
                      (sum(secondaryN == sec2) == 12)))
     out.chc2[2] <-  (((sum(primaryN == sec1) == 12) |
                      (sum(primaryN == sec2) == 12)) &
                     ((sum(secondaryN == prim1) == 12) |
                      (sum(secondaryN == prim2) == 12)))
     out.chc2[3] <-  (((sum(primaryN == prim1) == 12) |
                      (sum(primaryN == prim2) == 12)) &
                     ((sum(secondaryN == sec1) != 12) &
                      (sum(secondaryN == sec2) != 12)))
     out.chc2[4] <-  (((sum(primaryN == sec1) == 12) |
                      (sum(primaryN == sec2) == 12)) &
                     ((sum(secondaryN == prim1) != 12) &
                      (sum(secondaryN == prim2) != 12)))
     out.chc2[5] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12))&
                     ((sum(secondaryN==prim1) == 12) |
                      (sum(secondaryN==prim2) == 12)))
     out.chc2[6] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12)) &
                     ((sum(secondaryN==sec1) == 12) |
                      (sum(secondaryN==sec2) == 12)))
     out.chc2[7] <- (((sum(primaryN == prim1) != 12) &
                      (sum(primaryN == prim2) != 12) &
                      (sum(primaryN == sec1) != 12) &
                      (sum(primaryN == sec2) != 12)) &
                     ((sum(secondaryN==prim1) != 12) &
                      (sum(secondaryN==prim2) != 12) &
                      (sum(secondaryN==sec1) != 12) &
                      (sum(secondaryN==sec2) != 12)))

     hpermout <- HierarchicalSparseCluster.permute(t(x), nperms=5)
     sparsehc <- HierarchicalSparseCluster(dists=hpermout$dists,
                                           wbound=hpermout$bestw,
                                           method="complete",
                                           silent=TRUE)
     sparsehc2 <- HierarchicalSparseCluster(dists=hpermout$dists,
                                            wbound=hpermout$bestw,
                                            method="complete",
                                            uorth=sparsehc$u,
                                            silent=TRUE)
     primaryWN = cutree(sparsehc$hc, k=2)
     secondaryWN = cutree(sparsehc2$hc, k=2)
     out.wchc <- rep(NA, 7)
     out.wchc[1] <-  (((sum(primaryWN == prim1) == 12) |
                       (sum(primaryWN == prim2) == 12)) &
                      ((sum(secondaryWN == sec1) == 12) |
                       (sum(secondaryWN == sec2) == 12)))
     out.wchc[2] <-  (((sum(primaryWN == sec1) == 12) |
                       (sum(primaryWN == sec2) == 12)) &
                      ((sum(secondaryWN == prim1) == 12) |
                       (sum(secondaryWN == prim2) == 12)))
     out.wchc[3] <-  (((sum(primaryWN == prim1) == 12) |
                       (sum(primaryWN == prim2) == 12)) &
                      ((sum(secondaryWN == sec1) != 12) &
                       (sum(secondaryWN == sec2) != 12)))
     out.wchc[4] <-  (((sum(primaryWN == sec1) == 12) |
                       (sum(primaryWN == sec2) == 12)) &
                      ((sum(secondaryWN == prim1) != 12) &
                       (sum(secondaryWN == prim2) != 12)))
     out.wchc[5] <- ((sum(primaryWN == prim1) != 12) &
                     (sum(primaryWN == prim2) != 12) &
                     (sum(primaryWN == sec1) != 12) &
                     (sum(primaryWN == sec2) != 12))
     out.wchc[5] <- (((sum(primaryWN == prim1) != 12) &
                       (sum(primaryWN == prim2) != 12) &
                       (sum(primaryWN == sec1) != 12) &
                          (sum(primaryWN == sec2) != 12))&
                      ((sum(secondaryWN==prim1) == 12) |
                       (sum(secondaryWN==prim2) == 12)))
     out.wchc[6] <- (((sum(primaryWN == prim1) != 12) &
                       (sum(primaryWN == prim2) != 12) &
                       (sum(primaryWN == sec1) != 12) &
                       (sum(primaryWN == sec2) != 12)) &
                      ((sum(secondaryWN==sec1) == 12) |
                       (sum(secondaryWN==sec2) == 12)))
     out.wchc[7] <- (((sum(primaryWN == prim1) != 12) &
                       (sum(primaryWN == prim2) != 12) &
                       (sum(primaryWN == sec1) != 12) &
                       (sum(primaryWN == sec2) != 12)) &
                      ((sum(secondaryWN==prim1) != 12) &
                       (sum(secondaryWN==prim2) != 12) &
                       (sum(secondaryWN==sec1) != 12) &
                       (sum(secondaryWN==sec2) != 12)))
     c(out.psc, out.chc1, out.chc2, out.wchc)
  }
  results <- results[!is.na(results[,1]),]
  results <- results[1:nkeep,]
  return(colSums(results))
}

