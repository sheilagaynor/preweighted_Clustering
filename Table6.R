set.seed(12121)
library(sparcl)
library(survival)
library(compHclust)
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
#Preweighted Sparse Clustering
##################
##################
perm <- KMeansSparseCluster.permute2(OPPERA, yval, K=2, ws0=weights, nperms=5)
sspcl.out <- KMeansSparseCluster2(OPPERA, K=2, ws0=weights, wbounds=perm$bestw)

oppera.raw$psclust1 <- sparse$Cs
oppera.raw$psclust2 <- sspcl.out$Cs

##################
#Complementary Hierarchical Clustering
##################
##################

op.hc <- hclust(as.dist(1-cor(t(OPPERA))))
op.comp <- compHclust(t(OPPERA), op.hc)
op.chc <- hclust(as.dist(1-cor(op.comp$x.prime)))

oppera.raw$chclust1 <- cutree(op.hc, k=2)
oppera.raw$chclust2 <- cutree(op.chc, k=2)

##################
#Complementary Sparse Hierarchical Clustering
##################
##################

#hpermout <- HierarchicalSparseCluster.permute(as.matrix(OPPERA), nperms=5)
#sparsehc <- HierarchicalSparseCluster(dists=hpermout$dists,
#                                      wbound=hpermout$bestw,
#                                      method="complete",
#                                      silent=TRUE)
#sparsehc2 <- HierarchicalSparseCluster(dists=hpermout$dists,
#                                       wbound=hpermout$bestw,
#                                       method="complete",
#                                       uorth=sparsehc$u,
#                                       silent=TRUE)

#oppera.raw$cshclust1 <- cutree(sparsehc$hc, k=2)
#oppera.raw$cshclust2 <- cutree(sparsehc2$hc, k=2)

survdiff(Surv(followupyears,incidentcase)~psclust1,data=oppera.raw)
psclust1cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(psclust1),data=oppera.raw))

survdiff(Surv(followupyears,incidentcase)~psclust2,data=oppera.raw)
psclust2cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(psclust2),data=oppera.raw))

survdiff(Surv(followupyears,incidentcase)~chclust1,data=oppera.raw)
chclust1cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(chclust1),data=oppera.raw))

survdiff(Surv(followupyears,incidentcase)~chclust2,data=oppera.raw)
chclust2cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(chclust2),data=oppera.raw))

#survdiff(Surv(followupyears,incidentcase)~cshclust1,data=oppera.raw)
#cshclust1cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(cshclust1),data=oppera.raw))

#survdiff(Surv(followupyears,incidentcase)~cshclust2,data=oppera.raw)
#cshclust2cox <- summary(coxph(Surv(followupyears,incidentcase)~factor(cshclust2),data=oppera.raw))

out <- matrix(nrow=4, ncol=2)
out[1,1] <- exp(-1*psclust1cox$coefficients[1])
out[1,2] <- psclust1cox$coefficients[5]
out[2,1] <- exp(psclust2cox$coefficients[1])
out[2,2] <- psclust2cox$coefficients[5]
out[3,1] <- exp(chclust1cox$coefficients[1])
out[3,2] <- chclust1cox$coefficients[5]
out[4,1] <- exp(-1*chclust2cox$coefficients[1])
out[4,2] <- chclust2cox$coefficients[5]
#out[5,1] <- exp(-1*cshclust1cox$coefficients[1])
#out[5,2] <- cshclust1cox$coefficients[5]
#out[6,1] <- exp(cshclust2cox$coefficients[1])
#out[6,2] <- cshclust2cox$coefficients[5]

colnames(out) <- c("HR", "p-value")
write.csv(out, "Table6.csv")

