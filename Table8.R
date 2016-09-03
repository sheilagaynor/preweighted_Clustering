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
#  	x1 <- permx[[k]][y==0,]
#	x2 <- permx[[k]][y==1,]
#  	n1 <- nrow(x1)
#  	n2 <- nrow(x2)
#  	x1.bar <- t(x1) %*% rep(1/n1, n1)
#  	x2.bar <- t(x2) %*% rep(1/n2, n2)
#  	x1.sd <- sqrt((t(x1)-as.vector(x1.bar))^2 %*% rep(1/(n1-1),n1))
#  	x2.sd <- sqrt((t(x2)-as.vector(x2.bar))^2 %*% rep(1/(n2-1),n2))
#  	sd.pool <- sqrt(x1.sd^2/n1 + x2.sd^2/n2)
#  	tstat <- (x1.bar-x2.bar) / sd.pool
#  	ws.star <- rep(0, ncol(x))
#  	toptstat <- order(abs(tstat), decreasing=TRUE)
#  	ws.star[toptstat[1:nws]] <- 1/sqrt(nws)
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




set.seed(4321)

library(sparcl)
library(pamr)
library(rms)
library(survival)
load("../data/sspcl_funcs.RData")


#Read the data:
pollack.data <- read.csv("pollack.dat", header=FALSE)
pollack.surv <- read.csv("pollack.surv", header=FALSE)
pollack.x <- t(pollack.data)
pollack.y <- pollack.surv[,2]


#Create a training set and a test set:
n.case <- sum(pollack.y==1)
n.control <- sum(pollack.y==0)
ftr <- rep(FALSE, nrow(pollack.x))
ftr[pollack.y==0][sample(1:n.control, size=floor(n.control/2))] <- TRUE
ftr[pollack.y==1][sample(1:n.case, size=floor(n.case/2))] <- TRUE
ftst <- !ftr


#Scale the data set:
ftr.center <- colMeans(pollack.x[ftr,])
ftr.scale <- apply(pollack.x[ftr,], 2, sd)
pollack.x <- scale(pollack.x, center=ftr.center, scale=ftr.scale)


#Make table for output
fsummary <- matrix(-10,4,2)
rownames(fsummary) <- c("SupSparse", "Sparse", "Supervised", "PCA")
colnames(fsummary) <- c("HR", "Pval")


####################
##Supervised sparse

####################

#Apply supervised sparse clustering to the training set:
coxstat = rep(NA,ncol(pollack.x))
for (i in c(1:ncol(pollack.x)))
{coxstat[i] = abs(as.numeric(summary(coxph(Surv(pollack.surv[ftr,1],pollack.surv[ftr,2])~pollack.x[ftr,i]))[[7]][4]))}
weights = rep(0, ncol(pollack.x))
topcoxstat <- order(coxstat, decreasing=TRUE)
n.ws <- floor(sqrt(ncol(pollack.x)))
weights[topcoxstat[1:n.ws]] <- 1/sqrt(n.ws)


perm <- KMeansSparseCluster.permute2(pollack.x[ftr,], pollack.y[ftr], K=2, ws0=weights, nperms=5)
sspcl.ftr <- KMeansSparseCluster2(pollack.x[ftr,], K=2, ws0=weights, wbounds=perm$bestw, maxiter=10)

sspcl.pamr <- pamr.train(list(x=t(pollack.x[ftr,]), y=factor(sspcl.ftr$Cs)))
sspcl.pamr.cv <- pamr.cv(sspcl.pamr,
                         list(x=t(pollack.x[ftr,]), y=factor(sspcl.ftr$Cs)))
cur.ndx <- max(which(sspcl.pamr.cv$loglik==max(sspcl.pamr.cv$loglik)))
sspcl.yhat <- pamr.predict(sspcl.pamr, t(pollack.x[ftst,]),
                           sspcl.pamr.cv$threshold[cur.ndx])


#Evaluate the association between the clusters and the outcome variable on the test data:
sspcl.cox <- summary(coxph(Surv(pollack.surv[ftst,1],pollack.surv[ftst,2]) ~ sspcl.yhat))
fsummary[1,1] <- exp(-1*sspcl.cox$coefficients[1])
fsummary[1,2] <- sspcl.cox$coefficients[5]






#################
##Sparse

#################

#Apply sparse clustering to the training set:

permout <- KMeansSparseCluster.permute(pollack.x[ftr,], K = 2, nperms=5)
sparse.ftr <- KMeansSparseCluster(pollack.x[ftr,], K = 2, wbounds=permout$bestw,maxiter=10)[[1]]

sparse.pamr <- pamr.train(list(x=t(pollack.x[ftr,]), y=factor(sparse.ftr$Cs)))
sparse.pamr.cv <- pamr.cv(sparse.pamr,
                         list(x=t(pollack.x[ftr,]), y=factor(sparse.ftr$Cs)))
cur.ndx <- max(which(sparse.pamr.cv$loglik==max(sparse.pamr.cv$loglik)))
sparse.yhat <- pamr.predict(sparse.pamr, t(pollack.x[ftst,]),
                           sparse.pamr.cv$threshold[cur.ndx])

#Evaluate the association between the clusters and the outcome variable on the test data:
sparse.cox <- summary(coxph(Surv(pollack.surv[ftst,1],pollack.surv[ftst,2]) ~ sparse.yhat))
fsummary[2,1] <- exp(-1*sparse.cox$coefficients[1])
fsummary[2,2] <- sparse.cox$coefficients[5]


#sparse.yhat <- apply(cbind(cl1b.dist, cl2b.dist), 1, which.min)
##This only returns one group

#################
##Supervised
#################

#Apply supervised clustering to the training set:

supervised.ftr <- kmeans(pollack.x[ftr,topcoxstat[1:n.ws]], 2, nstart=10)

supervised.pamr <- pamr.train(list(x=t(pollack.x[ftr,]), y=factor(supervised.ftr$cluster)))
supervised.pamr.cv <- pamr.cv(supervised.pamr,
                         list(x=t(pollack.x[ftr,]), y=factor(supervised.ftr$cluster)))
cur.ndx <- max(which(supervised.pamr.cv$loglik==max(supervised.pamr.cv$loglik)))
supervised.yhat <- pamr.predict(supervised.pamr, t(pollack.x[ftst,]),
                           supervised.pamr.cv$threshold[cur.ndx])

#Evaluate the association between the clusters and the outcome variable on the test data:
supervised.cox <- summary(coxph(Surv(pollack.surv[ftst,1],pollack.surv[ftst,2]) ~ supervised.yhat))
fsummary[3,1] <- exp(supervised.cox$coefficients[1])
fsummary[3,2] <- supervised.cox$coefficients[5]


#################
##PCA
#################

#Apply PCA clustering to the training set:
pca <- prcomp(~., data=data.frame(pollack.x[ftr,]), scale=TRUE)
scores.ftr <- predict(pca)
pca.ftr <- kmeans(scores.ftr[,1:10], 2, nstart=10)

pca.pamr <- pamr.train(list(x=t(pollack.x[ftr,]), y=factor(pca.ftr$cluster)))
pca.pamr.cv <- pamr.cv(pca.pamr,
                         list(x=t(pollack.x[ftr,]), y=factor(pca.ftr$cluster)))
cur.ndx <- max(which(pca.pamr.cv$loglik==max(pca.pamr.cv$loglik)))
pca.yhat <- pamr.predict(pca.pamr, t(pollack.x[ftst,]),
                           pca.pamr.cv$threshold[cur.ndx])

#Evaluate the association between the clusters and the outcome variable on the test data:
pca.cox <- summary(coxph(Surv(pollack.surv[ftst,1],pollack.surv[ftst,2]) ~ pca.yhat))
fsummary[4,1] <- exp(pca.cox$coefficients[1])
fsummary[4,2] <- pca.cox$coefficients[5]





write.csv(fsummary, file="Table8.csv")

