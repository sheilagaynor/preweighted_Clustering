set.seed(12121)
library(sparcl)
library(glmnet)
library(survival)
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
oppera.x <- oppera.raw[,12:127]
oppera.y <- Surv(oppera.raw$followupyears, oppera.raw$incidentcase)
oppera.y[is.na(oppera.y[,1]),1] <- 0
oppera.y[is.na(oppera.y[,2]),2] <- 0

#Create a training set and a test set:
n.case <- sum(oppera.y[,2]==1)
n.control <- sum(oppera.y[,2]==0)
otr <- rep(FALSE, nrow(oppera.x))
otr[oppera.y[,2]==0][sample(1:n.control, size=floor(n.control/2))] <- TRUE
otr[oppera.y[,2]==1][sample(1:n.case, size=floor(n.case/2))] <- TRUE
otst <- !otr


#Scale the data set:
otr.center <- colMeans(oppera.x[otr,])
otr.scale <- apply(oppera.x[otr,], 2, sd)
oppera.x <- scale(oppera.x, center=otr.center, scale=otr.scale)

####################
##Supervised sparse

####################

#Apply supervised sparse clustering to the training set:

zstat = rep(NA,ncol(oppera.x))
for (i in c(1:ncol(oppera.x)))
{zstat[i] = abs(summary(coxph(oppera.y[otr,]~oppera.x[otr,i]))$coefficients[4])}
weights = rep(0, ncol(oppera.x))
toptstat <- order(zstat, decreasing=TRUE)
n.ws <- floor(sqrt(ncol(oppera.x)))
weights[toptstat[1:n.ws]] <- 1/sqrt(n.ws)



perm <- KMeansSparseCluster.permute2(oppera.x[otr,], oppera.y[otr], K=2, ws0=weights, nperms=5)
sspcl <- KMeansSparseCluster2(oppera.x[otr,], K=2, ws0=weights, wbounds=perm$bestw, maxiter=10)

sspcl.glmnet <- glmnet(oppera.x[otr,], factor(sspcl$Cs), family="binomial")
sspcl.glmnet.cv <- cv.glmnet(oppera.x[otr,], factor(sspcl$Cs), family="binomial")
sspcl.yhat <- predict(sspcl.glmnet, newx=oppera.x[otst,], type="class",
                      s=sspcl.glmnet.cv$lambda.1se)

#################
##Sparse

#################

#Apply sparse clustering to the training set:

permout <- KMeansSparseCluster.permute(oppera.x[otr,], K = 2, nperms=5)

sparse <- KMeansSparseCluster(oppera.x[otr,], K = 2, wbounds=permout$bestw,
                              maxiter=10)[[1]]

sparse.glmnet <- glmnet(oppera.x[otr,], factor(sparse$Cs), family="binomial")
sparse.glmnet.cv <- cv.glmnet(oppera.x[otr,], factor(sparse$Cs), family="binomial")
sparse.yhat <- predict(sparse.glmnet, newx=oppera.x[otst,], type="class",
                      s=sparse.glmnet.cv$lambda.1se)

#################
##Supervised
#################

#Apply supervised clustering to the training set:

supervised <- kmeans(oppera.x[otr,toptstat[1:n.ws]], 2, nstart=10)

supervised.glmnet <- glmnet(oppera.x[otr,], factor(supervised$cluster), family="binomial")
supervised.glmnet.cv <- cv.glmnet(oppera.x[otr,], factor(supervised$cluster), family="binomial")
supervised.yhat <- predict(supervised.glmnet, newx=oppera.x[otst,], type="class",
                      s=supervised.glmnet.cv$lambda.1se)

#################
##PCA
#################

#Apply PCA clustering to the training set:

pca <- prcomp(~., data=data.frame(oppera.x[otr,]), scale=TRUE)
scores <- predict(pca)
pca <- kmeans(scores[,1:5], 2, nstart=10)

pca.glmnet <- glmnet(oppera.x[otr,], factor(pca$cluster), family="binomial")
pca.glmnet.cv <- cv.glmnet(oppera.x[otr,], factor(pca$cluster), family="binomial")
pca.yhat <- predict(pca.glmnet, newx=oppera.x[otst,], type="class",
                      s=pca.glmnet.cv$lambda.1se)

sspcl.cox <- summary(coxph(oppera.y[otst,]~factor(sspcl.yhat)))

sparse.cox <- summary(coxph(oppera.y[otst,]~factor(sparse.yhat)))

sclust.cox <- summary(coxph(oppera.y[otst,]~factor(supervised.yhat)))

pca.cox <- summary(coxph(oppera.y[otst,]~factor(pca.yhat)))

out <- matrix(nrow=4, ncol=2)
out[1,1] <- exp(sspcl.cox$coefficients[1])
out[1,2] <- sspcl.cox$coefficients[5]
out[2,1] <- exp(-1*sparse.cox$coefficients[1])
out[2,2] <- sparse.cox$coefficients[5]
out[3,1] <- exp(sclust.cox$coefficients[1])
out[3,2] <- sclust.cox$coefficients[5]
out[4,1] <- exp(-1*pca.cox$coefficients[1])
out[4,2] <- pca.cox$coefficients[5]

colnames(out) <- c("HR", "p-value")
write.csv(out, "Table7.csv", row.names=c("SupSparse", "Sparse", "Supervised",
                             "PCA"))
