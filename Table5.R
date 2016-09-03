set.seed(12221)
library(sparcl)
library(doMC)

registerDoMC(cores=multicore:::detectCores())

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



#Tuning parameter
tune <- 42.55639

mc <- 0.3

misclass <- foreach(k=1:1000, .combine=rbind) %dopar% {

#Start entire matrix
genes =matrix(rnorm(1000000, 0), 5000, 200)
y.true <- c(rep(2,100), rep(1,100))
yval <- c(rep(2,100), rep(1,100))

#Make genes 1-50 split between mean 1, 2
genes[1:50, 1:100] = rnorm(5000, 1)
genes[1:50, 101:200] = rnorm(5000, 2)

rnum <- runif(100)
ind1a <- which(rnum<mc)
yval[1:100][ind1a] <- 1

rnum <- runif(100)
ind1b <- which(rnum<mc)
yval[101:200][ind1a] <- 2

yval <- yval-1
#to make the permute function work

#Make 40% of genes 51-100 have mean 2
rnum <- runif(200)
ind1 <- which(rnum<0.4)
genes[51:100,ind1] <- rnorm(50*length(ind1), 2)

#Make 70% of genes 101-200 have mean 0.5
rnum <- runif(200)
ind2 <- which(rnum<0.7)
genes[101:200,ind2] <- rnorm(100*length(ind2), 0.5)

#Make 30% of genes 201-300 have mean 1.5
rnum <- runif(200)
ind3 <- which(rnum<0.3)
genes[201:300,ind3] <- rnorm(100*length(ind3), 1.5)

#Make the matrix into a data frame
gene = data.frame(genes)
#To transpose the matrix with genes/patients on opposite sides
tran=t(gene)
trans=data.frame(tran)


###############################
#To run supervised sparse 2means
tstat = rep(NA,5000)
for (i in c(1:5000))
{tstat[i] = abs(as.numeric(t.test(trans[,i]~yval)$statistic))}

weights = rep(0, 5000)
toptstat <- order(tstat, decreasing=TRUE)
n.ws <- floor(sqrt(ncol(trans)))
for (i in toptstat[1:n.ws])
{weights[i] = (1/sqrt(n.ws))}

sspcl.out <- KMeansSparseCluster2(trans, K=2, ws0=weights, wbounds=tune)
sspcl.table <- table(y.true,sspcl.out$Cs)



###############################
#To run sparse 2means
permout <- KMeansSparseCluster.permute(trans, K = 2, nperm=5)

sparse <- KMeansSparseCluster(trans, K = 2 ,wbounds=permout$bestw)
sparse.table <- table(y.true,sparse[[1]]$Cs)


###############################
#To run supervised 2means
trans2 <- trans[,toptstat[1:n.ws]]
supervised <- kmeans(trans2, 2, nstart=10)
supervised.table <- table(y.true,supervised$cluster)


###############################
#Based on principal components
pca <- prcomp(~., data=trans, scale=TRUE)
scores <- predict(pca)
pcaclust <- kmeans(scores[,1:n.ws], 2, nstart=10)
pca.table <- table(y.true,pcaclust$cluster)

out <- rep(NA, 4)

#misclass[t,1] <- min((sspcl.table[1,1]+sspcl.table[2,2]), (sspcl.table[1,2]+sspcl.table[2,1]))
#misclass[t,2] <- min((sparse.table[1,1]+sparse.table[2,2]), (sparse.table[1,2]+sparse.table[2,1]))
#misclass[t,3] <- min((supervised.table[1,1]+supervised.table[2,2]), (supervised.table[1,2]+supervised.table[2,1]))
#misclass[t,4] <- min((pca.table[1,1]+pca.table[2,2]), (pca.table[1,2]+pca.table[2,1]))
out[1] <- min((sspcl.table[1,1]+sspcl.table[2,2]), (sspcl.table[1,2]+sspcl.table[2,1]))
out[2] <- min((sparse.table[1,1]+sparse.table[2,2]), (sparse.table[1,2]+sparse.table[2,1]))
out[3] <- min((supervised.table[1,1]+supervised.table[2,2]), (supervised.table[1,2]+supervised.table[2,1]))
out[4] <- min((pca.table[1,1]+pca.table[2,2]), (pca.table[1,2]+pca.table[2,1]))
out
}

colnames(misclass) <- c("SupSparse", "Sparse", "Supervised", "PCA")


#write.csv(misclass, file="Sim1000_Summary.csv")
a <- 2














###############################
##Summarizing all the simulations
###############################

misclassSUM <- matrix(NA,2,4)
colnames(misclassSUM) <- c("SupSparse", "Sparse", "Supervised", "PCA")
rownames(misclassSUM) <- c("Mean", "SD")
misclassSUM[1,] <- colMeans(misclass)
misclassSUM[2,1] <- sd(misclass[,1])
misclassSUM[2,2] <- sd(misclass[,2])
misclassSUM[2,3] <- sd(misclass[,3])
misclassSUM[2,4] <- sd(misclass[,4])
write.csv(misclassSUM, file="Table5.csv")






##Standard error
#misclassSUM <- read.csv("Sim1000_Summary.csv")
