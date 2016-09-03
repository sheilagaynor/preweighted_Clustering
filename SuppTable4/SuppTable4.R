set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=5, ncol=16)

for (i in 1:5) {
  print(i)
  bv <- c(20, 30, 50, 100, 200)
  out[i,] <- sim7.test(bv[i])
}

write.csv(out, "SuppTable4.csv", row.names=FALSE)
