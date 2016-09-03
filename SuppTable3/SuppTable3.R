set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=5, ncol=16)

for (i in 1:5) {
  print(i)
  bv <- c(0.5, 0.75, 1, 2, 3)
  out[i,] <- sim6.test(bv[i])
}

write.csv(out, "SuppTable3.csv", row.names=FALSE)
