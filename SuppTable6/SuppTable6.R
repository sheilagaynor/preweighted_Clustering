set.seed(12221)

source("Nowak_Sims2.R")

out <- matrix(nrow=6, ncol=20)

for (i in 1:6) {
  print(i)
  bv <- c(0, 5, 10, 15, 20, 25)
  out[i,] <- sim1.test(bv[i])
}

write.csv(out, "SuppTable6.csv", row.names=FALSE)
