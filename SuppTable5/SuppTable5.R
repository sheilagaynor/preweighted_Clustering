set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=6, ncol=20)

for (i in 1:6) {
  print(i)
  bv <- c(0.5, 0.75, 1, 2, 3, 6)
  out[i,] <- sim1.test(bv[i])
}

write.csv(out, "SuppTable5.csv", row.names=FALSE)
