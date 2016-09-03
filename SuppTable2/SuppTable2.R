set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=5, ncol=20)

for (i in 1:5) {
  print(i)
  rv <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  out[i,] <- sim5.test(rv[i])
}

write.csv(out, "SuppTable2.csv", row.names=FALSE)
