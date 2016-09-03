set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=6, ncol=20)

for (i in 1:6) {
  print(i)
  pe <- c(4, 8, 12, 16, 20, 24)
  out[i,] <- sim4.test(pe[i])
}

write.csv(out, "Table4.csv", row.names=FALSE)
