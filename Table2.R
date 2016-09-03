set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=7, ncol=20)

for (i in 1:7) {
  print(i)
  sig <- c(1, 2, 2.5, 3, 4, 5, 6)
  out[i,] <- sim2.test(sig[i])
}

write.csv(out, "Table2.csv", row.names=FALSE)
