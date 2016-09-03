set.seed(12221)

source("Nowak_Sims.R")

out <- matrix(nrow=3, ncol=20)

for (i in 1:3) {
  print(i)
  na <- c(6, 8, 10)
  out[i,] <- sim3.test(na[i])
}

write.csv(out, "Table3.csv", row.names=FALSE)
