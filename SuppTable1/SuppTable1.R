set.seed(12221)

source("Nowak_Sims_Supp1.R")

out <- matrix(nrow=7, ncol=28)

for (i in 1:7) {
  print(i)
  ali <- c(0.1,0.25,0.5,1,2,4,6)
  out[i,] <- sim2.test(ali[i])
}

write.csv(out, "SuppTable1.csv", row.names=FALSE)
