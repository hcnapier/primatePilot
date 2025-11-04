testNBin <- rnbinom(100, size = 10, mu = 8)
hist(testNBin, breaks = 10)
z <- mean(testNBin)
num <- sum(testNBin)
s <- num/z
testNBin_norm <- log((testNBin/s) + 1)
hist(testNBin_norm)
var(testNBin)
var(testNBin_norm)

z2 <- 14
s2 <- num/z2
testNBin_norm2 <- log((testNBin/s2)+1)
hist(testNBin_norm2, breaks = 10)
var(testNBin_norm2)

var(testNBin_norm)/mean(testNBin_norm)
var(testNBin_norm2)/mean(testNBin_norm2)
sd(testNBin_norm)/mean(testNBin_norm)
sd(testNBin_norm2)/mean(testNBin_norm2)

testNBin2 <- rnbinom(200, size = 2, mu = 2)
hist(testNBin2)
z_dist2 <- mean(testNBin2)
num_dist2 <- sum(testNBin2)
s_dist2 <- num_dist2/z_dist2
testNBin2_norm <- log(testNBin/s_dist2 + 1)
hist(testNBin2_norm)
var(testNBin2)
var(testNBin2_norm)

testNBin2_norm2 <- log(testNBin2/s + 1)
hist(testNBin2_norm2)
var(testNBin2_norm2)

z_avg <- mean(c(testNBin2, testNBin))
num_avg <- sum(testNBin, testNBin2)
s_avg <- num_avg/z_avg
testNBin2_normAvg <- log(testNBin2/s_avg + 1)
testNBin_normAvg <- log(testNBin/s_avg + 1)
hist(testNBin2_normAvg)
hist(testNBin_normAvg)
