args <- commandArgs(T)

cat("make sure:\n    y and x must have the same sample id\n    train_y must have only two levels and will change to 0 and 1\n    test_y which not in train_y will change to 2\n    set marker_num 0 to compute automaticity\nnote:if x and y have different levels could lead to errors\n")

if (length(args) != 9) {
  stop("Rscript *.R [train_x] [train_y] [test_x] [test_y] [cv_fold] [cv_step] [cv_time] [marker_num] [prefix]\n")
}

train.x <- args[1]
train.y <- args[2]
test.x <- args[3]
test.y <- args[4]
cv.fold <- as.numeric(args[5])
cv.step <- as.numeric(args[6])
cv.time <- as.numeric(args[7])
marker.num <- as.numeric(args[8])
prefix <- args[9]

# package
library(randomForest)

args <- commandArgs(F)
## SD <- dirname(sub("--file=", "", args[grep("--file=", args)]))
SD <- 
# function
source(paste0(SD, "/rfcv1.R"))
source(paste0(SD, "/ROC.R"))

# data
train.x <- t(read.table(train.x))
train.y <- as.factor(read.table(train.y)[, 1])
train.l <- levels(train.y)
levels(train.y) <- 0:1

test.x <- t(read.table(test.x))
test.y <- as.factor(read.table(test.y)[, 1])
test.l <- levels(test.y)
levels(test.y) <- pmatch(test.l, train.l) - 1
test.y <- factor(test.y, 0:2)

# crossvalidation
pdf.dir <- paste0(prefix, "_randomForest.pdf")
pdf(pdf.dir, width = 35, height = 7)
par(mfrow = c(1, 5))

set.seed(0)
train.cv <- replicate(cv.time, rfcv1(train.x, train.y, cv.fold = cv.fold, step = cv.step), simplify = F)
error.cv <- sapply(train.cv, "[[", "error.cv")
error.cv.rm <- rowMeans(error.cv)
# id <- error.cv.rm < min(error.cv.rm) + diff(range(error.cv.rm))/20
id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
error.cv[id, ]
if (marker.num == 0) {
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
}
matplot(train.cv[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cv.time), main = paste("select", marker.num, "Vars"), xlab = "Number of vars", 
  ylab = "CV Error", lty = 1)
lines(train.cv[[1]]$n.var, error.cv.rm, lwd = 2)
abline(v = marker.num, col = "pink", lwd = 2)

# pick marker by corossvalidation
marker.t <- table(unlist(lapply(train.cv, function(x) {
  lapply(x$res, "[", 1:marker.num)
})))
marker.t <- sort(marker.t, d = T)
names(marker.t) <- colnames(train.x)[as.numeric(names(marker.t))]
marker.dir <- paste0(prefix, "_marker.txt")
write.table(marker.t, marker.dir, col.names = F, sep = "\t", quote = F)
marker.p <- names(marker.t)[1:marker.num]

# train model
set.seed(0)
train.rf <- randomForest(train.x[, marker.p], train.y, importance = T)
train.p <- predict(train.rf, type = "prob")
boxplot(train.p[, 2] ~ train.y, col = 2:3, main = "Probability", names = train.l)
pr.dir <- paste0(prefix, "_train_probability.txt")
write.table(train.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)

# train ROC
plot_roc(train.y, train.p[, 2])

# test predict
test.p <- predict(train.rf, test.x, type = "prob")
pr.dir <- paste0(prefix, "_test_probability.txt")
write.table(test.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)

# predict plot
p.col <- ifelse(is.na(test.y), 4, as.numeric(test.y) + 1)
plot(rank(test.p[, 2]), test.p[, 2], col = p.col, pch = 16, xlab = "", ylab = "Probability", main = "Testset")
txt <- train.l
if (length(test.l) > 2) {
  txt <- c(txt, "the rest")
}
legend("bottomright", txt, col = 2:4, pch = 16)
abline(h = 0.5)

# test ROC
plot_roc(test.y, test.p[, 2])
dev.off() 
