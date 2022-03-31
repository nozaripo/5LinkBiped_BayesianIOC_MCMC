install.packages("R.matlab")

library("R.matlab")
library(glmnet)
library(corrplot)

setwd('C:/Users/nozaripo/OneDrive - University of Southern California/OptimTraj-master/demo/fiveLinkBiped - IOC')


Data = readMat('Fukuchi_Features_DTW_RMS.mat', fixNames=TRUE, drop=c("singletonLists"),
        sparseMatrixClass=c("Matrix"), verbose=FALSE)


Features.Matrix = Data$Cost.Components.Subjects
DTW.Vector      = Data$DTW.Sum.Subjects
RMSE.Vector     = Data$RMSE.Sum.Subjects

Features.scaled = Data$X.scaled

Corr_Matrix = cor(Features.scaled, method = c("pearson", "kendall", "spearman"))

corrplot(Corr_Matrix, method = "circle")

colnames(Corr_Matrix) <- c(Criteria.Names)





cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(Features.scaled, 0.95)

corrplot(Corr_Matrix, p.mat = res1[[1]], method = 'circle', type = 'lower', insig='blank',
         order = 'AOE', diag = TRUE)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))



corrplot(Corr_Matrix, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=TRUE)

corrplot.mixed(Corr_Matrix, order = 'AOE')




X = Data$X.scaled
#X = Data$Cost.Components
y = Data$RMSE.Sum
#X = scale(X, center = TRUE, scale = TRUE)

# perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(X, y, alpha = 1)
# find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
# Lasso on X and y
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda, standardize = FALSE)
Lasso_Coefficients_DTW_Sum = matrix(best_model$beta)


y = Data$RMSE.Sum
cv_model <- cv.glmnet(X, y, alpha = 1)
# find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
# Lasso on X and y
best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda, standardize = FALSE)
Lasso_Coefficients_RMS_Sum = matrix(best_model$beta)



Lasso_Coefficients_DTW = matrix(0L, nrow = 16, ncol = 10)
Lasso_Coefficients_RMS = matrix(0L, nrow = 16, ncol = 10)

for (i in 1:10){
  ## LASSO FOR DTW MEASURE ##
  X = Features.Matrix[[i]][[1]]
  y = DTW.Vector[[i]][[1]]
  
  X = scale(X, center = TRUE, scale = TRUE)
#  y1 = scale(y, center = FALSE, scale = TRUE)
  
  
#  library(dplyr)
#  as.matrix(X) %>% summarise_if(is.numeric, max)
  
  
  
  # perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(X, y, alpha = 1)
  # find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  # Lasso on X and y
  best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda, standardize = FALSE)
  Lasso_Coefficients_DTW[,i] = matrix(best_model$beta)
  
#  grid=10^seq(10,-4, length =100)
#  out=glmnet (X,y,alpha=1, lambda=grid)
#  lasso.coef=predict (out ,type=" coefficients",s= best_lambda)
  
  
  
  ## LASSO FOR RMS MEASURE ##
  y = RMSE.Vector[[i]][[1]]
#  y = scale(y, center = TRUE, scale = TRUE)
  
  # perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(X, y, alpha = 1)
  # find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  # Lasso on X and y
  best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda, standardize = FALSE)
  Lasso_Coefficients_RMS[,i] = matrix(best_model$beta)
  
  
}















########################  PCA FOR sUBJECTS  ###############################

V=var.explained.spca.15 = matrix(0L, nrow = numSteps, ncol = num.Component)

#####################################################################
for (i in seq(0, numSteps, by=1)){
  
  sparsity = 5e-2*i/numSteps
  
  aa6 = robspca(X, k = num.Component, alpha = sparsity, beta = 1e-04, gamma = 100,
                center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-05,
                verbose = TRUE)  
  
  for (nComp in seq(1, num.Component, by=1)){
    
    var.explained.spca.15[i,nComp] = sum(aa6$sdev[1:nComp]^2)
    
  }
}


var.exp.normal.15 = var.explained.spca.15/max(var.explained.spca.15)

sp = seq(0, ((numSteps-1)/numSteps)*5e-2, by=5e-2/numSteps)


matplot(sp, var.exp.normal.15, type = "l", lty = 1:1, lwd = 2, pch = NULL,
        col = 1:16, cex = NULL, bg = NA,
        main = "Variance Explained for Different PC Numbers and Sparsity", xlab = "Sparsity", ylab = "VAF", 
        xlim = c(0,.06), ylim = c(.4,1), add = FALSE)
legend("topright",                              # Add legend to plot
       legend = c("#PCs = 1", "#PCs = 2", 
                  "#PCs = 3" , "#PCs = 4" , "#PCs = 5" , "#PCs = 6" , "#PCs = 7"),
       col = 1:16,
       lwd = 2)

View(var.exp.normal.15)
























View(best_model)
View(coef(best_model))
Features.Matrix[[31]][[1]]

length(unlist(Data$Cost.Components.Subjects))





dim(Data)
Last_Data = unlist(Data$Cost.Components.Subjects [[31]][[1]])
Last_Data = unlist(Data$DTW.Sum.Subjects [[31]][[1]])
Last_Data = Data$Cost.Components

Last_Data = matrix(Last_Data,length(Last_Data)/16,16)
length(Last_Data)/16
dim(Last_Data)[1]
