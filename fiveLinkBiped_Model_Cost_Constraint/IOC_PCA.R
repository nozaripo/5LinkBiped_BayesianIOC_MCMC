library(sparsepca)

setwd('C:/Users/nozaripo/OneDrive - University of Southern California/OptimTraj-master/demo/fiveLinkBiped - IOC')

X = read.csv('Cost_Comp_Eval.txt',header=F)
X = read.csv('Cost_Comp_Eval_Traj.txt',header=F)
X = read.csv('Cost_Components_Eval_Fakuchi.txt',header=F)
X = read.csv('Cost_Components_Eval_Fakuchi_AllData_NoDynamics.txt',header=F)
X = read.csv('Cost_Components_Eval_Fakuchi_AllData.txt',header=F)

#X <- subset (X, select = -5)


#y = read.csv('C:/Users/nozaripo/OneDrive - University of Southern California/PCA_Y.txt',header=F)
y = read.csv('Cost_Deviation_Eval_y_forPLSR.txt',header=F)


X.scaled = scale(X, center = TRUE, scale = TRUE)
y.scaled = scale(y, center = TRUE, scale = TRUE)


num.Component = dim(X)[2]
#TSS = sum((y-mean(t(y)))^2)


numSteps =1000


resid.test  <- matrix(0L, nrow = numSteps, ncol = 10)
resid.train <- matrix(0L, nrow = numSteps, ncol = 10)
mean.resid  <- matrix(0L, nrow = numSteps, ncol = 1)
Var.Exp     <- matrix(0L, nrow = numSteps, ncol = 1)

var.explained = matrix(0L, nrow = numSteps, ncol = 1)
total.var     = matrix(0L, nrow = numSteps, ncol = 1)

var.explained.spca = matrix(0L, nrow = numSteps, ncol = num.Component)







############# STOPPPPPPPPPPPPPPPPPPPPP

#####################
for (i in seq(0, numSteps, by=1)){

sparsity = 5e-1*i/numSteps


for (nComp in seq(1, num.Component, by=1)){
  
aa6 = robspca(X, k = nComp, alpha = sparsity, beta = 1e-04, gamma = 100,
          center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-04,
          verbose = TRUE)  
# aa6=spca(X, k = nComp, alpha = sparsity, beta = 1e-04, center = TRUE,
#          scale = TRUE, max_iter = 3000, tol = 1e-06, verbose = TRUE)

#View(aa6)

#View(aa6$loadings)

var.explained.spca[i,nComp] = sum(aa6$sdev^2)
}
}

var.exp.normal = var.explained.spca/max(var.explained.spca)


# sp = seq(1e-5, 1e-5+((numSteps-1)/numSteps)*(2e-1 - 1e-5), by=(2e-1 - 1e-5)/numSteps)
sp = seq(0, ((numSteps-1)/numSteps)*5e-1, by=5e-1/numSteps)


matplot(sp, var.exp.normal, type = "l", lty = 1:1, lwd = .5, pch = NULL,
        col = 1:5, cex = NULL, bg = NA,
        xlab = NULL, ylab = NULL, xlim = c(0,.11), ylim = c(.4,1), add = FALSE)



sparsity = 5e-2*153/numSteps

aa6 = robspca(X, k = 5, alpha = sparsity, beta = 1e-04, gamma = 100,
              center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-05,
              verbose = TRUE)
sum(aa6$sdev[1:5]^2)/max(var.explained.spca)
sum(aa6$sdev^2)/max(var.explained.spca)



# with bigger sparsity
sparsity = 5e-2*300/numSteps

aa6 = robspca(X, k = 5, alpha = sparsity, beta = 1e-04, gamma = 100,
              center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-05,
              verbose = TRUE)
sum(aa6$sdev[1:5]^2)/max(var.explained.spca)
sum(aa6$sdev^2)/max(var.explained.spca)










var.explained.spca.15 = matrix(0L, nrow = numSteps, ncol = num.Component)

#####################################################################
for (i in seq(0, numSteps, by=1)){
  
  sparsity = 5e-2*i/numSteps
  
  aa6 = robspca(X, k = num.Component, alpha = sparsity, beta = 1e-04, gamma = 100,
                center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-04,
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

##################################################################

sparsity = 5e-2*118/numSteps
sparsity = 5e-2*132/numSteps

sparsity = 0

sparsity = 5e-2*200/numSteps

sparsity = 5e-2*150/numSteps
sparsity = 5e-2*0/numSteps

aa6.15 = robspca(X, k = nComp, alpha = sparsity, beta = 1e-04, gamma = 100,
              center = TRUE, scale = TRUE, max_iter = 100000, tol = 1e-05,
              verbose = TRUE)
sum(aa6.15$sdev[1:4]^2)/max(var.explained.spca.15)

View(aa6.15$loadings)

sum(aa6.15$sdev^2)/max(var.explained.spca.15)

# attach(aa6)
# 
# weights = loadings[,1:2]
# 
# 
# Reg.data = as.matrix(X.scaled)%*%as.matrix(weights)
# 
# 
# 
# for (j in seq(1, 10, by=1)){
#   
#   test = seq(40*(j-1)+1,40*j, by=1)
#   
#   X.test  = Reg.data[test,]
#   X.train = Reg.data[-test,]
#   
#   y.test  = y[test,1]
#   y.train = y[-test,1]
#   
#   
#   feat1 = X.train[,1]
#   feat2 = X.train[,2]
#   
#   feat = cbind(feat1,feat2)
#   
#   
#   lm.fit = lm(y.train~feat-1)
#   
#   resid.train[i,j] = sum(lm.fit$residuals^2)
#   
#   feat1 = X.test[,1]
#   feat2 = X.test[,2]
#   
#   feat = cbind(feat1,feat2)
#   
#   yy = predict(lm.fit, data.frame(feat=feat))
#   
#   resid.test[i,j] = sum((yy-y.test)^2)
# }
# 
# mean.resid[i,1] = mean(resid.test[i,])
# 
# Var.Exp[i,1] = 1-mean.resid[i,1]/TSS
# }
  



library(caret)


mean.resid  <- matrix(0L, nrow = numSteps, ncol = 1)

numSteps=100
for (i in seq(0, numSteps, by=1)){
  
  sparsity = 5e-2*i/numSteps
  
  aa6 = robspca(X, k = num.Component, alpha = sparsity, beta = 1e-04, gamma = 100,
                center = TRUE, scale = TRUE, max_iter = 1000, tol = 1e-05,
                verbose = TRUE)  
  
  
  X.Principal = as.matrix(X.scaled)%*%as.matrix(aa6$loadings)

  Results.lm <- lm(y.scaled ~ X.Principal)
  
  
#  train_control <- trainControl(method="cv", number=10)
#  # train the model
#  model <- train(y.scaled ~ X.Principal, trControl=train_control)
#  # summarize results
  
#  print(model)
  
  
#  summary(Results.lm)
  
  mean.resid[i,1] = sqrt(mean(Results.lm$residuals^2))
  
#  for (nComp in seq(1, num.Component, by=1)){
    
#    var.explained.spca.15[i,nComp] = sum(aa6$sdev[1:nComp]^2)
    
#  }

}

min(mean.resid)

sd(mean.resid)


View(mean.resid)



sparsity = 5e-2*6/numSteps

aa6 = robspca(X, k = num.Component, alpha = sparsity, beta = 1e-04, gamma = 100,
              center = TRUE, scale = TRUE, max_iter = 1000, tol = 1e-05,
              verbose = TRUE)

View(aa6$loadings)

#####################
# for (i in seq(1, numSteps, by=1)){
#   
#   sparsity = 1e-5+((i-1)/numSteps)*(1e-2 - 1e-5)  
#   
#   
#   # for (Ncomp in seq(1, 15, by=1)){
#   # Ncomp = 2
#   
#   aa6=spca(X.scaled, k = NULL, alpha = sparsity, beta = 1e-04, center = TRUE,
#            scale = FALSE, max_iter = 3000, tol = 1e-06, verbose = TRUE)
#   
#   #View(aa6)
#   
#   #View(aa6$loadings)
#   
#   var.temp = aa6$sdev^2
#   total.var[i,1]= sum(aa6$sdev^2)
#   var.explained[i,1] = ((var.temp[1]+var.temp[2])/total.var[i,1])*100
#   # }
#   
# 
# }









## install BiocManager if not installed if (!requireNamespace("BiocManager", quietly = TRUE))     
install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics')
######################################################
## alternatively: 
install.packages("devtools")
# then load
library(devtools)
install_github("mixOmicsTeam/mixOmics")



detach("package:spls", unload=TRUE)
library(mixOmics)

################################################ spls
N.Comp = 16
rss = matrix(0L, nrow = N.Comp , ncol = N.Comp )
MSPE= matrix(0L, nrow = N.Comp , ncol = N.Comp )

library(MASS)
library(lattice)
library(ggplot2)
library(mixOmics)

X = read.csv('Cost_Components_Eval_Fakuchi_AllData.txt',header=F)
y = read.csv('Cost_Deviation_Eval_y_forPLSR.txt',header=F)


N.Comp = 16




for (i.Comp in seq(1, N.Comp, by=1)){
  print(i.Comp)
  for (i.keep in seq(1, N.Comp, by=1)){
    spls_output = spls(X,y,ncomp = i.Comp,
                       mode = c("regression", "canonical", "invariant", "classic"),
                       keepX=rep(i.keep, i.Comp),
                       scale = TRUE,
                       tol = 1e-06,
                       max.iter = 3000,
                       near.zero.var = FALSE,
                       logratio="none",
                       multilevel=NULL,
                       all.outputs = TRUE)
    
    CV = perf(spls_output, validation = "Mfold", folds = 10, nrepeat = 10, progressBar = FALSE )
    
    # rss[i.keep,i.Comp] = CV$RSS[-(1:i.Comp)]
    # rss[i.keep,i.Comp]  = sum(CV$measures$RSS$summary$mean)
    MSPE[i.keep,i.Comp] = CV$measures$MSEP$summary$mean[i.Comp]
    sprintf('i.keep = %i', i.keep)
  }
}










########
for (i.Comp in seq(1, N.Comp, by=1)){
  print(i.Comp)
  for (i.keep in seq(1, N.Comp, by=1)){
    spls_output = spls(X,y,ncomp = i.Comp,
                       mode = c("regression", "canonical", "invariant", "classic"),
                       keepX=rep(i.keep, i.Comp),
                       scale = TRUE,
                       tol = 1e-06,
                       max.iter = 3000,
                       near.zero.var = FALSE,
                       logratio="none",
                       multilevel=NULL,
                       all.outputs = TRUE)
    
    CV2 = perf(spls_output, validation = "Mfold", folds = 10, nrepeat = 10, progressBar = FALSE )
    
    # rss[i.keep,i.Comp] = CV$RSS[-(1:i.Comp)]
    # rss[i.keep,i.Comp]  = sum(CV$measures$RSS$summary$mean)
    MSPE[i.keep,i.Comp] = CV$measures$MSEP$summary$mean[i.Comp]
    sprintf('i.keep = %i', i.keep)
  }
}
#################







rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}


rand_vect(4, 3)





i.index = 0
keep.vector = matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 16 )
MSPE.vector = matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )
i.Comp.index=matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )
i.keep.index=matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )

########
for (i.Comp in seq(3, N.Comp, by=1)){
  print(i.Comp)
  for (i.keep in seq(1, N.Comp-1, by=1)){
    for (i.rnd in seq(1,50, by=1)){
      i.index= i.index+1
      
      #keep.2comps = rep(i.keep,min(i.Comp,2))
      
      keep.vector[i.index,] = c(rep(i.keep,2), (rand_vect(N.Comp-2,N.Comp-i.keep-1)+1))
      i.Comp.index[i.index,1] = i.Comp
      i.keep.index[i.index,1] = i.keep
      spls_output2 = spls(X,y,ncomp = i.Comp,
                         mode = c("regression", "canonical", "invariant", "classic"),
                         keepX=keep.vector[i.index,],
                         scale = TRUE,
                         tol = 1e-06,
                         max.iter = 3000,
                         near.zero.var = FALSE,
                         logratio="none",
                         multilevel=NULL,
                         all.outputs = TRUE)
      
      CV2 = perf(spls_output2, validation = "Mfold", folds = 5, nrepeat = 5, progressBar = FALSE )
        
      # rss[i.keep,i.Comp] = CV$RSS[-(1:i.Comp)]
      # rss[i.keep,i.Comp]  = sum(CV$measures$RSS$summary$mean)
      #MSPE[i.keep,i.Comp] = CV$measures$MSEP$summary$mean[i.Comp]
      #MSPE.vector[i.index,1] = CV$measures$MSEP$summary$mean[i.Comp]
      MSPE.vector[i.index,1] = sum(CV2$measures$MSEP$summary$mean)
      
      sprintf('i.keep = %i', i.keep)
    }
  }
}
#################







i.index = 0
keep.vector3 = matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 16 )
MSPE.vector3 = matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )
i.Comp.index3=matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )
i.keep.index3=matrix(0L, nrow = N.Comp*(N.Comp-1)*50 , ncol = 1 )
########
for (i.Comp in seq(3, N.Comp, by=1)){
  print(i.Comp)
  for (i.keep in seq(3, N.Comp-1, by=1)){
    for (i.rnd in seq(1,50, by=1)){
      i.index= i.index+1
      
      #keep.2comps = rep(i.keep,min(i.Comp,2))
      
      keep.vector3[i.index,] = c(i.keep, (rand_vect(i.Comp-1,N.Comp-i.keep-(i.Comp-1))+1), matrix(0L, nrow = 1 , ncol = N.Comp-1-(i.Comp-1) ))
      i.Comp.index3[i.index,1] = i.Comp
      i.keep.index3[i.index,1] = i.keep
      spls_output3 = spls(X,y,ncomp = i.Comp,
                          mode = c("regression", "canonical", "invariant", "classic"),
                          keepX=keep.vector3[i.index,1:i.Comp],
                          scale = TRUE,
                          tol = 1e-06,
                          max.iter = 3000,
                          near.zero.var = FALSE,
                          logratio="none",
                          multilevel=NULL,
                          all.outputs = TRUE)
      
      CV3 = perf(spls_output3, validation = "Mfold", folds = 5, nrepeat = 5, progressBar = FALSE )
      
      # rss[i.keep,i.Comp] = CV$RSS[-(1:i.Comp)]
      # rss[i.keep,i.Comp]  = sum(CV$measures$RSS$summary$mean)
      #MSPE[i.keep,i.Comp] = CV$measures$MSEP$summary$mean[i.Comp]
      #MSPE.vector[i.index,1] = CV$measures$MSEP$summary$mean[i.Comp]
      MSPE.vector3[i.index,1] = sum(CV2$measures$MSEP$summary$mean)/i.Comp
      
      sprintf('i.keep = %i', i.keep)
    }
  }
}
#################









MSPE.mean = sum(MSPE)/ (N.Comp*N.Comp/2 + N.Comp/2)

MSPE.mean.sqrd = (MSPE-MSPE.mean)^2
MSPE.mean.sqrd[MSPE.mean.sqrd>1e-2] = 0

MSPE.sd = sqrt( sum(MSPE.mean.sqrd)/(N.Comp*N.Comp/2 + N.Comp/2 -1) )


# ggplot(MSE.SEP.frame, aes(x=components, y = sparsity, fill = SE)) + geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(MSE.SEP.frame[,3]), 
#                                                                                                                      limit = c(min(MSE.SEP.frame[,3]), max(MSE.SEP.frame[,3])), space = "Lab" ) + 
#  labs(x = "Latent Variables", y =  "Sparsity", title = "SEs") +  scale_x_continuous(breaks = c(1:10)) + scale_y_continuous(breaks = c(1:10))

components = c(1:N.Comp)
sparsity   = c(1:N.Comp)


# ggplot(as.data.frame(MSPE), aes(x=components, y = sparsity, fill = SE)) + geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = MSPE.mean, 
#                                                                                                                      limit = c(min(MSPE[1:8,1:8]), max(MSPE[,16])), space = "Lab" ) + 
#  labs(x = "Latent Variables", y =  "Sparsity", title = "SEs") +  scale_x_continuous(breaks = c(1:N.Comp)) + scale_y_continuous(breaks = c(1:N.Comp))
 ggplot(as.data.frame(MSPE), aes(x=components, y = sparsity)) + geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = MSPE.mean, 
                                                                                                                             limit = c(min(MSPE), max(MSPE)), space = "Lab" ) + 
   labs(x = "Latent Variables", y =  "Sparsity", title = "SEs") +  scale_x_continuous(breaks = c(1:N.Comp)) + scale_y_continuous(breaks = c(1:N.Comp))
 
 values <- runif(16, 0, 5)
 
 ggplot(as.data.frame(MSPE), aes(x = components, y = sparsity) , fill=values) + 
   geom_tile(aes(fill = values), color = "white", size = 1) + 
   scale_fill_gradient(low = "gray95", high = "tomato") + 
   xlab("characteristics") + 
   theme_grey(base_size = 20) + 
   ggtitle("Heatmap (ggplot)") + 
   theme(axis.ticks = element_blank(), 
         panel.background = element_blank(), 
         plot.title = element_text(size = 12, colour = "gray50")) 

 
 
 ggplot(as.data.frame(MSPE), aes(x = components, y = sparsity) , fill=values) + 
   geom_tile(color = "white", size = 1) + 
   scale_fill_gradient(low = "gray95", high = "tomato") + 
   xlab("characteristics") + 
   theme_grey(base_size = 16) + 
   ggtitle("Heatmap (ggplot)") + 
   theme(axis.ticks = element_blank(), 
         panel.background = element_blank(), 
         plot.title = element_text(size = 12, colour = "gray50")) 
 
 ggplot(as.data.frame(MSPE), aes(x = components, y = sparsity) , fill=values)
 

i.keep = 7
i.Comp = 6

   spls_output = spls(X,y,ncomp = i.Comp,
                      mode = c("regression", "canonical", "invariant", "classic"),
                      keepX=rep(i.keep, i.Comp),
                      scale = TRUE,
                      tol = 1e-06,
                      max.iter = 3000,
                      near.zero.var = FALSE,
                      logratio="none",
                      multilevel=NULL,
                      all.outputs = TRUE)

View(spls_output$loadings$X)



install.packages("spls")

#####################################
X = read.csv('Cost_Components_Eval_Fakuchi_AllData_NoDynamics.txt',header=F)
y = read.csv('Cost_Deviation_Eval_y_forPLSR.txt',header=F)

X.scaled = scale(X, center = TRUE, scale = TRUE)
y.scaled = scale(y, center = TRUE, scale = TRUE)
X = X.scaled
y = y.scaled


N.Comp = 16



library(spls)

#set.seed(2)
cv = cv.spls(X,y , eta=seq(.1,.9,.02) , K=c(3:N.Comp-1))

f = spls( X , y ,  eta = cv$eta.opt  ,  K = cv$K.opt )
f = spls( X , y ,  eta = cv$eta.opt  ,  K = N.Comp )

print(f)

View(f$projection)

save(MSPE, file = "saveddf.RData")
dput(as.data.frame(MSPE), "saveddf.txt")

MSPE2 <- dget("saveddf.txt")


write.csv(MSPE, file = "saveddf.csv")








## SPCR 

X = read.csv('Cost_Components_Eval_Fakuchi_AllData.txt',header=F)
y = read.csv('Cost_Deviation_Eval_y_forPLSR.txt',header=F)



Data = readMat('Fukuchi_Features_DTW_RMS.mat', fixNames=TRUE, drop=c("singletonLists"),
               sparseMatrixClass=c("Matrix"), verbose=FALSE)



Features.Matrix = Data$Cost.Components.Subjects
DTW.Vector      = Data$DTW.Sum.Subjects
RMSE.Vector     = Data$RMSE.Sum.Subjects





X = Data$Cost.Components
y = Data$DTW.Sum

X.scaled = scale(X, center = TRUE, scale = TRUE)
y.scaled = scale(y, center = TRUE, scale = TRUE)

mean.resid  <- matrix(0L, nrow = numSteps, ncol = 1)


numSteps=10

num.Component=16

for (i in seq(0, numSteps, by=1)){
  sparsity = 5e-2*i/numSteps
  aa6 = robspca(X, k = num.Component, alpha = sparsity, beta = 1e-04, gamma = 100,
                center = TRUE, scale = TRUE, max_iter = 1000, tol = 1e-05, verbose = TRUE)  
  X.Principal = as.matrix(X.scaled)%*%as.matrix(aa6$loadings)
  Results.lm <- lm(y.scaled ~ X.Principal)
  mean.resid[i,1] = sqrt(mean(Results.lm$residuals^2))
}






N.Comp = 16

library(spcr)

lambda.B = 1
lambda.gamma = 1
k = N.Comp

#Results.SPCR = spcr(as.matrix(X), as.vector(y), k, lambda.B, lambda.gamma, w=0.1, xi=0.01, 
#     adaptive=FALSE, center=TRUE, scale=FALSE)



setwd('C:/Users/nozaripo/OneDrive - University of Southern California/OptimTraj-master/demo/fiveLinkBiped - IOC')



library("R.matlab")

Data = readMat('Fukuchi_Features_DTW_RMS.mat', fixNames=TRUE, drop=c("singletonLists"),
               sparseMatrixClass=c("Matrix"), verbose=FALSE)



Features.Matrix = Data$Cost.Components.Subjects
DTW.Vector      = Data$DTW.Sum.Subjects
RMSE.Vector     = Data$RMSE.Sum.Subjects





X = Data$Cost.Components
y = Data$RMSE.Sum

X.scaled = scale(X, center = TRUE, scale = TRUE)
y.scaled = scale(y, center = TRUE, scale = TRUE)


X_scaled = Data$X.scaled
X_scaled = X.scaled


#lambda_values = cv.spcr(as.matrix(X), as.vector(y), k, w=0.1, xi=0.01, nfolds=10, adaptive=FALSE,
#                       center=FALSE, scale=FALSE, lambda.B.length=100, lambda.gamma.length=100,
#                        lambda.B=NULL, lambda.gamma=NULL)

lambda_values = cv.spcr(as.matrix(X_scaled), as.vector(y), k, w=0.1, xi=0.01, nfolds=10, adaptive=FALSE,
                        center=FALSE, scale=FALSE, lambda.B.length=100, lambda.gamma.length=100,
                        lambda.B=NULL, lambda.gamma=NULL)
lambda_values$
lambda.B = 70
lambda.gamma = 1000
#lambda.B = 33.6
#lambda.gamma = 33.6
#Results.SPCR = spcr(as.matrix(X), as.vector(y), k, lambda.B, lambda.gamma, w=0.1, xi=0.01, 
#                    adaptive=FALSE, center=TRUE, scale=TRUE)
Results.SPCR = spcr(as.matrix(X_scaled), as.vector(y), k, lambda.B, lambda.gamma, w=0.1, xi=0.01, 
                    adaptive=FALSE, center=TRUE, scale=FALSE)

View(Results.SPCR$loadings.A)




Variance_Data = matrix(0L, nrow = 1, ncol = N.Comp)
Variance_spcr = matrix(0L, nrow = 1, ncol = N.Comp)
for(i in 1:N.Comp){
  Variance_Data[i] = var(X_scaled[,i])
  Variance_spcr[i] = var(X_scaled%*%Results.SPCR$loadings.A[,i])
  
}

Variance_Data_total = sum(Variance_Data)

VAF_spcr = Variance_spcr/Variance_Data_total


VAF_spcr_sorted = cumsum(sort(VAF_spcr,decreasing=TRUE))

pcr_components = sum(VAF_spcr_sorted<.90)
pcr_components = 8

idx = order(VAF_spcr,decreasing=TRUE)

Loadings_spcr = Results.SPCR$loadings.A[,idx[1:pcr_components]]

as.matrix(VAF_spcr_sorted)

#####
View(Loadings_spcr)
View(as.array(VAF_spcr_sorted))
View(VAF_spcr_sorted)
View(t(as.matrix(VAF_spcr_sorted))[,1:7])


Loadings_spcr.cut = Loadings_spcr
Loadings_spcr.cut[abs(Loadings_spcr.cut)<.3]=0


## Plot the loadings

Criteria.Names = c("ang accel sqr", "ang jerk sqr", "cartes accel sqr", "cartes jerk sqr", 
                   "torqs sqr", "torqs abs", "torqs rate sqr","torqs rate abs", "torq 2nd rate sqr", "torq 2nd rate abs",
                   "yGRF rate sqr", "xGRF rate sqr", "pos work", "abs work",
                   "peak kinetic energy", "pk2pk ang momentum")

par(mfrow=c(1,8),mar=c(8, 10, 2, 2))
i=1
bar.plt <- barplot(Loadings_spcr[,i],
                   main = paste("PC ", i, "\n Cum VAF=", format(VAF_spcr_sorted[i], digits = 2)),
                   ylab = "Criteria",
                   names.arg = Criteria.Names,
                   col = "darkred",
                   las = 2,
                   horiz = TRUE,
                   cex.lab = 1.7,
                   cex.axis = 1.2,
                   cex.names = 1.0)

for(i in seq(2,8,1)){
  
  bar.plt <- barplot(Loadings_spcr[,i],
                     main = paste("PC ", i, "\n Cum VAF=", format(VAF_spcr_sorted[i], digits = 2)),
                     col = "darkred",
                     las = 2,
                     horiz = TRUE,
                     cex.lab = 1.7,
                     cex.axis = 1.2,
                     cex.names = 1.0)
  
  #subplot(bar.plt, nrows = 2, margin = 0.05)
}

mtext("Loadings of criteria on each PC",1)




xlab = "Loadings of criteria on each PC"

i=8
bar.plt <- barplot(matrix(0L, nrow = N.Comp, ncol = 1),
                   main = paste("PC ", i, "\n Cum VAF=", format(VAF_spcr_sorted[i], digits = 2)),
                   ylab = "Criteria",
                   names.arg = Criteria.Names,
                   col = "darkred",
                   las = 2,
                   horiz = TRUE,
                   cex.lab = 1.7,
                   cex.axis = 1.2,
                   cex.names = 1.0)
mtext("Number of links",4)





################## Sparse visualization loading<.30 -> 0
par(mfrow=c(1,8),mar=c(8, 10, 2, 2))
i=1
bar.plt <- barplot(Loadings_spcr.cut[,i],
                   main = paste("PC ", i, "\n Cum VAF=", format(VAF_spcr_sorted[i], digits = 2)),
                   ylab = "Criteria",
                   names.arg = Criteria.Names,
                   col = "darkred",
                   las = 2,
                   horiz = TRUE,
                   cex.lab = 1.7,
                   cex.axis = 1.2,
                   cex.names = 1.0)

for(i in seq(2,8,1)){
  
  bar.plt <- barplot(Loadings_spcr.cut[,i],
                     main = paste("PC ", i, "\n Cum VAF=", format(VAF_spcr_sorted[i], digits = 2)),
                     col = "darkred",
                     las = 2,
                     horiz = TRUE,
                     cex.lab = 1.7,
                     cex.axis = 1.2,
                     cex.names = 1.0)
  
  #subplot(bar.plt, nrows = 2, margin = 0.05)
}

mtext("Loadings of criteria on each PC",1)






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




