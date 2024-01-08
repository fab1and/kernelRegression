#rm(list = ls())

library(np)

#Value Initializing

beta.0 <- 1
beta.2 <- 2
beta.3 <- 3

b.seq     <- seq(from = -0.4, to = 0.4, by = 0.05)

alpha.seq <- c(0.01, 0.05, 0.1)

n.seq <- c(50,100)

S <- 100
B <- 100

#Result Arrays
  I.hat.i.summary <- array(data = NA, dim = c(B,1))
  H0.reject       <- array(data = NA, dim = c(S,length(b.seq)), dimnames = list(NULL, b.seq))

set.seed(42)
(start.time <- Sys.time())
for(a in alpha.seq){
  for (n in n.seq){
    for(beta.1 in b.seq){
      for(s in 1:S){
        x       <- rnorm(n = n, mean = 0, sd = 1)
        z.1     <- sample(x = c(1,0), size = n, replace = TRUE)
        z.2     <- sample(x = c(1,0), size = n, replace = TRUE)
        epsilon <- rnorm(n = n, mean = 0, sd = 1)
        y       <- beta.0 + beta.1*z.1*(1+x^2) + beta.2 * z.2 + beta.3*x + epsilon
        
        dat.sample <- data.frame(y, x, z.1 = as.factor(z.1), z.2 = as.factor(z.2))
        
        np.bw.overall <- npregbw( #using the same cross-validated smoothing parameters for both regressions (Bootstrap method I)
          y = dat.sample[,1],
          x = dat.sample[,-1],
          regtype = "lc", #local constant estimator
          bwmethod = "cv.ls",
          ckertype = "gaussian",
          ckerorder = 2,
          ukertype = "liracine",
          nmulti = 1	
        )
        
        dat.sample.z11 <- dat.sample
        dat.sample.z11[,"z.1"] <- factor(1, levels = levels(dat.sample$z.1))
        
        dat.sample.z10 <- dat.sample
        dat.sample.z10[,"z.1"] <- factor(0, levels = levels(dat.sample$z.1))
        
        m.hat.z11 <- fitted(npreg(
          np.bw.overall,
          txdat		= dat.sample[, -1],
          tydat		= dat.sample[, 1],
          exdat		= dat.sample.z11[, -1]
        ))
        
        m.hat.z10		<- fitted(npreg(
          np.bw.overall,
          txdat		= dat.sample[, -1],
          tydat		= dat.sample[, 1],
          exdat		= dat.sample.z10[, -1]
        ))
        
        I.hat.n <- mean((m.hat.z11 - m.hat.z10)^2)
        
        
        for(i in 1:B){
          
          dat.sample.i <- dat.sample
          
          #difference to above:
          dat.sample.i[,"z.1"] <- sample(x = dat.sample[, "z.1"], size = nrow(dat.sample.i), replace = TRUE)
          
          m.hat.z11.i		<- fitted(npreg(
            np.bw.overall,
            txdat		= dat.sample.i[, -1],
            tydat		= dat.sample.i[, 1],
            exdat		= dat.sample.z11[, -1]
          ))
          
          m.hat.z10.i		<- fitted(npreg(
            np.bw.overall,
            txdat		= dat.sample.i[, -1],
            tydat		= dat.sample.i[, 1],
            exdat		= dat.sample.z10[, -1]
          ))
          
          I.hat.n.i	<- mean((m.hat.z11.i - m.hat.z10.i)^2)
          
          
          I.hat.i.summary[i,] <- I.hat.n.i 
        }
          #Test decision
          H0.reject[s, as.character(beta.1)] <- I.hat.n > quantile(x = I.hat.i.summary, 1-a)
      }
          assign(paste("H0.reject.n", n, ".a", a * 100, sep = ""), H0.reject)
    }
  }
}

### Power curves

X11()
par(mfrow = c(1,2))
plot(b.seq,colSums(H0.reject.n50.a5)/B, type = "l", ylim = c(0,1), main = "alpha = 5%, n = 50", ylab = expression(paste("P(reject  ", H[0], ")", sep = "")))
plot(b.seq,colSums(H0.reject.n100.a5)/B, type = "l", ylim = c(0,1), main = "alpha = 5%, n = 100", ylab = expression(paste("P(reject  ", H[0], ")", sep = "")))
par(mfrow = c(1,1))



########################################################
### Practical Application with the npsigtest function###
########################################################


#install.packages("np")
library(np)

#Data import (please adjust your working directory)

setwd("C:/Users/Fabia/OneDrive/Desktop/Studium/Master/Semester II/Topics in Applied Econometrics")
dat	<- read.table(file = "Pooled.txt", header = TRUE, sep = "\t")

dat.sample <- dat[,c("rent","size","age","heat","warm","bath","green")]


###Calculate the Bandwiths

np.bw <- npregbw( #using the same cross-validated smoothing parameters for both regressions (Bootstrap method I)
  y = dat.sample[,1],
  x = dat.sample[,-1],
  regtype = "lc", #local constant estimator
  bwmethod = "cv.ls",
  ckertype = "gaussian",
  ckerorder = 2,
  ukertype = "liracine",
  nmulit = 1
)

###Use the npsigtest function

np.sig <- npsigtest(np.bw,
                    xdat = dat.sample[,-1],
                    ydat = dat.sample[,1],
                    boot.num = 299,
                    boot.method = "wild",
                    boot.type = "I", #In the Simulation Bootstrap Method I is conducted
                    pivot = TRUE,
                    joint = FALSE,
                    index = seq(1,ncol(dat.sample[,-1])),
                    random.seed = 42)


area_poly <- function(cur, cutoff, side=c(1,-1), col = "grey", border=NA, ...)
{
  if (side[1]>0 )# on the right
  {
    pos <- min(which(cur$x > cutoff))
    end <- length(cur$x)
  }
  else # on the left
  {
    pos <- max(which(cur$x < cutoff))
    end <- 1
  }
  polygon(x=c(cur$x[end:pos], cur$x[pos], cur$x[end]),
          y=c(cur$y[end:pos], 0, 0), col=col, border=border, ...)
}

X11()
par(mfrow = c(2,1))
cc <- density(np.sig$In.bootstrap[,4])
plot(density(np.sig$In.bootstrap[,4]), main = "variable heat with 299 bootstrap samples", xlab = "I.n.bootstrap")
area_poly(cc, cutoff = quantile(np.sig$In.bootstrap[,4],0.99), side = 1, col = "orange", density = 80)
abline(v = np.sig$In[4], lty = 2, col = "black")
legend(x = "topright", legend = c("I.n","99%-quantile"), col = c("black","orange"),lty = c(2,1))


cc <- density(np.sig$In.bootstrap[,5])
plot(density(np.sig$In.bootstrap[,5]), main = "variable bath with 299 bootstrap samples", xlab = "I.n.bootstrap")
area_poly(cc, cutoff = quantile(np.sig$In.bootstrap[,5],0.90), side = 1, col = "orange", density = 80)
abline(v = np.sig$In[5], lty = 2, col = "black")
legend(x = "topright", legend = c("I.n","90%-quantile"), col = c("black","orange"),lty = c(2,1))
par(mfrow = c(1,1))

 

par(mfrow = c(1,2))
boxplot(np.sig$In.bootstrap[,4], main = "variable warm")
points(np.sig$In[4], col = "darkorange", type = "b", pch = 16)
abline(h = quantile(np.sig$In.bootstrap[,4],0.99), lty = 2, col = "red")
legend(x = "bottomleft", legend = c("I.n","99%-quantile","Median"), col = c("darkorange","red","black"), pch = c(16,NA,NA), lty = c(NA,2,1))


boxplot(np.sig$In.bootstrap[,5], main = "variable bath")
points(np.sig$In[5], col = "darkorange", type = "b", pch = 16)
abline(h = quantile(np.sig$In.bootstrap[,5],0.90), lty = 2, col = "red")
legend(x = "bottomleft", legend = c("I.n","90%-quantile","Median"), col = c("darkorange","red","black"), pch = c(16,NA,NA), lty = c(NA,2,1))

par(mfrow = c(1,2))
