simulateThresholdAutoRegressiveTS <- function(length=500, m=2, sigma=0.2, noRegimes=1, gamma=0, th=0, omega=0, phi=c(0.8, -0.5, 0.3)) {
  
  e <- sigma * rnorm(length+501)
  y <- rnorm(2)
  
  f <- data.frame()
  
  for(t in 3:(length+500)) {
    
    if(noRegimes == 2){
      f[t,1] <- 1 / (1 + exp(- gamma[1] *
                    (omega[1] * y[t-1] + omega[2] * y[t-2] - th[1])))
      
    } else if(noRegimes > 2) {
      for(reg in seq(length=(noRegimes-1)))  {
        f[t,reg] <- 1 / (1 + exp(- gamma[reg] *
                      (omega[1,reg] * y[t-1] + omega[2,reg] * y[t-2] - th[reg])))
        
      }      
    }
    
    x <- c(1, y[(t-1):(t-2)])
    
    if(noRegimes > 1) {
      
      y[t] <- phi[1,] %*% x + e[t]
      
      for(reg in seq(length=(noRegimes-1)))  {
        y[t] <- y[t]+ phi[(reg+1),] %*% x * f[t,reg]      
      }
      
    } else {
      y[t] <- phi %*% x + e[t]
    }
    
  }
  
  y <- y[501:(length+500)]  
  
  return(y)
}


########################################################
########################################################

runNcstarTest <- function(T=500,m=2,sigma=0.2,noRegimes=1,gamma=0,th=0,omega=0,phi=c(0.8, -0.5, 0.3), 
    datasize=3, svIter=100, alg="BFGS", noRegimesEstimate = Inf, cluster=c("localhost","localhost"), verbosity=1) {
  
library(tsDyn)
ncstarModels <- array(list(), 2*datasize)
nR <- array(NA, datasize)
nR_acc <- NA

phi1_acc <- NULL
phi2_acc <- NULL
phi1_median_acc <- NULL
phi2_median_acc <- NULL

plot <- FALSE
if(verbosity > 0) plot <- TRUE

for(i in 1:datasize) {
  cat("\n\n-------------------------------------------------------------\n",
      "DATASET NUMBER ", i, 
      "\n-------------------------------------------------------------\n")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # BUILD UP A TIME SERIES
  
  cat("Building time series with ", noRegimes, "regimes...\n")
  
  set.seed(i)
  y <- simulateThresholdAutoRegressiveTS(T,m,sigma,noRegimes,gamma,th,omega,phi)
  #plot(y,type="l")

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # ESTIMATE THE MODEL
  
  ncstarModels[i] <-
      try(list(ncstar(y, m, noRegimes = noRegimesEstimate, alg = alg,
                  cluster = cluster, svIter = svIter, trace = TRUE)));
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # COMPUTE MEDIANS
  
  nR[i] <- ncstarModels[[i]]$noRegimes
  nR_acc[i] <- mean(nR[1:i])
  cat("\n*** Mean of the number of regimes so far: ", mean(nR[1:i]))
  
  if(ncstarModels[[i]]$model.specific$noRegimes == noRegimes) { 
    
      cat("\n*** Median of the linear parameters:\n    ")
      phi1_acc <- rbind(phi1_acc,
          as.vector(ncstarModels[[i]]$model.specific$phi1))
      phi1_median <- apply(phi1_acc, 2, median)
      dim(phi1_median) <- c(noRegimes, m+1)
      print(phi1_median)
      phi1_median_acc <- rbind(phi1_median_acc, phi1_median)
      
      cat("\n*** Median of the nonlinear parameters:\n")
      if(noRegimes < 2) {
        cat("only one regime, so no nonlinear parameters")
      } else {
        phi2_acc <- rbind(phi2_acc,
            as.vector(ncstarModels[[i]]$model.specific$phi2))
        phi2_median <- apply(phi2_acc, 2, median) 
        phi2_median_acc <- rbind(phi2_median_acc, phi2_median)

        cat("gamma = ", phi2_median[1:(noRegimes - 1)])
        cat("; th = ", phi2_median[noRegimes:(2*(noRegimes - 1))])
        
        omega_median <- phi2_median[(2 * (noRegimes - 1) + 1):length(phi2_median)]
        dim(omega_median) <- c(m, noRegimes - 1)
        cat("\nomega = ", omega_median)
        
        if(plot) {
          par(mfrow=c((noRegimes-1), 2), pty = "m")
          
          for(reg in seq(length=2*(noRegimes-1)))  {
            
            if(reg < noRegimes) {
              yLabel <- paste("Gamma_",reg,sep="")
              plot(phi2_acc[,reg], type="l", ylab=yLabel);
              lines(phi2_median_acc[,reg], lty=2)
              abline(h=gamma[reg], col="red")
              
            } else {
              yLabel <- paste("Th_",reg-noRegimes+1,sep="")   
              plot(phi2_acc[,reg], type="l", ylab=yLabel);
              lines(phi2_median_acc[,reg], lty=2)
              abline(h=th[reg-noRegimes+1], col="red")
              
            }

            
          }
        }
      }
  }
}

#ncstarModels <- array(list(), 2*datasize)
#nR <- array(NA, datasize)
#nR_acc <- NA
#
#phi1_acc <- NULL
#phi2_acc <- NULL
#phi1_median_acc <- NULL
#phi2_median_acc <- NULL
#return list()

}


###########################################################
# EXAMPLE 1
###########################################################

datasize <- 3; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 100; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
cluster <- NULL

noRegimes <- 1
m <- 2 
sigma <- 1; # Varianza del modelo

phi <- c(0.8, -0.5, 0.3)

noRegimesEstimate <- Inf
verbosity <- 1

runNcstarTest(T,m,sigma,noRegimes,gamma=0,th=0,omega=0,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

###########################################################
# EXAMPLE 2
###########################################################
library(minpack.lm)

datasize <- 3; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 10; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "LM"
cluster <- NULL

noRegimes <- 2
m <- 2 
sigma <- 0.02; # Varianza del modelo

gamma <- c(20) 
th <- c(0.02)
omega <- c(1, 0)
phi <- rbind(c(0.0, 1.8, -1.06),
    c(0.02, -0.9, 0.795))

noRegimesEstimate <- Inf
verbosity <- 1

runNcstarTest(T,m,sigma,noRegimes,gamma,th,omega,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

###########################################################
# EXAMPLE 3
###########################################################

datasize <- 3; # Número de instancias del modelo que serán generadas
svIter <- 100; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
#cluster <- NULL
cluster <- c("localhost","localhost")
library(snow)

T <- 500; # Tamaño de cada instancia
noRegimes <- 3 
m <- 2 
sigma <- 0.5; # Varianza del modelo

gamma <- c(20, 20) 
th <- c(-0.6, 0.6)
omega <- cbind(c(1, 0), c(1, 0))
phi <- rbind(c(-0.1, 0.3, 0.2),
    c(0.0, -1.2, 0.5),
    c(0.0, 1.8, -1.2))

noRegimesEstimate <- Inf
verbosity <- 1

runNcstarTest(T,m,sigma,noRegimes,gamma,th,omega,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

###########################################################
# EXAMPLE 4
###########################################################

datasize <- 3; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 1000; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
cluster <- c("localhost","localhost")
library(snow)

# PARÁMETROS DEL MODELO ORIGINAL
noRegimes <- 2
m <- 2 
sigma <- 0.5; # Varianza del modelo

gamma <- c(11.31) 
th <- c(0.1414)
omega <- c(0.7071, -0.7071)
phi <- rbind(c(0.5, 0.8, -0.2),
    c(-0.5, -1.2, 0.8))

noRegimesEstimate <- Inf
verbosity <- 1

runNcstarTest(T,m,sigma,noRegimes,gamma,th,omega,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

###########################################################
# EXAMPLE 5
###########################################################

datasize <- 3; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 1000; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "BFGS"
cluster <- c("localhost","localhost")
library(snow)

# PARÁMETROS DEL MODELO ORIGINAL
noRegimes <- 3 
m <- 2 
sigma <- 1; # Varianza del modelo

gamma <- c(8.49, 8.49) 
th <- c(-1.0607, 1.0607)
omega <- cbind(c(0.7071, -0.7071), c(0.7071,-0.7071))
phi <- rbind(c(0.5, 0.8, -0.2),
    c(1.5, -0.6, -0.3),
    c(-0.5, -1.2, 0.7))

noRegimesEstimate <- noRegimes
verbosity <- 1

runNcstarTest(T,m,sigma,noRegimes,gamma,th,omega,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

###########################################################
# EXAMPLE 6
###########################################################

library(rgenoud)

datasize <- 1000; # Número de instancias del modelo que serán generadas
T <- 500; # Tamaño de cada instancia
svIter <- 1000; # Número de iteraciones de la búsqueda de p. iniciales

alg <- "GAD"
#cluster <- NULL
cluster <- c("localhost","localhost")
library(snow)

# PARÁMETROS DEL MODELO ORIGINAL
noRegimes <- 5
m <- 2 
sigma <- 0.2; # Varianza del modelo

gamma <- c(8.49, 8.49, 14.23, 14.23) 
th <- c(-1.0607, -0.59, 0.59, 1.0607)
omega <- cbind(c(0.7071, -0.7071), c(0.7071,-0.7071),
    c(0.7071,-0.7071), c(0.7071,-0.7071))
phi <- rbind(c(0.5, 0.8, -0.2),
    c(1.5, -0.6, -0.3),
    c(0.2, 0.3, -0.9),
    c(-1.2, 0.6, 0.8),
    c(-0.5, -1.2, 0.7))

noRegimesEstimate <- Inf
verbosity <- 0

#runNcstarTest(T,m,sigma,noRegimes,gamma,th,omega,phi,datasize,svIter,alg,noRegimesEstimate,cluster,verbosity)

