require(smooth)  
require(Mcomp)  
library(chron)  
library(RColorBrewer) 
library(lattice) 
library(smoother) 
library(zoo)  #si
library(imputeTS) 
#aditional library
require(RTOMO) 

#Selection   FDP for the calculate of SPI
library(eva) # extrem values analysis
library(MASS) # fit statistics distributions
library(fitdistrplus) # MORES STATISTICAL DISTRIBUTIONS OPTIONS
library(hydroGOF) # for RMSE
library(extRemes) # ALTERNATIVE TO EVA

#Lmoments
library(lmomco) #estimation
library(FAdist) #Distributions that are Sometimes Used in Hydrology
library(dgof)

#ROUTE of THE REGIONS FILES
ruta_carpeta <- '/home/Physics/pro-seq/1_getSPIindex/3_CAL-SPI1/1_calculo_SPI/3meses/est-prcp-65-percent/'

print('%%%%%%%--Inserte la Region correpondintes---- %%%%%%')
region <- readline()
num_reg <- as.integer(gsub('reg','',region))
print('%%%%%%%--En proceso--- %%%%%%')
ruta_carpeta <- paste(ruta_carpeta,region,'/',sep='')
ruta_salida0 <- '/home/Physics/pro-seq/1_getSPIindex/3_CAL-SPI1/1_calculo_SPI/3meses/est-SPI-65-percent/'
ruta_salida <- paste(ruta_salida0,region,sep='')
dir.create(ruta_salida)
ruta_salida <- paste(ruta_salida,'/',sep='')


#AS LIST the elementes in the region
list_est <- dir(ruta_carpeta, pattern = '_prcp.scv')

#null vectors
c_names  <- NULL
c_distr  <- NULL
c_rmse   <- NULL
c_aic    <- NULL
c_kstest <- NULL
c_best   <- NULL
c_numreg    <- NULL
#LOOP in lIST (list_est) ###
for (estation in list_est) {
  ruta_est <- paste(ruta_carpeta,estation,sep = '')
  name_est <- gsub('_prcp.scv','',estation)
  df_prec <- read.table(ruta_est, sep ='\t')
  prec <- df_prec$V2
  times <- c(1:length(prec)) 
     
  ######MOving Averange ########### 
  nfil = 3 #STEP
  #################################
  if (nfil  == 3) {
    k=1 #3-month low-pass
    phantomVector <- NaN*numeric(k)
    precPhantom <-c(phantomVector,prec,phantomVector)
    filtPrec<-numeric(length(prec))
    for (i in seq((1+k),length(precPhantom)-k)){ 
      filtPrec[i-k]= mean(precPhantom[(i-k):(i+k)],na.rm = TRUE)}
  } else if (nfil == 6) {
    filtPrec = filter(prec, sides=2, filter=rep(1/nfil,nfil)) # moving average sides = 2 centered
    filtPrec = filtPrec[3:417]
    init = c(NA,NA)
    fin =  c(NA,NA,NA)
  } else if (nfil == 12) {
    filtPrec = filter(prec, sides=2, filter=rep(1/nfil,nfil)) # moving average sides = 2 centered
    filtPrec = filtPrec[6:414]
    init = c(NA,NA ,NA , NA, NA)
    fin =  c(NA,NA,NA, NA, NA,NA)
  }


  ###################################################################
  ########################calculate SPI##############################
  ###################################################################
   
  #FUNCTION TO CALCULATE SPI
  calSPI <- function(filtPrec , fitCDF,nonZeros, ProbZeros){
    Hx<-numeric(length(filtPrec))
    Hx[which(!filtPrec == 0)] = ProbZeros + (1-ProbZeros)*fitCDF
    Hx[which(filtPrec == 0)]  = ProbZeros

    #SPI VECTOR CALC SPI
    spiFin <- numeric(length(filtPrec))
    Hx05 <- Hx[((Hx > 0.0) & (Hx <= 0.5))]
    Hx01 <- Hx[((Hx > 0.5) & (Hx <= 1.0))]
    #Calculate SPI 
    #Constantes
    c <- c(2.515517, 0.802853, 0.010328)
    d <- c(1.432788, 0.189269, 0.001308)

    t5 <- sqrt(log(1/(Hx05)**2))
    t1 <- sqrt(log((1/(1.0-Hx01)**2)))

    spiFin[which((Hx > 0.0) & (Hx <= 0.5))]   = - (t5-(c[1]+c[2]*t5+c[3]*t5**2)/(1+d[1]*t5+d[2]*t5**2+d[3]*t5**3))
    spiFin[which((Hx > 0.5) & (Hx <= 1.0))]   = + (t1-(c[1]+c[2]*t1+c[3]*t1**2)/(1+d[1]*t1+d[2]*t1**2+d[3]*t1**3))
    return(spiFin)}


  ############## first find zeros and non-zeros in prcp series ######
  nonZeros<-filtPrec[which(!filtPrec == 0)]
  zeros<-filtPrec[which(filtPrec == 0)]
  ProbZeros<-length(zeros)/length(filtPrec)

  #FIT Gamma , Weibull, normal and Pearson type3.  calculate CDF by one
  fitgamma       <- fitdistr(matrix(nonZeros),"gamma") 
  fitgamma.cdf   <- pgamma(nonZeros, shape = fitgamma$estimate[1], rate = fitgamma$estimate[2], lower.tail = TRUE, log.p = FALSE)

  fitweibull     <- fitdistr(matrix(nonZeros),"weibull")
  fitweibull.cdf <- pweibull(nonZeros, shape = fitweibull$estimate[1], scale = fitweibull$estimate[2], lower.tail = TRUE, log.p = FALSE)

  fitnormal <-      fitdistr(matrix(nonZeros),"normal")
  fitnormal.cdf<-   pnorm(nonZeros, mean = fitnormal$estimate[1], sd = fitnormal$estimate[2], lower.tail = TRUE, log.p = FALSE)

  fitlognormal <-   fitdistr(matrix(nonZeros),"log-normal")
  fitlognormal.cdf<-plnorm(nonZeros, meanlog = fitlognormal$estimate[1], sdlog = fitlognormal$estimate[2], lower.tail = TRUE, log.p = FALSE)

  lmr <- lmoms(nonZeros)
  fitpear3 = parpe3(lmr)
  #pdf3lm <- pdfpe3(nonZeros,pe3) 
  fitpear3.cdf <- cdfpe3(nonZeros,fitpear3)

  #primera opción in equal
  MU    <- fitpear3$para[1] # product moment mean
  SIGMA <- fitpear3$para[2] # product moment standard deviation
  GAMMA <- fitpear3$para[3] # product moment skew
  L <- fitpear3$para[1] - 2*SIGMA/GAMMA # location
  S <- (1/2)*SIGMA*abs(GAMMA)       # scale
  A <- 4/GAMMA^2                    # shape
  
  params <-c(A,S,L)
  #prob <- pgamma3(nonZeros,shape=A,scale=S,thres=L,lower.tail=TRUE,log.p=FALSE)

  #method MLE
  #mle <- mle2par(nonZeros, type = 'pe3')
  #pdfP3mle <- pdfpe3(nonZeros, mle)
  #cdfP3mle <- cdfpe3(nonZeros, mle)
  

  ######################### Which one is better?   #########################
  ## METHOD 1: COMPARISON OF EMPIRICAL CDF AGAINST FITTED CDF
  Fn<-ecdf(filtPrec) # Empirical CDF en filprec argmax(filtprec)
  sort(nonZeros)
 
  RMSE =c(numeric(5))
  RMSE[1] <- rmse(Fn(nonZeros),fitgamma.cdf)
  RMSE[2] <- rmse(Fn(nonZeros),fitweibull.cdf)
  RMSE[3] <- rmse(Fn(nonZeros),fitnormal.cdf)
  RMSE[4] <- rmse(Fn(nonZeros),fitlognormal.cdf)
  RMSE[5] <- rmse(Fn(nonZeros),fitpear3.cdf )

  #METHOD 2 : AKAIKE INFORMATION CRITERIUM
  #log likelihood for Gamma(alpha,scale=beta), X Data input
  LLabGamma <- function(pars,X){
    alpha <- pars[1] #parametro 1
    beta  <- pars[2]  #parametro 2
    return(-2*sum(log(dgamma(X,alpha,rate=beta))) + 2*2 )
   }

  LLabWeibull <- function(pars,X){
    alpha <- pars[1]
    beta <- pars[2]
    return(-2*sum(log(dweibull(X,alpha,scale=beta))) + 2*2)
  }

  LLabNormal <- function(pars,X){
    media <- fitnormal$estimate[1]
    dstd  <- fitnormal$estimate[2]
    return(-2*sum(log(dnorm(X,mean = media, sd = dstd))) + 2*2)
  }

  LLablogNormal <- function(pars,X){
    media <- fitlognormal$estimate[1]
    dstd  <- fitlognormal$estimate[2]
    return(-2*sum(log(dlnorm(X,meanlog = media, sdlog = dstd))) + 2*2)
  }

  LLabpear3 <- function(pars , X){
    return(-2*sum(log(dgamma3(X,shape=params[1],scale=params[2],thres=params[3],log=FALSE)) )+3*2)
  }


  #Call FUNCTIONS---------------------------------------
  AIC =c(numeric(5))
  AIC[1]    <- LLabGamma(fitgamma$estimate,    nonZeros)                 
  AIC[2]    <- LLabWeibull(fitweibull$estimate,nonZeros) 
  AIC[3]    <- LLabNormal(fitnormal$estimate,  nonZeros) 
  AIC[4]  <- LLablogNormal(fitlognormal$estimate, nonZeros)
  AIC[5]      <- LLabpear3(params, nonZeros)

  ###K-S TEST #####
  Alternativo <- c('two.sided')
  Dgamma    <- ks.test(x = sort(Fn(nonZeros)), y = sort(fitgamma.cdf),alternative=Alternativo, exact = FALSE)
  Dweibull  <- ks.test(x = sort(Fn(nonZeros)), y = sort(fitweibull.cdf), alternative=Alternativo, exact = FALSE)
  Dnormal   <- ks.test(x = sort(Fn(nonZeros)), y = sort(fitnormal.cdf), alternative=Alternativo, exact = FALSE)
  Dlgnormal <- ks.test(x = sort(Fn(nonZeros)), y = sort(fitlognormal.cdf), alternative=Alternativo, exact = NULL)
  Dpear3    <- ks.test(x = sort(Fn(nonZeros)), y = sort(fitpear3.cdf), alternative=Alternativo, exact = NULL)
#
  KSTEST <- c(Dgamma$statistic[[1]], Dweibull$statistic[[1]], Dnormal$statistic[[1]], Dlgnormal$statistic[[1]] , Dpear3$statistic[[1]])
  KSTEST <- abs(KSTEST)
  #MINIMUN VALUES
  minRMSE<- which.min(RMSE)
  minKS <-  which.min(KSTEST)
  minAIC <- which.min(AIC)

 #select the best
 indexSelected <- which.min(abs((RMSE-RMSE[minRMSE]) /RMSE[minRMSE])+ abs((KSTEST-KSTEST[minKS])/KSTEST[minKS]) +  abs((AIC-AIC[minAIC]) /AIC[minAIC])  )


  #SECUENCIA EN X
  x <- seq(0.00001, 400.0, by = 0.02)
  #lista legend
  legends <- c('Gamma','Weibull','Normal', 'log-Normal' ,'Pearson III')

  #SE LECCIONAMOS EL MEJOR INDICE #############################################################################################
  if (indexSelected==1){
     fitDistr <- fitgamma 
     idajuste <- 'Gamma'
     fitCDF <- fitgamma.cdf  
     legends[1] <- paste( idajuste,'(*)', sep=' ')
     BEST_FIT = c('YES','NO','NO','NO', 'NO')
     print('ES GAMMA!!!')

  } else if (indexSelected==2) {
     fitDistr <- fitweibull
     fitCDF <-fitweibull.cdf 
     idajuste <- 'Weibull'
     legends[2] <- paste( idajuste,'(*)', sep=' ')
     BEST_FIT = c('NO','YES','NO','NO', 'NO')
     print('CDF WEIBULL')
  } else if (indexSelected==3) {
     fitDistr <- fitnormal
     fitCDF <-fitnormal.cdf
     idajuste<- 'Normal'
     legends[3] <- paste( idajuste,'(*)', sep=' ')
     BEST_FIT = c('NO','NO','YES','NO', 'NO')
     print('ES NORMAL!!!')

  } else if (indexSelected==4) {
     fitDistr <- fitlognormal
     fitCDF <-fitlognormal.cdf
     idajuste<- 'log-Normal'
     legends[4] <- paste( idajuste,'(*)', sep=' ')
     BEST_FIT = c('NO','NO','NO','YES', 'NO')
     print('Es logNormal')

  } else if (indexSelected==5) {
     fitDistr <- params
     fitCDF <-fitpear3.cdf 
     idajuste<- 'Pearson III'
     legends[5] <- paste( idajuste,'(*)', sep=' ')
     BEST_FIT = c('NO','NO','NO','NO','YES')
     print('Es Pearson III')
  }
    

  parms_out <- capture.output(fitDistr)
  #PRINT best fits in txt
  cat(paste("Es Mejor ajuste es-------------:", idajuste, sep =' '), parms_out, ' ' ,  file=paste(ruta_salida,name_est, '_summary.txt', sep = ''), sep="\n", append=FALSE, fill = TRUE)
   #cat(' ' ,paste('Resultados de copula--------------:',name_cop, sep =''), cop_out , file=paste(ruta_salida,name_est, '_summary.txt', sep = ''), sep="\n", append=TRUE, fill = TRUE)

  distrs = c('Gamma','Weibull','Normal', 'lgNormal', 'Pearson III')
  ##ALmacenado valores
  c_names  <- append(c_names, rep(name_est, times = 5))
  c_distr  <- append(c_distr, distrs)
  c_rmse   <- append(c_rmse, RMSE)  
  c_aic    <- append(c_aic, AIC)   
  c_kstest <- append(c_kstest, KSTEST)
  c_best   <- append(c_best, BEST_FIT)    
  c_numreg    <- append(c_numreg, rep(num_reg, times=5))  

  #PLOT HISTOGRAM AND BEST-FIT DISTRIBUTION 
  pdf(paste(ruta_salida, name_est,'_plots.pdf', sep = ''), width=9.0, height=4.1)
  par(mfrow = c(1, 2))
  hist(nonZeros,ylab="Densidad de Probabilidad",xlab="Precipitación [mm]", main = '',  breaks = 12, prob = TRUE)
  #hist(pweibull, prob=TRUE, breaks = 30, ylim = c(0,270))
  curve(dgamma(x, shape = fitgamma$estimate[1], rate = fitgamma$estimate[2], log=FALSE),add=TRUE, col = 'black', lwd = 1.8, lty = 4)
  curve(dweibull(x, shape = fitweibull$estimate[1], scale = fitweibull$estimate[2], log=FALSE),add=TRUE, lwd = 1.8, col = 'black', lty = 5)
  curve(dnorm(x,  mean = fitnormal$estimate[1], sd = fitnormal$estimate[2], log = FALSE), add = TRUE,col = 'black', lwd = 1.8, lty = 3)
  curve(dlnorm(x, meanlog  = fitlognormal$estimate[1], sdlog  = fitlognormal$estimate[2], log = FALSE), add = TRUE,col = 'black', lwd = 1.8, lty = 1)
  curve(dgamma3(x,shape=params[1],scale=params[2],thres=params[3],log=FALSE), add = TRUE,col = 'black', lwd = 1.8, lty = 2)
  
  legend('topright', legend=legends, col= c('black','black', 'black','black', 'black'), lty = c(4,5,3,1,2),
        lwd =   c(1.8,1.8,1.8,1.8, 1.8) ,box.lty=0, cex = 0.8)

 
  #PLOT comulative ECD Y CDF comparation ###########################################################################################
  plot(sort(nonZeros),sort(Fn(nonZeros)),type = 's' , ylab = 'Probabilidad Acumulada', xlab = 'Precipitación[mm]', cex =0.5, col = 'purple4', lwd=1.8)
  lines(sort(nonZeros), sort(fitgamma.cdf), col = 'black', lwd = 1.3, lty = 5)
  lines(sort(nonZeros), sort(fitweibull.cdf), col = 'black', type  = 'l', lty = 2, lwd = 1.3)
  lines(sort(nonZeros), sort(fitnormal.cdf), col = 'black', type  = 'l', lty = 4, lwd = 1.3)
  lines(sort(nonZeros), sort(fitlognormal.cdf), col = 'black', type  = 'l', lty = 6, lwd = 1.3)
  lines(sort(nonZeros), sort(fitpear3.cdf ), col = 'black', type  = 'l', lty = 3, lwd = 1.7)
  #legend('bottomright', legend = '   CDF Empirico', col = 'black', pch = 4,  ,box.lty=0, cex = 0.8)
  legend('bottomright', legend=c(' Empirico',legends ), col= c('purple4','black', 'black','black', 'black', 'black'), lty = c(1,5,2,4,6,3), 
          lwd = c(1.9,1.2,1.3,1.3, 1.3, 1.8) ,box.lty=0,  cex = 0.8)
  #legend('topleft', legend= paste(name_est,': Mejor Ajuste ',idajuste,sep =''), col= 'black', cex = 0.75,box.lty=0)
  dev.off()

  #CALCULATE SPI  WITH pOLINOMIAL APPROXIMATION
  SPI <- calSPI(filtPrec , fitCDF,nonZeros, ProbZeros)
  if (nfil  == 3) {
    spiComplete = SPI
  }else{
    spiComplete <- c(init, spi,fin)  
  } 
  
  fechas <-seq(as.Date('1980/1/1'), as.Date('2014/12/31'), 'month')
  #SAVE SPICOMPLETE
  df_spi <- data.frame(fechas, spiComplete)
  write.table(df_spi, file = paste(ruta_salida,name_est,'_spi.scv', sep = ''), row.names = FALSE, col.names = FALSE, sep = '\t')
  rm(df_spi)

  #LIMits
  max_spi <- max(spiComplete) + 0.2
  min_spi <- min(spiComplete) - 0.1
  ranges <- c(min_spi, max_spi)
  x <- seq(1980.0, 2015.0,by = 0.0834)

  #PLot SPI COMPLETE
  pdf(paste(ruta_salida,name_est,'_SPI.pdf', sep = ''), width=9.0, height=4.0)
  plot(x,spiComplete,ylim = ranges ,type='h',col='black', lwd = 1.5, xlab = '', ylab = 'SPI',xaxt = 'n')
  #lines(fechas,c(init,spiGamma2,fin),type='h',col='red', lwd = 0.7, xlab = '', ylab = 'SPI', xaxt = 'n')
  #axis.Date(1, at  = fechas, cex.axis=1, format = '%Y')
  axis(1, at = seq(1980, 2014, by = 2), las=2)
  #legend('bottomright', legend='SPI(3 meses)', col='black', box.lty=0,cex = 0.75)
  legend('topright', legend= paste( 'SPI(3 meses)','de la estación :',name_est, sep=' '), col='black', box.lty=0,cex = 0.75)
  par
  abline(h = 0, untf = FALSE)
  abline(h = 1, untf = FALSE, lwd=0.4, lty=2)
  abline(h = -1, untf = FALSE, lwd=0.4, lty=2)
  dev.off()
    
  }


#CREATE DATA FRAME FITS and save #################################################################################################
df <- data.frame(c_numreg,c_names, c_distr, c_rmse, c_aic, c_kstest, c_best)
write.table(df,file  = paste(ruta_salida0,region,'_allNfo.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)

as <- df[df[[7]]=='YES', ]
write.table(as,file  = paste(ruta_salida0,region,'_bestNfo.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)
print('###############################################')
print('#              DONE RIGHT ThinkPad            #')
print('###############################################')

 #chigamma    <- chisq.test(x= sort(Fn(nonZeros)), y =  sort(fitgamma.cdf),  correct = TRUE, rescale.p=TRUE)
 #chiweibull  <- chisq.test(x = sort(Fn(nonZeros)),p = sort(fitweibull.cdf),  correct = TRUE, rescale.p=TRUE)
 #chinormal   <- chisq.test(x =sort(Fn(nonZeros)), p= sort(fitnormal.cdf),   correct = TRUE, rescale.p=TRUE)
 #chilgnormal <- chisq.test(x = sort(Fn(nonZeros)), p=sort(fitlognormal.cdf),  correct = TRUE, rescale.p=TRUE)
 #chipear3    <- chisq.test(x = sort(Fn(nonZeros)),p=  sort(fitpear3.cdf),   correct = TRUE, rescale.p=TRUE)



