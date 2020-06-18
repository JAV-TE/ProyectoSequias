require(smooth)  #si
require(Mcomp)  #si
library(chron)  #si
library(RColorBrewer) #si
library(lattice) #si
#library(ncdf4)  #si
library(smoother) #si
library(zoo)  #si
library(imputeTS) #si univariante series inputation
library(copula) #si
#aditional library
require(RTOMO) #si

#Seleccion  de la FDP para el calculo de SPI
library(eva) # analisis de valores extremos
library(MASS) # ajuste de distribuciones estadisticas
library(fitdistrplus) # MORES STATISTICAL DISTRIBUTIONS OPTIONS
library(hydroGOF) # for RMSE
library(extRemes) # ALTERNATIVE TO EVA

###open file###
#RUTA DE LA CARPETA DE ESTACIONES


#ruta_carpeta  <- paste('/home/Physics/TESIS/Agrupacion/Estaciones_Agrupadas/',region,'/',sep='')
ruta_carpeta <- '/home/Physics/TESIS/Agrupacion/ests_reagroup/'
print('%%%%%%%--Inserte la Region correpondintes---- %%%%%%')
#region <- readline()
region <- ''
print('%%%%%%%--En proceso--- %%%%%%')
#ruta_carpeta <- paste(ruta_carpeta,region,'/',sep='')
ruta_carpeta <- '//home/Physics/TESIS/Agrupacion/medias_cluster-65/'

ruta_salida <- '/home/Physics/pro-seq/1_getSPIindex/3_CAL-SPI1/1_calculo_SPI/3meses/SPI-al-65-percent/'
ruta_salida <- paste(ruta_salida,region,sep='')
#dir.create(ruta_salida)
#ruta_salida <- paste(ruta_salida,'/',sep='')

list_est <- dir(ruta_carpeta, pattern = '_only.scv')

#VeCTORES VACIOS PARA ALMACENAR
c_names  <- NULL
c_distr  <- NULL
c_rmse   <- NULL
c_aic    <- NULL
c_kstest <- NULL
c_best   <- NULL

#LOOP in lIST (list_est) #########################################################################
for (estation in list_est) {
     ruta_est <- paste(ruta_carpeta,estation,sep = '')
     name_est <- gsub('_only.scv','',estation)
     df_prec <- read.table(ruta_est, sep ='\t')
     prec <- df_prec$V2
     times <- c(1:length(prec)) 
     
     ##############################MOving Averange ###############################################    
     nfil = 3 #STEP
     #############################################################################################
     if (nfil  == 3) {
         k=1 #3-month low-pass
         phantomVector <- NaN*numeric(k)
         precPhantom <-c(phantomVector,prec,phantomVector)
         filtPrec<-numeric(length(prec))
         for (i in seq((1+k),length(precPhantom)-k)){ 
           filtPrec[i-k]= mean(precPhantom[(i-k):(i+k)],na.rm = TRUE)
           }
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

     ###############################################################################################
     ########################################Calculo del SPI########################################
     ###############################################################################################
     ################## Primero encuentra ceros y no ceros en las series temporales ################
     nonZeros<-filtPrec[which(!filtPrec == 0)]
     zeros<-filtPrec[which(filtPrec == 0)]
     ProbZeros<-length(zeros)/length(filtPrec)

     # Una histograma exploratorio
     pdf(paste(ruta_salida, name_est,'_hist.pdf', sep = ''), width=5.0, height=4.5)
     hist(filtPrec,ylab="Frecuencia",xlab="Precitación [mm]", col = 'Gray', main = '',  breaks = 15)
     dev.off()

     ################################FIT GAMM AND WEIBULL PDF#########################################
     ##Ajustar a distribucion gamma---------------------------------------------------
     fitgamma <- fitdistr(matrix(nonZeros),"gamma") #Maximum-likelihood fit, pgamma()vector de probabilities
     fitgamma.cdf<- ProbZeros + (1-ProbZeros)*pgamma(nonZeros, shape = fitgamma$estimate[1], rate = fitgamma$estimate[2], lower.tail = TRUE, log.p = FALSE)

     ## Calculate SPI index if Gamma-------------------------------------------------------------------------------
     HxGamma<-numeric(length(filtPrec))
     #calculamos H(x) en un array de distribucion
     HxGamma[which(!filtPrec == 0)] = fitgamma.cdf
     HxGamma[which(filtPrec == 0)]  = ProbZeros

     #SPI VECTOR CALC SPI
     spiGamma <- numeric(length(filtPrec))
     Hx05 <- HxGamma[((HxGamma > 0.0) & (HxGamma <= 0.5))]
     Hx01 <- HxGamma[((HxGamma > 0.5) & (HxGamma <= 1.0))]

     #Calculate SPI  with GAMMA distribution Gma METHOD 1 ----------------------------
     #Constantes
     c <- c(2.515517, 0.802853, 0.010328)
     d <- c(1.432788, 0.189269, 0.001308)
     t5 <- sqrt(log(1/(Hx05)**2))
     t1 <- sqrt(log((1/(1-Hx01)**2)))

     spiGamma[which((HxGamma > 0.0) & (HxGamma <= 0.5))]   = - (t5-(c[1]+c[2]*t5+c[3]*t5**2)/(1+d[1]*t5+d[2]*t5**2+d[3]*t5**3))
     spiGamma[which((HxGamma > 0.5) & (HxGamma <= 1.0))] = + (t1-(c[1]+c[2]*t1+c[3]*t1**2)/(1+d[1]*t1+d[2]*t1**2+d[3]*t1**3))

     #SPI Gma METHOD 2 ----------------------------------------------------------------
     fitgamma.cdf2<- (1-ProbZeros)*pgamma(nonZeros, shape = fitgamma$estimate[1], rate = fitgamma$estimate[2], lower.tail = TRUE, log.p = FALSE)  
     ## Calculate SPI index if Gamma
     spiGamma2 <- numeric(length(filtPrec))
     spiGamma2[which(!filtPrec == 0)]=(fitgamma.cdf2-mean(fitgamma.cdf2))/sd(fitgamma.cdf2)
     spiGamma2[which(filtPrec == 0)]=(0-mean(fitgamma.cdf2))/sd(fitgamma.cdf2)
  


    #Calculate SPI  with WEIBULL distribution  --------------------------------------
    ## Ajustar con la distribucion Weibull
    fitweibull <- fitdistr(matrix(nonZeros),"weibull")
    fitweibull.cdf<- (1-ProbZeros)*pweibull(nonZeros, shape = fitweibull$estimate[1], scale = fitweibull$estimate[2], lower.tail = TRUE, log.p = FALSE)

    ## Calculate SPI index if Weibull
    spiWeibull<-numeric(length(filtPrec))
    spiWeibull[which(!filtPrec == 0)]=(fitweibull.cdf-mean(fitweibull.cdf))/sd(fitweibull.cdf)
    spiWeibull[which(filtPrec == 0)]=(0-mean(fitweibull.cdf))/sd(fitweibull.cdf)


    ## Ajustar con la distribucion NORMAL---------------------------------------------
    fitnormal <- fitdistr(matrix(nonZeros),"normal")
    fitnormal.cdf<- (1-ProbZeros)*pnorm(nonZeros, mean = fitnormal$estimate[1], sd = fitnormal$estimate[2], lower.tail = TRUE, log.p = FALSE)

    ## Calculate SPI index if NORMAL
    spiNormal<-numeric(length(filtPrec))
    spiNormal[which(!filtPrec == 0)]=(fitnormal.cdf-mean(fitnormal.cdf))/sd(fitnormal.cdf)
    spiNormal[which(filtPrec == 0)]=(0-mean(fitnormal.cdf))/sd(fitnormal.cdf) 



    ## Ajustar con la distribucion NORMAL---------------------------------------------
    fitlognormal <- fitdistr(matrix(nonZeros),"log-normal")
    fitlognormal.cdf<- (1-ProbZeros)*plnorm(nonZeros, meanlog = fitlognormal$estimate[1], sdlog = fitlognormal$estimate[2], lower.tail = TRUE, log.p = FALSE)

    ## Calculate SPI index if NORMAL
    spilogNormal<-numeric(length(filtPrec))
    spilogNormal[which(!filtPrec == 0)]=(fitlognormal.cdf-mean(fitlognormal.cdf))/sd(fitlognormal.cdf)
    spilogNormal[which(filtPrec == 0)]=(0-mean(fitlognormal.cdf))/sd(fitlognormal.cdf) 

    ##################################################################################################################
    ############################################### Which one is better?   ###########################################
    ##################################################################################################################
    ## METHOD 1: COMPARISON OF EMPIRICAL CDF AGAINST FITTED CDF
    Fn<-ecdf(filtPrec) # Empirical CDF en filprec argmax(filtprec)
    sort(nonZeros)

    rmseGamma    <- rmse(Fn(nonZeros),fitgamma.cdf2)
    rmseWeibull  <- rmse(Fn(nonZeros),fitweibull.cdf)
    rmseNormal   <- rmse(Fn(nonZeros),fitnormal.cdf)
    rmselgNormal <- rmse(Fn(nonZeros),fitlognormal.cdf)

    #Save RMSE in a vector
    RMSE<-c(rmseGamma,rmseWeibull, rmseNormal,  rmselgNormal )

    #METHOD 2 : AKAIKE INFORMATION CRITERIUM
    #log likelihood for Gamma(alpha,scale=beta), X Data input
    LLabGamma <- function(pars,X){
      alpha <- pars[1] #parametro 1
      beta  <- pars[2]  #parametro 2
      return(-2*sum(log10(dgamma(X,alpha,rate=beta))) + 2*2 )
    }

    LLabWeibull <- function(pars,X){
      alpha <- pars[1]
      beta <- pars[2]
      return(-2*sum(log10(dweibull(X,alpha,scale=beta))) + 2*2)
    }

    LLabNormal <- function(pars,X){
       media <- fitnormal$estimate[1]
       dstd  <- fitnormal$estimate[2]
       return(-2*sum(log10(dnorm(X,mean = media, sd = dstd))) + 2*2)
    }


    LLablogNormal <- function(pars,X){
      media <- fitlognormal$estimate[1]
      dstd  <- fitlognormal$estimate[2]
      return(-2*sum(log10(dlnorm(X,meanlog = media, sdlog = dstd))) + 2*2)
    }

    #Call FUNCTIONS---------------------------------------
    akaiGamma   <- LLabGamma(fitgamma$estimate,    nonZeros)                 
    akaiWeibull <- LLabWeibull(fitweibull$estimate,nonZeros) 
    akaiNormal  <- LLabNormal(fitnormal$estimate,  nonZeros) 
    akaiLognormal <- LLablogNormal(fitlognormal$estimate, nonZeros)
    #Akaike Criterium informatop call in a vector
    AIC <-c(akaiGamma ,akaiWeibull,akaiNormal, akaiLognormal)  
    
    ###K-S TEST #####
    Dgamma    <- ks.test(x = Fn(nonZeros), y = fitgamma.cdf2)
    Dweibull  <- ks.test(x = Fn(nonZeros), y = fitweibull.cdf)
    Dnormal   <- ks.test(x = Fn(nonZeros), y = fitnormal.cdf)
    Dlgnormal <- ks.test(x = Fn(nonZeros), y = fitlognormal.cdf)
    #almacenamos en un array
    KSTEST <- c(Dgamma$statistic[1], Dweibull$statistic[1], Dnormal$statistic[1], Dlgnormal$statistic[1] )

    #MINIMUN VALUES
    minRMSE<- which.min(RMSE)
    minKS <-  which.min(KSTEST)
    minAIC <- which.min(AIC)

    #select the best
    indexSelected <- which.min(abs((RMSE-RMSE[minRMSE]) /RMSE[minRMSE])+ abs((KSTEST-KSTEST[minKS])/KSTEST[minKS]) + 
                                  abs((AIC-AIC[minAIC]) /AIC[minAIC])  )

    if (indexSelected==1){
     idajuste <- 'Gamma'
     spi<-spiGamma
     BEST_FIT = c('YES','NO','NO','NO')
     print('ES GAMMA!!!')

    } else if (indexSelected==2) {
     spi<-spiWeibull
     idajuste <- 'Weibull'
     BEST_FIT = c('NO','YES','NO','NO')
     print('CDF WEIBULL')
    } else if (indexSelected==3) {
     spi<-spiNormal
     idajuste<- 'Normal'
     BEST_FIT = c('NO','NO','YES','NO')
     print('ES NORMAL!!!')

    } else if (indexSelected==4) {
     spi<-spilogNormal
     idajuste<- 'logNormal'
     BEST_FIT = c('NO','NO','NO','YES')
     print('ES logNORMAL!!!')
    }
    
    
    distrs = c('Gamma','Weibull','Normal', 'lgNormal')
    ##ALmacenado valores
    c_names  <- append(c_names, rep(name_est, times = 4))
    c_distr  <- append(c_distr, distrs)
    c_rmse   <- append(c_rmse, RMSE)  
    c_aic    <- append(c_aic, AIC)   
    c_kstest <- append(c_kstest, KSTEST)
    c_best   <- append(c_best, BEST_FIT)    

    #PLOT comulative ECD Y CDF comparation ###########################################################################################
    pdf(paste(ruta_salida, name_est,'_CDFs.pdf', sep = ''), width=6.3, height=5.0)
    plot(sort(nonZeros),sort(Fn(nonZeros)),type = 's' , ylab = 'Probabilidad Acumulada', xlab = 'Precipitación[mm]', cex =0.5, col = 'black', lwd=1.3)
    lines(sort(nonZeros), sort(fitgamma.cdf2), col = 'blue', lwd = 1.8, lty = 3)
    lines(sort(nonZeros), sort(fitweibull.cdf), col = 'red', type  = 'l', lty = 2, lwd = 1.8)
    lines(sort(nonZeros), sort(fitnormal.cdf), col = 'purple', type  = 'l', lty = 4, lwd = 1.8)
    lines(sort(nonZeros), sort(fitlognormal.cdf), col = 'orange', type  = 'l', lty = 6, lwd = 1.8)
    #legend('bottomright', legend = '   CDF Empirico', col = 'black', pch = 4,  ,box.lty=0, cex = 0.8)
    legend('bottomright', legend=c('CDF Empirico','CDF Gamma','CDF Weibull','CDF Normal', 'CDF Lognormal'), col= c('black','blue','red','purple', 'orange'), lty = c(1,3,2,4,6), 
          lwd = c(1.3,2.0,2.0,2.0, 2.0) ,box.lty=0,  cex = 0.8)
    legend('topleft', legend= paste(name_est,': Mejor Ajuste ',idajuste,sep =''), col= 'black', cex = 0.75,box.lty=0)
    dev.off()

    #PLOT BEST SPI SELECTED #######################################################################################
    #Creamos data time series
    fechas <-seq(as.Date('1980/1/1'), as.Date('2014/12/31'), 'month')
    #complete array
    if (nfil  == 3) {
      spiComplete = spi
     }else{
      spiComplete <- c(init, spi,fin)  
     } 

    #SAVE SPICOMPLETE
    df_spi <- data.frame(fechas, spiComplete)
    write.table(df_spi, file = paste(ruta_salida,name_est,'_spi.scv', sep = ''), row.names = FALSE, col.names = FALSE, sep = '\t')
    rm(df_spi)

    #PLot SPI COMPLETE
    pdf(paste(ruta_salida,name_est,'_SPI.pdf', sep = ''), width=9.0, height=4.0)
    plot(fechas,spiComplete,type='h',col='black', lwd = 1.5, xlab = '', ylab = 'SPI', xaxt = 'n')
    #lines(fechas,c(init,spiGamma2,fin),type='h',col='red', lwd = 0.7, xlab = '', ylab = 'SPI', xaxt = 'n')
    axis.Date(1, at  = fechas, cex.axis=1, format = '%Y')
    legend('bottomright', legend='3 MESES', col='black', box.lty=0,cex = 0.8)
    par
    abline(h = 0, untf = FALSE)
    abline(h = 1, untf = FALSE, lwd=0.4, lty=2)
    abline(h = -1, untf = FALSE, lwd=0.4, lty=2)
    dev.off()
}

#CREATE DATA FRAME FITS and save #################################################################################################
df <- data.frame(c_names, c_distr, c_rmse, c_aic, c_kstest, c_best)
write.table(df,file  = paste(ruta_salida,'reg_info.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)
print('###############################################')
print('#              DONE RIGHT ThinkPad            #')
print('###############################################')



