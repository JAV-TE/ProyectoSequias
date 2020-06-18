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
library(fitdistrplus)
require(RTOMO) #si

#Seleccion  de la FDP para el calculo de SPI
library(eva) # analisis de valores extremos
library(MASS) # ajuste de distribuciones estadisticas
library(fitdistrplus) # MORES STATISTICAL DISTRIBUTIONS OPTIONS
library(hydroGOF) # for RMSE
library(extRemes) # ALTERNATIVE TO EVA

#sedebe especificar
library(nsRFA)
library(lmomco)
library(scatterplot3d)
#plot(lattice objects)
library(latticeExtra) # sumar graficas
library(ggpubr) #varias graficas en uno

#ggplot and ggpmisc
library(ggplot2)
library(ggpmisc)


#extrem values
library(ismev)
library(FAdist)


################################################open file#################################################

#ruta_carpeta  <- paste('/home/Physics/TESIS/Agrupacion/Estaciones_Agrupadas/',region,'/',sep='')
ruta_salida0 <- '/home/Physics/pro-seq/1_getSPIindex/3_CAL-SPI1/OUTCOPULA/'
ruta_carpeta <- '/home/Physics/pro-seq/1_getSPIindex/3_CAL-SPI1/1_calculo_SPI/SPI/est-SPI-65-percent/'

print('%%%%%%%--Inserte la Region correpondintes---- %%%%%%')
region <- readline()
ruta_carpeta <- paste(ruta_carpeta,region,'/',sep='')

ruta_salida  <- paste(ruta_salida0,region,sep='')

dir.create(ruta_salida)

ruta_salida <- paste(ruta_salida,'/',sep='')


list_est <- dir(ruta_carpeta, pattern = '_spi.scv')

#VeCTORES VACIOS PARA ALMACENAR
c_ests    <- NULL
c_corrk   <- NULL
c_corrp   <- NULL
c_corrsp  <- NULL

c_marg_D  <- NULL
c_marg_S  <- NULL

c_copula  <- NULL

c_pearson <- NULL
c_fitline <- NULL
c_rmsecop <- NULL

c_msevRad   <- NULL
c_msevGrade <- NULL
 

for (estation in list_est) {
   name_est <- gsub('_spi.scv','', estation)
   df_spi <- read.table(paste(ruta_carpeta, estation ,sep = '') ,sep ='\t')
   spi <- df_spi$V2
   times <- c(1:length(spi)) 

   #VeCTORES VACIOS PARA ALMACENAR
   c_names   <- NULL
   c_distr   <- NULL
   c_Drmse   <- NULL
   c_Daic    <- NULL
   c_Dkstest <- NULL
   c_Dbest   <- NULL
   c_Srmse   <- NULL
   c_Saic    <- NULL
   c_Skstest <- NULL
   c_Sbest   <- NULL

   ###############3#####################################################
   ########################### Selection of Drughts Events #############
   #####################################################################
   idx<-which(spi < 0)
   dif = idx[-1]-idx[1:length(idx)-1] #dalculate diferences betwen steps (duration months)
   timeBetweenEvents=mean(dif[which(dif>1)]) ## Average Time Between Drought Events

   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULO DE DURACION, SEVERIDAD Y RECURRENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   index2 <- which(dif>1)
   index3 <- append(1,index2)

   intera_time <- NULL #tiempo entre eventos en meses

   for (i in seq(length(index3)-1)){
      if (i ==1){
        sum_init <- sum(dif[index3[i]:index3[i+1]])
        intera_time <- append(intera_time, sum_init)
        print(sum_init)
      }else {
        sum_rot <- sum(dif[(index3[i]+1):index3[i+1]]) 
        intera_time <- append(intera_time, sum_rot)
      }
     }

   #MEAN  INRETARIIVAL TIME %
   mean_interar <- mean(intera_time)

   rm(dif)  

   #Compute Duration & Severity
   itS = 1
   itD = 0
   #white arrays
   severity<- 0*numeric(length(idx)) 
   duration<- 0*numeric(length(idx)) 
   
   #idx indices de durability caculate Duration Severity events
   for (i in seq(1,length(idx)-1)){
     itD = itD + 1 
     if(idx[i+1]-idx[i]==1 & i < (length(idx)-1)){
       severity[itS] = severity[itS] - spi[idx[i]]
     }else if (idx[i+1]-idx[i]==1 & i==(length(idx)-1)){
       duration[itS] = itD + 1
       severity[itS] = severity[itS] - spi[idx[i]]-spi[idx[i+1]]
     }else{ 
       duration[itS] = itD
       severity[itS] = severity[itS] - spi[idx[i]]
       itS = itS + 1 
     itD = 0
    }}

   #seleccionamos solo  mayores a cero
   idx1<-which(duration > 0)
   idx2<-which(severity > 0)
  
   duration = duration[idx1] #######33
   severity = severity[idx2] #######


   #CALCULO DE LINEA DE TENDENCIA PAra SEVERIDAD
   #anulas severidades 1
   indexsev <- which(duration>3)
   severity0 = severity[indexsev]

   nfil = 6 #STEP
   func_smooth <- function(severity, nfil){
     filtSev = filter(severity, sides=2, filter=rep(1/nfil,nfil)) # moving average sides = 2 centered
     #filtSev = filtSev[3:417]
     #init = c(NA,NA)
     #fin =  c(NA,NA,NA)
     return(filtSev)
    }

   filtSev <- func_smooth(severity0, nfil)
   times <- seq(1, length(filtSev))
   fit     <- lm(filtSev   ~   times )
   m = round(fit$coefficients[[2]], 4)
   m_str = as.character(m)
   m_grade = atan(m)*180/pi
   #plot trend lines
   pdf(paste(ruta_salida,name_est ,'_TRENDSev.pdf', sep = ''), width= 7.5, height= 4.5)
   plot(times, severity0,ylim =c((min(severity)-0.5), (max(severity)+0.5)), xlab = 'Tiempo', ylab= 'Severidad', lwd = 1.3, type='o')
   lines(times, filtSev, lwd = 1.3, col = 'red', lty = 5)
   abline(fit,  lwd =1.3 ,col ='black', lty = 1) 
   legend('topleft', legend=c('Tendencia','Media Movil 6'), col= c('black','red'), lty = c(1,5),lwd = c(1.3,1.2) ,box.lty=0)
   legend('bottomright', legend = paste('Estación ', name_est,': m = ',m_str, sep = '' ), box.lty=0)
   dev.off()


   ############Imprimir datos observados de DURAcion SEVERIDAD T RECURRENCIA

   rm(idx1)
   rm(idx2)

   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIN DE CALCULO DE DURACION, SEVERIDAD Y RECURRENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ##CORRELACION ++++++++++++++++++++++++++++++++++++++++++++++++++++
   cor_pearson  <- round(cor(duration, severity, method='pearson'), 3)
   cor_kendall  <- round(cor(duration, severity, method='kendall'), 3)
   cor_spear <- round(cor(duration, severity, method='spearman'), 3)
   cor_per  <- as.character(cor_pearson,2)
   cor_ken  <- as.character(cor_kendall,2)
   cor_sper <- as.character(cor_spear,2)

   #plots de DUR vs Severidad	
   pdf(paste(ruta_salida,name_est ,'_dur_sev.pdf', sep = ''), width=6.3, height=5.0)
   plot(duration,severity,xlab="Duración [Meses]", ylab="Severidad",pch = 16, type ='p', xlim=c(1, 10), ylim=c(1, 10))
   axis(1, at=seq(1,10, by=1), labels=seq(1,10,1), lwd=0, lwd.ticks=2)
   axis(2, at=seq(1,10, by=1), labels=seq(1,10,1), lwd=0, lwd.ticks=2)
   legend('topleft', legend= paste('Estación :',name_est,sep =' '), col= 'black', cex = 0.80,box.lty=0)
   legend('bottomright', legend=paste('per = ',cor_per, ', ken =', cor_ken), col= 'black', cex = 0.75,box.lty=0)
   dev.off()

   #stop('HASTA AQUI TODO BIEN')
   #rm(symbol1)
   #rm(cor_char)

   #########################################################################################
   ###########################  Fit of severity and duration  ##############################
   #########################################################################################
   duration_art <- seq(0.0001,14.0,length.out=100)
   severity_art <- seq(0.0001,14.0,length.out=100)


   ##DEFINE gumbel distr 
   dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
   pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
   qgumbel <- function(p, a, b) a-b*log(-log(p))
   ########################################AJUSTES PARA EVENTOS DE DURACION  ########################################################################################
   ## Fitting duration to a Gamma distribution
   fitgamma.duration <- fitdistr(duration,densfun = "gamma")
   fitgamma_durationCDF<-pgamma(duration, shape = fitgamma.duration$estimate[1], rate = fitgamma.duration$estimate[2], lower.tail = TRUE, log.p = FALSE)
   fitgamma_durCDFart <-pgamma(duration_art, shape = fitgamma.duration$estimate[1], rate = fitgamma.duration$estimate[2], lower.tail = TRUE, log.p = FALSE)  

   ## Fitting duration to a Weibull distribution
   fitweibull.duration <- fitdistr(duration,"weibull")
   fitweibull_durationCDF<-pweibull(duration, shape = fitweibull.duration$estimate[1], scale = fitweibull.duration$estimate[2],   lower.tail = TRUE, log.p = FALSE)
   fitweibull_durCDFart<-pweibull(duration_art, shape = fitweibull.duration$estimate[1], scale = fitweibull.duration$estimate[2],   lower.tail = TRUE, log.p = FALSE)

   ## Fitting duration to a Exponential distribution
   fitexp.duration <- fitdistr(duration,"exponential")
   fitexp_durationCDF<-pexp(duration, rate = fitexp.duration$estimate[1], lower.tail = TRUE, log.p = FALSE)
   fitexp_durCDFart  <-pexp(duration_art, rate = fitexp.duration$estimate[1], lower.tail = TRUE, log.p = FALSE)

   ## Fitting duration to a Lognormal distribution
   fitlogn.duration <- fitdistr(duration,"lognormal")
   fitlogn_durationCDF<- plnorm(duration, meanlog  = fitlogn.duration$estimate[1], sdlog  = fitlogn.duration$estimate[2],lower.tail = TRUE, log.p = FALSE)
   fitlogn_durCDFart  <- plnorm(duration_art, meanlog  = fitlogn.duration$estimate[1], sdlog  = fitlogn.duration$estimate[2],lower.tail = TRUE, log.p = FALSE)

   ## Fitting duration to a gumbel distribution
   fitgum.duration <- fitdist(duration, "gumbel", start=list(a=10, b=10), method = 'mle')
   fitgum_durationCDF <- pgumbel(duration,fitgum.duration$estimate[1], fitgum.duration$estimate[2] )
   fitgum_durCDFart <- pgumbel(duration_art,fitgum.duration$estimate[1], fitgum.duration$estimate[2] )  

   ## Fitting duration to a GEV distribution
   lmrgevD <- lmoms(duration)
   fitgev.duration = pargev(lmrgevD) #mle fit
   fitgev_durationCDF <- cdfgev(duration,fitgev.duration)
   fitgev_durCDFart <- cdfgev(duration_art,fitgev.duration)

   ## Fitting duration to a Geometric distribution
   fitgeom.duration <- fitdistr(duration,densfun = "geometric")
   fitgeom_durationCDF<-pgeom(duration, prob = fitgeom.duration$estimate[1], lower.tail = TRUE, log.p = FALSE)
   fitgeom_durCDFart<-pgeom(duration_art, prob = fitgeom.duration$estimate[1], lower.tail = TRUE, log.p = FALSE)

   ## Fitting duration to a Normal distribution
   fitnorm.duration <- fitdistr(duration,"normal")
   fitnorm_durationCDF<- pnorm(duration, mean = fitnorm.duration$estimate[1], sd = fitnorm.duration$estimate[2], lower.tail = TRUE, log.p = FALSE)
   fitnorm_durCDFart  <- pnorm(duration_art, mean = fitnorm.duration$estimate[1], sd = fitnorm.duration$estimate[2], lower.tail = TRUE, log.p = FALSE)


   #####################################AJUSTAR A EVENTOS DE SEVERIDAD ################################################################################################
   ## Fitting severity to a Gamma distribution
   fitgamma.severity     <- fitdistr(severity,"gamma")
   fitgamma_severityCDF  <-pgamma(severity, shape = fitgamma.severity$estimate[1], rate = fitgamma.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)
   fitgamma_sevCDFart  <-pgamma(severity_art, shape = fitgamma.severity$estimate[1], rate = fitgamma.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)
  
   ## Fitting severity to a Weibull distribution
   fitweibull.severity   <- fitdistr(severity,"weibull")
   fitweibull_severityCDF<-pweibull(severity, shape = fitweibull.severity$estimate[1], scale = fitweibull.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)
   fitweibull_sevCDFart  <-pweibull(severity_art, shape = fitweibull.severity$estimate[1], scale = fitweibull.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)
    
   ## Fitting severity to a Exponential distribution
   fitexp.severity      <- fitdistr(severity,"exponential")
   fitexp_severityCDF   <-pexp(severity, rate = fitexp.severity$estimate[1], lower.tail = TRUE, log.p = FALSE)
   fitexp_sevCDFart     <-pexp(severity_art, rate = fitexp.severity$estimate[1], lower.tail = TRUE, log.p = FALSE)

   ## Fitting severity to a Lognormal distribution
   fitlogn.severity    <- fitdistr(severity ,"lognormal")
   fitlogn_severityCDF <-plnorm(severity, meanlog  = fitlogn.severity$estimate[1], sdlog  = fitlogn.severity$estimate[2],lower.tail = TRUE, log.p = FALSE)
   fitlogn_sevCDFart   <-plnorm(severity_art, meanlog  = fitlogn.severity$estimate[1], sdlog  = fitlogn.severity$estimate[2],lower.tail = TRUE, log.p = FALSE)


   fitgum.severity <- fitdist(severity, "gumbel", start=list(a=10, b=10), method = 'mle')
   fitgum_severityCDF <- pgumbel(severity,fitgum.severity$estimate[1], fitgum.severity$estimate[2] )
   fitgum_sevCDFart <- pgumbel(severity_art,fitgum.severity$estimate[1],fitgum.severity$estimate[2] ) 

   lmrgevS <- lmoms(severity)
   fitgev.severity     = pargev(lmrgevS) #mle fit
   fitgev_severityCDF <- cdfgev(severity,fitgev.severity)
   fitgev_sevCDFart   <- cdfgev(severity_art,fitgev.severity)

   ## Fitting severity to a Geometric distribution
   fitgeom.severity <- fitdistr(severity,"geometric")
   fitgeom_severityCDF<-pgeom(severity, prob = fitgeom.severity$estimate[1], lower.tail = TRUE, log.p = FALSE)
   fitgeom_sevCDFart<-pgeom(severity_art, prob = fitgeom.severity$estimate[1], lower.tail = TRUE, log.p = FALSE)

   #Fitting duration to a normal distribution
   fitnorm.severity  <- fitdistr(severity,"normal")
   fitnorm_severityCDF <- pnorm(severity, mean = fitnorm.severity$estimate[1], sd = fitnorm.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)
   fitnorm_sevCDFart   <- pnorm(severity_art, mean = fitnorm.severity$estimate[1], sd = fitnorm.severity$estimate[2], lower.tail = TRUE, log.p = FALSE)

   ############################################################### DETERMINAR CUAL ES EL MEJOR AJUSTE ##################################################################
   ## METHOD 1: COMPARISON OF EMPIRICAL CDF AGAINST FITTED CDF od duration and severity
   Fn.duration<-ecdf(duration)
   RMSE_du =c(numeric(8))
   RMSE_du[1]   <- rmse(Fn.duration(duration),fitgamma_durationCDF)
   RMSE_du[2]   <- rmse(Fn.duration(duration),fitweibull_durationCDF)
   RMSE_du[3]   <- rmse(Fn.duration(duration),fitexp_durationCDF)
   RMSE_du[4]   <- rmse(Fn.duration(duration),fitlogn_durationCDF)
   RMSE_du[5]   <- rmse(Fn.duration(duration),fitgum_durationCDF)
   RMSE_du[6]   <- rmse(Fn.duration(duration) ,fitgev_durationCDF)
   RMSE_du[7]   <- rmse(Fn.duration(duration),fitgeom_durationCDF)
   RMSE_du[8]   <- rmse(Fn.duration(duration) ,fitnorm_durationCDF)


   Fn.severity<-ecdf(severity)
   RMSE_se =c(numeric(8))
   RMSE_se[1]  <- rmse(Fn.severity(severity),fitgamma_severityCDF)
   RMSE_se[2] <- rmse(Fn.severity(severity),fitweibull_severityCDF)
   RMSE_se[3]    <- rmse(Fn.severity(severity),fitexp_severityCDF)
   RMSE_se[4]   <- rmse(Fn.severity(severity),fitlogn_severityCDF)
   RMSE_se[5]   <- rmse(Fn.severity(severity),fitgum_severityCDF)
   RMSE_se[6]   <- rmse(Fn.severity(severity) ,fitgev_severityCDF )
   RMSE_se[7]   <- rmse(Fn.severity(severity),fitgeom_severityCDF)
   RMSE_se[8]   <- rmse(Fn.severity(severity) ,fitnorm_severityCDF )

   ## METHOD 2 : AKAIKE INFORMATION CRITERIUM#########################################

   #log likelihood for Gamma(alpha,scale=beta), X Data input
   LLabGamma <- function(pars,X){
      alph <- pars[1] #parametro 1
      beta  <- pars[2]  #parametro 2
      return(-2*sum(log(dgamma(X,shape = alph,rate=beta))) + 2*2) 
      } 

   LLabWeibull <- function(pars,X){
      shap <- pars[1]
      scal <- pars[2]
      return(-2*sum(log(dweibull(X,shape = shap,scale=scal))) + 2*2)
      }

   LLabExp <- function(pars,X){
      rat <- pars[1]
      return(-2*sum(log(dexp(X,rate=rat)))  + 2*1)
      }

   LLabLogn <- function(pars,X){
      a <- pars[1]
      b <- pars[2]
      return(-2*sum(log(dlnorm(X,meanlog = a, sdlog = b))) + 2*2)
      }

   LLabGumbel <- function(pars,X){
      par1 <- pars[1]
      par2  <- pars[2]
      return(-2*sum(log(dgumbel(X,par1, par2 ))) + 2*2)
      }

   LLabGev <- function(fitgev.duration,X){
      return(-2*sum(log(pdfgev(X,fitgev.duration))) + 3*2)
      }

   LLabGeom <- function(pars,X){
    # log likelihood for Gamma(alpha,scale=beta)
    # X is data
    probb <- pars[1]
    return(-2*sum(log(dgeom(X,prob=probb)))  + 2*1)
    }


   LLabNormal <- function(pars,X){
      media <- pars[1]
      dstd  <- pars[2]
      return(-2*sum(log(dnorm(X,mean = media, sd = dstd))) + 2*2)
      }




   AIC_du <- c(numeric(6))
   AIC_du[1]   <- LLabGamma(fitgamma.duration$estimate,duration)
   AIC_du[2] <- LLabWeibull(fitweibull.duration$estimate,duration)
   AIC_du[3]     <- LLabExp(fitexp.duration$estimate,duration)
   AIC_du[4]   <- LLabLogn(fitlogn.duration$estimate,duration)
   AIC_du[5]    <- LLabGumbel(fitgum.duration$estimate,duration)
   AIC_du[6]   <- LLabGev(fitgev.duration,duration)
   AIC_du[7]    <- LLabGeom(fitgeom.duration$estimate,duration)
   AIC_du[8]    <- LLabNormal(fitnorm.duration$estimate,duration)


   AIC_se = c(numeric(6))
   AIC_se[1] <- LLabGamma(fitgamma.severity$estimate,severity)
   AIC_se[2] <- LLabWeibull(fitweibull.severity$estimate,severity)
   AIC_se[3] <- LLabExp(fitexp.severity$estimate,severity)
   AIC_se[4] <- LLabLogn(fitlogn.severity$estimate,severity)  
   AIC_se[5] <- LLabGumbel(fitgum.severity$estimate,severity)
   AIC_se[6] <- LLabGev(fitgev.severity ,severity)
   AIC_se[7] <- LLabGeom(fitgeom.severity$estimate,severity)  
   AIC_se[8] <- LLabNormal(fitnorm.severity$estimate,severity)  

   ###METHOD-3  K-S TEST #####
   #Alternativo <- c('greater')
   #DDgamma    <- ks.test(x = Fn.duration(duration), y = fitgamma_durationCDF, alternative=Alternativo, exact = FALSE)
   #DDweibull  <- ks.test(x = Fn.duration(duration), y = fitweibull_durationCDF,alternative=Alternativo, exact = FALSE)
   #DDexp      <- ks.test(x = Fn.duration(duration), y = fitexp_durationCDF, alternative=Alternativo, exact = FALSE)
   #DDlgnormal <- ks.test(x = Fn.duration(duration), y = fitlogn_durationCDF, alternative=Alternativo, exact = FALSE)
   #DDnormal   <- ks.test(x = Fn.duration(duration), y = fitnorm_durationCDF, alternative=Alternativo, exact = FALSE)

   #DSgamma    <- ks.test(x = Fn.severity(severity), y = fitgamma_severityCDF, alternative=Alternativo, exact = FALSE)
   #DSweibull  <- ks.test(x = Fn.severity(severity), y = fitweibull_severityCDF, alternative=Alternativo, exact = FALSE)
   #DSexp      <- ks.test(x = Fn.severity(severity), y = fitexp_severityCDF, alternative=Alternativo, exact = FALSE)
   #DSlgnormal <- ks.test(x = Fn.severity(severity), y = fitlogn_severityCDF, alternative=Alternativo, exact = FALSE)
   #DSnormal   <- ks.test(x = Fn.severity(severity), y = fitnorm_severityCDF, alternative=Alternativo, exact = FALSE)

   #almacenamos en un array TOOOOODOs los resultados
   #JOIN RMSE
   #RMSE_du<-c(rmseGamma.duration,rmseWeibull.duration, rmseExp.duration , rmseLogn.duration, rmseNorm.duration) 
   #RMSE_se<-c(rmseGamma.severity ,rmseWeibull.severity, rmseExp.severity, rmseLogn.severity, rmseNorm.severity)

   #JOIN AKAIKE
   #AIC_du <- c(LLabGamma.duration, LLabWeibull.duration,LLabExp.duration,LLabLogn.duration,LLabNormal.duration) 
   #AIC_se <- c(LLabGamma.severity, LLabWeibull.severity,LLabExp.severity,LLabLogn.severity,LLabNormal.severity)
     
   getInf <- is.infinite(AIC_du)
   AIC_du[getInf]=NaN
   getInf <- is.infinite(AIC_se)
   AIC_se[getInf]=NaN

   #JOIN KS_TEST
   #KStest_du <- c(DDgamma$statistic[1], DDweibull$statistic[1], DDexp$statistic[1], DDlgnormal$statistic[1], DDnormal$statistic[1])
   #KStest_se <- c(DSgamma$statistic[1], DSweibull$statistic[1], DSexp$statistic[1], DSlgnormal$statistic[1], DSnormal$statistic[1])
    
   #SELECT THE BEST ----
   #MINIMUN VALUES duration
   DminRMSE0<- which.min(RMSE_du)
   #DminKS <-  which.min(KStest_du)
   DminAIC <- which.min(AIC_du)
   DminRMSE <-  which.min(AIC_du) 
   #MINIMUN VALUES severity
   #SminRMSE<- which.min(RMSE_se)
   #SminKS <-  which.min(KStest_se
   SminRMSE0 <-  which.min(RMSE_se)
   SminRMSE <-  which.min(AIC_se)

   SminAIC <-  which.min(AIC_se)
   #select the best duration
   DindexSelected <- which.min(abs((RMSE_du-RMSE_du[DminRMSE]) /RMSE_du[DminRMSE])+ abs((AIC_du-AIC_du[DminAIC]) /AIC_du[DminAIC])  )
   #select the best severity
   SindexSelected <- which.min(abs((RMSE_se-RMSE_se[SminRMSE]) /RMSE_se[SminRMSE])+ abs((AIC_se-AIC_se[SminAIC]) /AIC_se[SminAIC])  )

   #los indices seleccionado estan mal pero ahora solo tomaremos encuenta el RMSE por la falta de datos

   #PARA DURATION------------------------------------------------------###################################
   print('Duracion-------------------------')
   if (DminRMSE==1){
      idajusteD <- 'Gamma'
      fit.distDUR <- fitgamma.duration
      durationOrdered    <-duration[order(fitgamma_durationCDF)] 
      durationOrderedCDF <-fitgamma_durationCDF[order(fitgamma_durationCDF)]
      durationCDF        <-fitgamma_durationCDF
      dur_distr <- dgamma(duration_art, shape = fitgamma.duration$estimate[1], rate = fitgamma.duration$estimate[2], log = FALSE)
      BEST_FIT = c('YES','NO','NO','NO','NO','NO','NO','NO')
      legendaD <- c('Observado','Gamma (*)','Weibull', 'Exponencial', 'logNormal', 'Gumbel','GEV','Geometrica','Normal')
      print('ES GAMMA!!!')

   } else if (DminRMSE==2) {
      fit.distDUR <- fitweibull.duration
      durationOrdered    <-duration[order(fitweibull_durationCDF)] 
      durationOrderedCDF <-fitweibull_durationCDF[order(fitweibull_durationCDF)]
      durationCDF        <-fitweibull_durationCDF 
      dur_distr <- dweibull(duration_art, shape = fitweibull.duration$estimate[1], scale = fitweibull.duration$estimate[2],    log = FALSE)
      idajusteD <- 'Weibull'
      BEST_FIT = c('NO','YES','NO','NO','NO','NO','NO','NO')
      legendaD <- c('Observado','Gamma','Weibull (*)', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal')
      print('CDF WEIBULL')

   } else if (DminRMSE==3) {
      fit.distDUR <- fitexp.duration
      durationOrdered    <-duration[order(fitexp_durationCDF)] 
      durationOrderedCDF <-fitexp_durationCDF[order(fitexp_durationCDF)]
      durationCDF        <-fitexp_durationCDF
      dur_distr  <-dexp(duration_art, rate = fitexp.duration$estimate[1], log = FALSE)
      idajusteD<- 'Exponencial'
      BEST_FIT = c('NO','NO','YES','NO','NO','NO','NO','NO')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial (*)', 'logNormal','Gumbel','GEV','Geometrica','Normal')
      print('ES exponential!!!')

   } else if (DminRMSE==4) {
      fit.distDUR <- fitlogn.duration
      durationOrdered    <-duration[order(fitlogn_durationCDF)] 
      durationOrderedCDF <-fitlogn_durationCDF[order(fitlogn_durationCDF)]
      durationCDF        <-fitlogn_durationCDF 
      dur_distr <- dlnorm(duration_art, meanlog  = fitlogn.duration$estimate[1], sdlog  = fitlogn.duration$estimate[2], log = FALSE)
      idajusteD<- 'logNormal'
      BEST_FIT = c('NO','NO','NO','YES','NO','NO','NO','NO')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal (*)', 'Gumbel','GEV','Geometrica','Normal')
      print('ES logNORMAL!!!')

   } else if (DminRMSE==5) {
      fit.distDUR <- fitgum.duration
      #durationOrdered    <-duration[order(fitnorm_durationCDF)] 
      #durationOrderedCDF <-fitnorm_durationCDF[order(fitnorm_durationCDF)]
      durationCDF        <-fitgum_durationCDF
      dur_distr <- dgumbel(duration_art,fitgum.duration$estimate[1], fitgum.duration$estimate[2] )
      idajusteD<- 'Gumbel'
      BEST_FIT = c('NO','NO','NO','NO', 'YES','NO','NO','NO')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal', 'Gumbel (*)','GEV','Geometrica','Normal')
      print('ES GUMBEL!!!')

   } else if (DminRMSE==6) {
      fit.distDUR <- fitgev.duration
      #durationOrdered    <-duration[order(fitnorm_durationCDF)] 
      #durationOrderedCDF <-fitnorm_durationCDF[order(fitnorm_durationCDF)]
      durationCDF        <-fitgev_durationCDF 
      dur_distr <- pdfgev(duration_art,fitgev.duration)
      idajusteD<- 'GEV'
      BEST_FIT = c('NO','NO','NO','NO','NO','YES','NO','NO')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal', 'Gumbel','GEV (*)','Geometrica','Normal')
      print('ES GEV!!!')


   } else if (DminRMSE==7) {
      fit.distDUR <- fitgeom.duration
      #durationOrdered    <-duration[order(fitexp_durationCDF)] 
      #durationOrderedCDF <-fitexp_durationCDF[order(fitexp_durationCDF)]
      durationCDF        <-fitgeom_durationCDF
      dur_distr  <-dgeom(duration_art, prob = fitgeom.duration$estimate[1], log = FALSE)
      idajusteD<- 'Geometrica'
      BEST_FIT = c('NO','NO','NO','NO','NO','NO','YES','NO')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica (*)','Normal')
      print('ES Geometrico!!!')
   } else if (DminRMSE==8) {
      fit.distDUR <- fitnorm.duration
      #durationOrdered    <-duration[order(fitexp_durationCDF)] 
      #durationOrderedCDF <-fitexp_durationCDF[order(fitexp_durationCDF)]
      durationCDF        <-fitnorm_durationCDF
      dur_distr  <-dnorm(duration_art, mean = fitnorm.duration$estimate[1], sd = fitnorm.duration$estimate[2], log = FALSE)
      idajusteD<- 'Normal'
      BEST_FIT = c('NO','NO','NO','NO','NO','NO','NO','YES')
      legendaD <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal (*)')
      print('ES Normal!!!')


   }






   # LISTA DE DISTRIBUCIONES
   distrs  <- c('Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal')

   #save data frame duration erros
   df <- data.frame(distrs, RMSE_du,AIC_du ,BEST_FIT)
   write.table(df,file  = paste(ruta_salida,name_est,'_DUR_inf.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)
 

   #######################################PARA SEVERITY----------------------------------------------############################################
   print('Severidad-------------------------')

   if (SminRMSE==1){
      idajusteS <- 'Gamma'
      fit.distSEV <- fitgamma.severity
      severityOrdered    <-severity[order(fitgamma_severityCDF)] 
      severityOrderedCDF <-fitgamma_severityCDF[order(fitgamma_severityCDF)]
      severityCDF        <-fitgamma_severityCDF
      sev_distr <- dgamma(severity_art, shape = fitgamma.severity$estimate[1], rate = fitgamma.severity$estimate[2], log = FALSE)
      BEST_FIT <- c('YES','NO','NO','NO','NO','NO','NO','NO')
      legendaS <- c('Observado','Gamma (*)','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal')
      print('ES GAMMA!!!')

   } else if (SminRMSE==2) {
      fit.distSEV <- fitweibull.severity
      severityOrdered    <-severity[order(fitweibull_severityCDF)] 
      severityOrderedCDF <-fitweibull_severityCDF[order(fitweibull_severityCDF)]
      severityCDF        <-fitweibull_severityCDF
      sev_distr <- dweibull(severity_art, shape = fitweibull.severity$estimate[1], scale = fitweibull.severity$estimate[2], log = FALSE)
      idajusteS <- 'Weibull'
      BEST_FIT = c('NO','YES','NO','NO','NO','NO','NO','NO')
      legendaS <- c('Observado','Gamma','Weibull (*)', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal')
      print('CDF WEIBULL')

   } else if (SminRMSE==3) {
      fit.distSEV <- fitexp.severity
      severityOrdered    <-duration[order(fitexp_severityCDF)] 
      severityOrderedCDF <-fitexp_severityCDF[order(fitexp_severityCDF)]
      severityCDF        <-fitexp_severityCDF 
      sev_distr <- dexp(severity_art, rate = fitexp.severity$estimate[1], log= FALSE)
      idajusteS<- 'Exponencial'
      BEST_FIT = c('NO','NO','YES','NO','NO','NO','NO','NO')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial (*)', 'logNormal','Gumbel','GEV','Geometrica','Normal')
      print('ES exponential!!!')

   } else if (SminRMSE==4) {
      fit.distSEV <- fitlogn.severity
      severityOrdered<-severity[order(fitlogn_severityCDF)] 
      severityOrderedCDF<-fitlogn_severityCDF[order(fitlogn_severityCDF)]
      severityCDF<-fitlogn_severityCDF
      sev_distr <- dlnorm(severity_art, meanlog  = fitlogn.severity$estimate[1], sdlog  = fitlogn.severity$estimate[2],log = FALSE)
      idajusteS<- 'log-Normal'
      BEST_FIT = c('NO','NO','NO','YES','NO','NO','NO','NO')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal (*)','Gumbel','GEV','Geometrica','Normal')
      print('ES logNORMAL!!!')


   } else if (SminRMSE==5) {
      fit.distSEV <- fitgum.duration
      #severityOrdered<-severity[order(fitnorm_severityCDF)] 
      #severityOrderedCDF<-fitnorm_severityCDF[order(fitnorm_severityCDF)]
      severityCDF<-fitgum_severityCDF
      sev_distr <- dgumbel(severity_art,fitgum.severity$estimate[1], fitgum.severity$estimate[2] )
      idajusteS<- 'Gumbel'
      BEST_FIT = c('NO','NO','NO','NO','YES','NO','NO','NO')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel (*)','GEV','Geometrica','Normal')
      print('ES Gumbel!!!')

   } else if (SminRMSE==6) {
      fit.distSEV <- fitgev.severity
      #severityOrdered<-severity[order(fitnorm_severityCDF)] 
      #severityOrderedCDF<-fitnorm_severityCDF[order(fitnorm_severityCDF)]
      severityCDF<-fitgev_severityCDF
      sev_distr <- pdfgev(severity_art,fitgev.severity)
      idajusteS<- 'GEV'
      BEST_FIT = c('NO','NO','NO','NO','NO','YES','NO','NO')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV (*)','Geometrica','Normal')
      print('ES GEV!!!')

   } else if (SminRMSE==7) {
      fit.distSEV <- fitgeom.severity
      #severityOrdered    <-severity[order(fitweibull_severityCDF)] 
      #severityOrderedCDF <-fitweibull_severityCDF[order(fitweibull_severityCDF)]
      severityCDF        <-fitgeom_severityCDF
      sev_distr <- dgeom(severity_art, prob = fitgeom.severity$estimate[1] , log = FALSE)
      idajusteS <- 'Weibull'
      BEST_FIT = c('NO','NO','NO','NO','NO','NO','YES','NO')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica (*)','Normal')
      print('CDF geometrico')


   } else if (SminRMSE==8) {
      fit.distSEV <- fitnorm.severity
      #severityOrdered    <-severity[order(fitweibull_severityCDF)] 
      #severityOrderedCDF <-fitweibull_severityCDF[order(fitweibull_severityCDF)]
      severityCDF        <-fitnorm_severityCDF
      sev_distr <- dnorm(severity_art, mean = fitnorm.severity$estimate[1], sd = fitnorm.severity$estimate[2], log = FALSE)
      idajusteS <- 'Normal'
      BEST_FIT = c('NO','NO','NO','NO','NO','NO','NO','YES')
      legendaS <- c('Observado','Gamma','Weibull', 'Exponencial', 'logNormal','Gumbel','GEV','Geometrica','Normal (*)')
      print('es Normal')

   }


   #SAVE in data frame severity erross
   df <- data.frame(distrs, RMSE_se,AIC_se ,BEST_FIT)
   write.table(df,file  = paste(ruta_salida,name_est,'_SEV_inf.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)

   #################Grafica de Histograma con LA DISTRIBUCION CON MEJOR AJUSTE########################
   maxis= max(sev_distr)+0.05
   maxid= max(dur_distr)+0.05
   #
   pdf(paste(ruta_salida,name_est ,'_histDURsev.pdf', sep = ''), width=9.0, height=4.5)
   par(mfrow = c(1, 2))
   hist(severity, prob=TRUE, ylim = c(0,maxis) ,xlab = 'Severidad', ylab='Densidad de Probabilidad', main = NA)
   lines(severity_art,sev_distr, col='red',lty = 1 ,lwd = 1.6)
   legend('topleft', legend= c( paste('Estación: ',name_est,sep = ''),paste('Distr.',idajusteS,sep = ' ')), col= c(NA,'red'), lty = c(NA, 1), lwd =   c(NA,1.6) ,box.lty=0, cex = 0.78)

   hist(duration, prob=TRUE,  ylim =c(0,maxid), xlab = 'Duración(Meses)' , ylab='Densidad de Probabilidad', main = NA)
   lines(duration_art,dur_distr, col='red',lty = 1 ,lwd = 1.6)
   legend('topleft', legend= c(paste('Estación: ',name_est,sep = ''),paste('Distr.',idajusteD,sep = ' ')), col= c(NA,'red'), lty = c(NA, 1), lwd =   c(NA,1.6) ,box.lty=0, cex = 0.78)
   dev.off()

   ##################PLOTS DURATION Empirical  EDF  vs CDF ##########################################
   Fnd<-ecdf(duration) 
   Fns<-ecdf(severity) 
   pdf(paste(ruta_salida,name_est ,'_plot_durSEV.pdf', sep = ''), width=9.0, height=4.0)
   par(mfrow = c(1, 2), mar=c(4.5,4.5,1.0,1.0 ))
   plot(sort(duration), sort(Fnd(duration)) , xlab = 'Duración (Meses)', ylab = 'Probabilidad Acumulada', 
       ylim=c(0, 1),type = 'p', col = 'red', pch = 4, cex=0.6)
   lines(sort(duration_art), sort(fitgamma_durCDFart),    col = 'black',    lty = 2 ,lwd = 0.8)
   lines(sort(duration_art), sort(fitweibull_durCDFart),  col = 'black', lty = 4, lwd = 0.8)
   lines(sort(duration_art), sort(fitexp_durCDFart),      col = 'black',   lty = 6, lwd = 0.8)
   lines(sort(duration_art), sort(fitlogn_durCDFart),     col = 'black' , lty = 5, lwd = 0.8)
   #lines(sort(duration_art), sort(fitnorm_durCDFart),     col = 'cyan4',  lty = 6, lwd = 1.3)
   lines(sort(duration_art), sort(fitgum_durCDFart),     col = 'black',  lty = 3, lwd = 1.2)
   lines(sort(duration_art), sort(fitgev_durCDFart),     col = 'black',  lty = 1, lwd = 0.8)

   lines(sort(duration_art), sort(fitgeom_durCDFart), type = 'o' ,pch = 1,   col = 'black',  lty = 1, lwd = 0.5, cex = 0.30)
   lines(sort(duration_art), sort(fitnorm_durCDFart), type = 'o' ,pch = 19  ,col = 'black',  lty = 1, lwd = 0.5, cex = 0.30)

   legend('bottomright', legend=legendaD, col= c('red','black','black','black', 'black','black','black','black','black'), lty = c(NA,2,4,6,5,3,1,1,1), pch  = c(4,NA,NA,NA,NA,NA,NA,1,19), 
            lwd =   c(NA,1.0, 1.0, 1.0, 1.0,1.3, 1.0,0.5, 0.5) ,box.lty=0, cex = 0.70)
   legend('topleft', legend= paste('Estación: ',name_est,sep = ''), col= 'black', cex = 0.70,box.lty=0)

   #################  PLOTS SEVERITY Empirical  EDF  vs CDF ##########################################
   plot(sort(severity), sort(Fns(severity)), xlab = 'Severidad', ylab = NA, 
        ylim=c(0, 1),type = 'p', col = 'red', pch = 4, cex = 0.6) 
   lines(sort(severity_art), sort(fitgamma_sevCDFart),   col = 'black', lty = 2, lwd = 0.80)
   lines(sort(severity_art), sort(fitweibull_sevCDFart), col = 'black',lty = 4, lwd = 0.80)
   lines(sort(severity_art), sort(fitexp_sevCDFart),     col = 'black',   lty = 6, lwd = 0.80)
   lines(sort(severity_art), sort(fitlogn_sevCDFart),    col = 'black',  lty = 5, lwd = 0.80)
   #lines(sort(severity_art), sort(fitnorm_sevCDFart),    col = 'cyan4', lty = 5, lwd = 1.2)
   lines(sort(severity_art), sort(fitgum_sevCDFart),    col = 'black', lty = 3, lwd = 1.2)
   lines(sort(severity_art), sort(fitgev_sevCDFart),    col = 'black', lty = 1, lwd = 0.80)

   lines(sort(severity_art), sort(fitgeom_sevCDFart), type = 'o' ,pch = 1,   col = 'black',  lty = 1, lwd = 0.5, cex = 0.30)
   lines(sort(severity_art), sort(fitnorm_sevCDFart), type = 'o' ,pch = 19  ,col = 'black',  lty = 1, lwd = 0.5, cex = 0.30)

   legend('bottomright', legend=legendaS, col= c('red','black','black' , 'black', 'black','black','black','black','black') , lty = c(NA,2,4,6,5,3,1,1,1), 
           lwd = c(NA,1.0, 1.0, 1.0, 1.0,1.3, 1.0,0.5, 0.5) ,box.lty=0 , pch  =c(4,NA,NA,NA,NA,NA,NA,1,19),cex = 0.70)
    #legend('topleft', legend= paste('Estación: ',name_est,sep = ''), col= 'black', cex = 0.76,box.lty=0)
    dev.off()


   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   #%HASTA Aqui se tiene 3 resultados fundamentales:severityOrdered ,severityOrderedCDF,severityCDF  %
   #%durationOrdered  ,durationOrderedCDF ,durationCDF:resultados fundamentales para las copulas     %
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ###########################################################################################################################################
   #################################################COPULA ESTIMATION ########################################################################
   ###########################################################################################################################################

   #el mejor CDFtanto severidad como duracion.
   concatenatedDataCDF<-cbind(durationCDF,severityCDF)  
   concatenatedData   <-cbind(matrix(duration,ncol = 1),matrix(severity,ncol = 1))
   copula.empirical  <- C.n(concatenatedDataCDF, concatenatedData,smoothing = "beta") #beta para datos pequenios

   ############################ 
   ###### Check best copula ###
   ############################
   #ITAU method
   #copula.cross<-list() #create a empty list
   #k <- 30
   #copula.cross[1]<-xvCopula(frankCopula(dim=2), x = matrix(concatenatedDataCDF,ncol = 2),k = k   )
   #copula.cross[2]<-xvCopula(gumbelCopula(dim=2),x = matrix(concatenatedDataCDF,ncol = 2),k = k   )
   #copula.cross[3]<-xvCopula(tCopula(dim = 2),    x = matrix(concatenatedDataCDF,ncol = 2),k = k  )
   #copula.cross[4]<-xvCopula(normalCopula(dim=2),x = matrix(concatenatedDataCDF,ncol = 2),k = k   )
   #copula.cross[5]<-xvCopula(joeCopula(dim=2),   x = matrix(concatenatedDataCDF,ncol = 2),k = k   )
   #copula.cross[6]<-xvCopula(claytonCopula(dim=2), x = matrix(concatenatedDataCDF,ncol = 2),k = k)
   #copula.cross[7]<-xvCopula(plackettCopula(),     x = matrix(concatenatedDataCDF,ncol = 2),k = k )
   #copula.cross[8]<-xvCopula(galambosCopula(), x =  matrix(concatenatedDataCDF, ncol =2),k= k )
   #copula.cross[1]<-xvCopula(frankCopula(dim=2),matrix(concatenatedDataCDF,ncol = 2),k = NULL,verbose = interactive(),method=c("itau"),optim.control = list(maxit=5000,deps=1e-4))

   ###################FIT COPULA method maximun log-likelihood########################### #

   cop.fitF<-fitCopula(copula=frankCopula(dim = 2), data=concatenatedDataCDF,
                     method = c("ml"),posDef = is(frankCopula(dim = 2), "ellipCopula"),
                     start = NULL, lower = NULL, upper = NULL,
                     optim.method = optimMeth(frankCopula(dim = 2), method, dim = d),
                     optim.control = list(maxit=9000),
                     estimate.variance = NA, hideWarnings = FALSE)
   paramFitF<-as.numeric(coef(cop.fitF))

   cop.fitG<-fitCopula(copula=gumbelCopula(dim = 2), concatenatedDataCDF,
                     method = c("ml"),posDef = is(gumbelCopula(dim = 2), "ellipCopula"),
                     start = NULL, lower = NULL, upper = NULL,
                     optim.method = optimMeth(gumbelCopula(dim = 2), method, dim = d),
                     optim.control = list(maxit=9000),
                     estimate.variance = NA, hideWarnings = FALSE)
   paramFitG<-as.numeric(coef(cop.fitG))
  
   cop.fitT<-fitCopula(copula=tCopula(dim = 2), concatenatedDataCDF,
                    method = c("ml"),posDef = is(tCopula(dim = 2), "ellipCopula"),
                     start = NULL, lower = NULL, upper = NULL,
                     optim.method=optimMeth(gumbelCopula(dim = 2), method, dim = d),
                     optim.control = list(maxit=3000),
                     estimate.variance = NA, hideWarnings = FALSE)
   paramFitT<-as.numeric(coef(cop.fitT))
   #"Nelder-Mead"
   cop.fitN<-fitCopula(copula=normalCopula(dim = 2), concatenatedDataCDF,
                     method = c("ml"),posDef = is(normalCopula(dim = 2), "ellipCopula"),
                     start = NULL, lower = NULL, upper = NULL,
                     optim.method = optimMeth(normalCopula(dim = 2), method, dim = d),
                     optim.control = list(maxit=9000),
                     estimate.variance = NA, hideWarnings = FALSE)
   paramFitN<-as.numeric(coef(cop.fitN))  
  
   cop.fitJ<-fitCopula(copula=joeCopula(dim = 2), concatenatedDataCDF,
                   method = c("ml"),posDef = is(joeCopula(dim = 2), "ellipCopula"),
                   start =  NULL, lower = NULL, upper = NULL,
                   optim.method = optimMeth(joeCopula(dim = 2), method = 'ml', dim = d),
                   optim.control = list(maxit=10000),
                   estimate.variance = NA, hideWarnings = FALSE)
   paramFitJ<-as.numeric(coef(cop.fitJ))

   cop.fitC<-fitCopula(copula=claytonCopula(dim = 2), concatenatedDataCDF,
                   method = c("ml"),posDef = is(claytonCopula(dim = 2), "ellipCopula"),
                   start = NULL, lower = 0.0001, upper = 30.0,
                   optim.method ="Brent",
                   optim.control = list(maxit=20000),
                   estimate.variance = NA, hideWarnings = FALSE)
   paramFitC<-as.numeric(coef(cop.fitC))

   cop.fitP<-fitCopula(copula=plackettCopula(), concatenatedDataCDF,
                   method = c("ml"),posDef = is(plackettCopula(), "ellipCopula"),
                   start = NULL, lower = NULL, upper = NULL,
                   optim.method = optimMeth(plackettCopula(), method, dim = d),
                   optim.control = list(maxit=9000),
                   estimate.variance = NA, hideWarnings = FALSE)
   paramFitP<-as.numeric(coef(cop.fitP))

   cop.fitGl<-fitCopula(copula=galambosCopula(), concatenatedDataCDF,
                   method = c("ml"),posDef = is(galambosCopula(), "ellipCopula"),
                   start = NULL, lower = NULL, upper = NULL,
                   optim.method = optimMeth(galambosCopula(), method, dim = d),
                   optim.control = list(maxit=9000),
                   estimate.variance = NA, hideWarnings = FALSE)
   paramFitGl<-as.numeric(coef(cop.fitGl))


   ########################## Copula result from Observations #############################################
   nCopulaObs<-pCopula(concatenatedDataCDF,copula = normalCopula(param=paramFitN[1])) 
   frankCopulaObs<-pCopula(concatenatedDataCDF,copula = frankCopula(param=paramFitF))
   gumbelCopulaObs<-pCopula(concatenatedDataCDF,copula = gumbelCopula(param=paramFitG))
   tCopulaObs<-pCopula(concatenatedDataCDF,copula = tCopula(param=paramFitT[1],df=as.integer(paramFitT[2])))
   joeCopulaObs<-pCopula(concatenatedDataCDF,copula = joeCopula(param=paramFitJ))
   claytonCopulaObs<-pCopula(concatenatedDataCDF,copula = claytonCopula(param=paramFitC[1]))
   plackettCopulaObs<-pCopula(concatenatedDataCDF,copula = plackettCopula(param=paramFitP[1]))
   galambosCopulaObs<-pCopula(concatenatedDataCDF,copula =galambosCopula(param=paramFitP[1]))

   #################################################################
   ############### Which Copula is the best suited #################
   #################################################################

   #QQ-plot betweeen empirical and fitted distributions
   # qqplot(copula.empirical,frankCopulaObs,plot.it = TRUE, xlab = "Empirical",
   #ylab = "Fitted")

   rmse=c(numeric(8))
   rmse[1]<-rmse(nCopulaObs,copula.empirical,       na.rm = TRUE)
   rmse[2]<-rmse(frankCopulaObs,copula.empirical,   na.rm = TRUE)
   rmse[3]<-rmse(gumbelCopulaObs,copula.empirical,  na.rm = TRUE)
   rmse[4]<-rmse(tCopulaObs,copula.empirical,       na.rm = TRUE)
   rmse[5]<-rmse(joeCopulaObs,copula.empirical,     na.rm = TRUE)
   rmse[6]<-rmse(claytonCopulaObs,copula.empirical, na.rm = TRUE)
   rmse[7]<-rmse(plackettCopulaObs,copula.empirical,na.rm = TRUE)
   rmse[8]<-rmse(galambosCopulaObs,copula.empirical,na.rm = TRUE)

   #CORRELACIONES DE EMPIRICO Y AJUSTADA
   corcop =c(numeric(8))
   corcop[1]<-cor(nCopulaObs,copula.empirical,        method = 'pearson')
   corcop[2]<-cor(frankCopulaObs,copula.empirical,    method = 'pearson')
   corcop[3]<-cor(gumbelCopulaObs,copula.empirical,   method = 'pearson')
   corcop[4]<-cor(tCopulaObs,copula.empirical,        method = 'pearson')
   corcop[5]<-cor(joeCopulaObs,copula.empirical,      method = 'pearson')
   corcop[6]<-cor(claytonCopulaObs,copula.empirical,  method = 'pearson')
   corcop[7]<-cor(plackettCopulaObs,copula.empirical, method = 'pearson')
   corcop[8]<-cor(galambosCopulaObs,copula.empirical, method = 'pearson')


   LLabFrank <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = frankCopula(param=pars),log = TRUE))-2*1)
     }
  
   LLabGumbel <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = gumbelCopula(param=pars),log = TRUE)) - 2*1)
     }
   
   LLabtStudent <- function(pars,X){
      #log likelihood for Gamma(alpha,scale=beta)
      #X is data
     return(2*sum(dCopula(X,copula = tCopula(param=pars[1],df=as.integer(pars[2])),log = TRUE))- 2*1)
     }
  
   LLabNormal <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = normalCopula(param=pars[1]),log = TRUE)) - 2*1)
     }

   LLabJoe <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = joeCopula(param=pars[1]),log = TRUE)) - 2*1)
     }  

   LLabClayton <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = claytonCopula(param=pars[1]),log = TRUE)) - 2*1)
     }  

   LLabPlackett <- function(pars,X){
    # log likelihood for Gamma(alpha,scale=beta)
    # X is data
    return(2*sum(dCopula(X,copula = plackettCopula(param=pars[1]),log = TRUE)) - 2*1)
    }  

   LLabGalambos <- function(pars,X){
     # log likelihood for Gamma(alpha,scale=beta)
     # X is data
     return(2*sum(dCopula(X,copula = galambosCopula(param=pars[1]),log = TRUE)) - 2*1)
     } 

   maxLL<-c(numeric(8))
   maxLL[1]=abs(LLabNormal(paramFitN,concatenatedDataCDF))
   maxLL[2]=abs(LLabFrank(paramFitF,concatenatedDataCDF))
   maxLL[3]=abs(LLabGumbel(paramFitG,concatenatedDataCDF))
   maxLL[4]=abs(LLabtStudent(paramFitT,concatenatedDataCDF))
   maxLL[5]=abs(LLabJoe(paramFitJ,concatenatedDataCDF))
   maxLL[6]=abs(LLabClayton(paramFitC,concatenatedDataCDF))
   maxLL[7]=abs(LLabPlackett(paramFitP,concatenatedDataCDF))
   maxLL[8]=abs(LLabGalambos(paramFitP,concatenatedDataCDF))

   ## Generate Synthetic Values of CDF for Duration and Severity , Valores###################################################################333333
   list_lim <- c(max(duration), max(severity))    
   lim <- max(list_lim) + 1

   x<-c(seq(0.0, round(lim), length.out=90)) #0.05
   y<-c(seq(0.0, round(lim), length.out=90))
   z<-meshgrid(x,y)

   durationInverse<-as.vector(z$x)
   severityInverse<-as.vector(z$y)

   #Create vectors
   if (SminRMSE==1){
     severityCDFCop<-pgamma( as.vector(z$y),shape = fitgamma.severity$estimate[1], rate = fitgamma.severity$estimate[2], lower.tail = TRUE,log.p = FALSE)
   } else if (SminRMSE==2){
     severityCDFCop<-pweibull( as.vector(z$y), shape = fitweibull.severity$estimate[1], scale = fitweibull.severity$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }else if (SminRMSE==3){
     severityCDFCop<-pexp( as.vector(z$y),rate = fitexp.severity$estimate[1], lower.tail = TRUE,log.p = FALSE)
   }else if (SminRMSE==4){ 
     severityCDFCop<-plnorm( as.vector(z$y), meanlog = fitlogn.severity$estimate[1], sdlog = fitweibull.severity$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }else if (SminRMSE==5){
     severityCDFCop<-pgumbel( as.vector(z$y), fitgum.severity$estimate[1], fitgum.severity$estimate[2])
   }else if (SminRMSE==6){
     severityCDFCop<-cdfgev(as.vector(z$y),fitgev.severity)
   }else if (SminRMSE==7){
     severityCDFCop<-pgeom( as.vector(z$y), prob = fitgeom.severity$estimate[1], lower.tail = TRUE,log.p = FALSE)
   }else if (SminRMSE==8){
     severityCDFCop<-pnorm( as.vector(z$y), mean = fitnorm.severity$estimate[1], sd = fitnorm.severity$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }


   if (DminRMSE==1){
     durationCDFCop<-pgamma(as.vector(z$x),shape = fitgamma.duration$estimate[1], rate = fitgamma.duration$estimate[2], lower.tail = TRUE,log.p = FALSE)
   } else if (DminRMSE==2){
     durationCDFCop<-pweibull(as.vector(z$x), shape = fitweibull.duration$estimate[1], scale = fitweibull.duration$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }else if (DminRMSE==3){
     durationCDFCop<-pexp( as.vector(z$x),    rate = fitexp.duration$estimate[1], lower.tail = TRUE,log.p = FALSE)
   }else if (DminRMSE==4){
     durationCDFCop<-pnorm(as.vector(z$x), mean = fitnorm.duration$estimate[1], sd = fitnorm.duration$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }else if (DminRMSE==5){
     durationCDFCop<-pgumbel( as.vector(z$x), fitgum.duration$estimate[1], fitgum.duration$estimate[2])
   }else if (DminRMSE==6){
     durationCDFCop<-cdfgev(as.vector(z$x),fitgev.duration)
   }else if (DminRMSE==7){
     durationCDFCop<-pgeom(as.vector(z$x), prob = fitgeom.duration$estimate[1], lower.tail = TRUE,log.p = FALSE)
   }else if (DminRMSE==8){
     durationCDFCop<-pnorm(as.vector(z$x), mean = fitnorm.duration$estimate[1], sd = fitnorm.duration$estimate[2], lower.tail = TRUE,log.p = FALSE)
   }
   
   ##################################################################################
   ####################### Select better Copula distribution  #######################
   ##################################################################################

   minRMSE<-which.min(rmse)
   getInf <- is.infinite(maxLL)
   maxLL[getInf]=NaN
   maxLL.max<-which.min(maxLL)
 
   maxcor <- which.max(corcop)

   ## A combination of both metrics is created to choose the most convinient distributions #####################
   indexSelected <- which.min(abs( (rmse-rmse[minRMSE])/rmse[minRMSE]) +abs( (maxLL-maxLL[maxLL.max])/maxLL[maxLL.max])     )
   #indexSelected2 <- which.min(copula.cross)

   print(paste('El indice de mejor ajuste no cruzada es : ',as.character(indexSelected[1]), sep = ''))
   #print(paste('El indice de mejor ajuste cruzada    es : ',as.character(indexSelected2), sep = ''))
   ###########################################################################################
   ####################### Generate probabilities with the selected copula #################
   ###########################################################################################

   ZZ<-cbind(as.vector(durationCDFCop),as.vector(severityCDFCop))
   ZZ2<-cbind(as.vector(durationCDF),as.vector(severityCDF))

   if (indexSelected==1){
     copulaGen     <-pCopula(ZZ,copula = normalCopula(param=paramFitN))
     empirical.cop <-pCopula(ZZ2,copula = normalCopula(param=paramFitN))
     cop.fit <- cop.fitN
     BEST_FIT = c('YES','NO','NO','NO', 'NO', 'NO', 'NO','NO')
     name_cop <- 'Normal Copula'
   
   }else if(indexSelected==2){
     copulaGen<-pCopula(ZZ,copula = frankCopula(param=paramFitF))
     empirical.cop  <-pCopula(ZZ2,copula = frankCopula(param=paramFitF))
     cop.fit <- cop.fitF
     BEST_FIT = c('NO','YES','NO','NO', 'NO', 'NO', 'NO', 'NO')
     name_cop <- 'frank Copula'

   }else if(indexSelected==3){
     copulaGen     <-pCopula(ZZ,copula = gumbelCopula(param=paramFitG))
     empirical.cop <-pCopula(ZZ2,copula = gumbelCopula(param=paramFitG))
     cop.fit       <- cop.fitG
     BEST_FIT = c('NO','NO','YES','NO', 'NO', 'NO', 'NO','NO')
     name_cop <- 'Gumbel Copula'

   }else if(indexSelected==4){
     copulaGen<-pCopula(ZZ,copula = tCopula(param=paramFitT))
     empirical.cop <-pCopula(ZZ2,copula = tCopula(param=paramFitT))
     cop.fit <- cop.fitT
     BEST_FIT = c('NO','NO','NO','YES', 'NO', 'NO', 'NO', 'NO')
     name_cop <- 'T. Copula' 

   }else if(indexSelected==5){
     copulaGen<-pCopula(ZZ,copula = joeCopula(param=paramFitJ))
     empirical.cop<-pCopula(ZZ2,copula = joeCopula(param=paramFitJ))
     cop.fit <- cop.fitJ
     name_cop <- 'Joe Copula'
     BEST_FIT = c('NO','NO','NO','NO', 'YES', 'NO', 'NO','NO')

   }else if(indexSelected==6){
     copulaGen<-pCopula(ZZ,copula = claytonCopula(param=paramFitC))
     empirical.cop <-pCopula(ZZ2,copula = claytonCopula(param=paramFitC))
     cop.fit <- cop.fitC
     name_cop <- 'Clayton Copula'
     BEST_FIT = c('NO','NO','NO','NO', 'NO', 'YES', 'NO','NO')

   }else if(indexSelected==7){
     copulaGen<-pCopula(ZZ,copula = plackettCopula(param=paramFitP))
     empirical.cop <- pCopula(ZZ2,copula = plackettCopula(param=paramFitP))
     cop.fit <- cop.fitP
     name_cop <- 'Plackett Copula'
     BEST_FIT = c('NO','NO','NO','NO', 'NO', 'NO', 'YES','NO')

   }else if(indexSelected==8){
     copulaGen<-pCopula(ZZ,copula = galambosCopula(param=paramFitGl))
     empirical.cop <- pCopula(ZZ2,copula = galambosCopula(param=paramFitGl))
     cop.fit <- cop.fitGl
     name_cop <- 'Galambos Copula'
     BEST_FIT = c('NO','NO','NO','NO', 'NO', 'NO', 'NO','YES')
   }

   #copula.empirical  #metodo empirico
   #empirical.cop     #metodo ajustado

   ###### PRINT FUNCIONES DE LOSS333METRICAS ##################3333
   nam_copulas <- c('Normal Cop.','Frank Cop.', 'Gumbel Cop.', 't Cop.', 'Joe Cop.', 'Clayton Cop.', 'Plackett Cop.','Galambos Cop.')
   #BEST_FIT, rmse, maxLL, copula.cross
   #copula_cross <- unlist(as.vector(matrix(copula.cross, )))
   dursev.kendall <- rep(cor_kendall, times = 8 )
   dursev.pearson <- rep(cor_pearson, times = 8 )

   df_copul <- data.frame(nam_copulas, rmse, maxLL, corcop , BEST_FIT, dursev.pearson, dursev.kendall )
   write.table(df_copul,file  = paste(ruta_salida,name_est,'_COP_inf.txt', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)

   rm(df_copul)
   rm(dursev.pearson)
   rm(dursev.kendall)

   #concatenatedDataCDF
   ########################PLOT BEST FITT AND empirical ######################################################################
   ########################PLOT BEST FITT AND empirical ######################################################################
   #make linear regression with empirical copula and Teoretical copula
   #copula.empirical  #metodo empirico
   #empirical.cop     #metodo ajustado
   df_emp <- data.frame(x = copula.empirical , y = empirical.cop)
   my.formula <- y ~ x
   p <- ggplot(data = df_emp, aes(x = x, y = y)) +theme_classic()+
     geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula, lwd = 0.7) +
     stat_poly_eq(formula = my.formula,
                eq.with.lhs = "italic(TC)~`=`~ ", eq.x.rhs = "italic(EC)",
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE) +   labs(x = 'Copula Empirica (EC)', y = paste(name_cop,' Teorica',' (TC) ', sep='' ))      +
     geom_point(color='black', pch = 4, cex = 1.6)

   #copula.empirical=  copula empirica beta
   corCopulas <- cor(copula.empirical, empirical.cop, method = 'pearson')
   fit        <- lm( empirical.cop  ~  copula.empirical)
   regr_sum   <- summary(fit)
   coef_fit   <- round(regr_sum[[9]],3)
   rmse_cop <- round(rmse(copula.empirical, empirical.cop),4)

   #corCopulas coef_fit 
   
   
   #plot(empirical.cop, copula.empirical, type = 'p', pch=4, cex = 1.2, xlab = 'Copula (Empirica)', ylab='Copula galambos (teórica) ')
   #abline(fit,  lwd =1.3, col ='blue')
   #legend('bottomright', legend=c(paste('CC = ', as.character(coef), sep = '')), col= c( NA) , lty = c(NA), lwd = c(NA) ,box.lty=0 , pch =c(NA),cex = c(NA))
   #make plots empiricas cop and teoretical
   cpFra <- contourplot2(cop.fit@copula, FUN = pCopula, region=FALSE, key= list(corner=c(0.04,0.04),lines = list(col=c(1,2), lwd= 1.9, lty = c(1,2)), 
              text = list(c(paste(name_cop, '  Ajustada', sep =''), 'Copula  Empírico')) ), lwd =1.6)
   #compute the ocntours of the empirical Copula
   U <- as.matrix(concatenatedDataCDF)
   u <- seq(0,1,length.out = 16)
   grid <- as.matrix(expand.grid (u1=u, u2 = u))
   val <- cbind(grid,z = C.n(grid, X=U, smoothing = 'beta'))
   cpCn <- contourplot2(val, region = FALSE , labels =FALSE, col = 2, lwd = 1.7, lty = 2)

   cop_finPLOT <- cpFra + cpCn

   pdf(paste(ruta_salida,name_est ,'_empiricos.pdf', sep = ''), width= 10.0, height= 5.0)
   plot(ggarrange(p, cop_finPLOT , ncol=2, nrow=1))
   dev.off()

   ###########################################################################################################################
   #N = length(as.vector(z$x))
   #cmap = rich.colors(N)
   #par(mfrow = c(1,1))
   #scatterplot3d(ZZ[,1],ZZ[,2],copulaGen, color = cmap, cex.symbols = 1.2, pch=19, type = "b")  
   #scatterplot3d(z,z,copulaGen, color = cmap, cex.symbols = 1.2, pch=19, type = "b")
   #plot(x,copulaGen)


   #CREAMOS UN DATA FRAME DE DURATION, SEVERITY AND COPULA
   #df.fitt <- data.frame(durationCDFCop,severityCDFCop,copulaGen)

   ##################################################################################################
   ########################### NOW ESTIMATE JOINT PROBABILITIES PARA PLOTS  #########################
   ##################################################################################################
   P = 1 + copulaGen - ZZ[,1]- ZZ[,2] ## Including marginal probabilities
   P_2 = 1 - copulaGen                ## Without marginal probabilities

   ## Compute probabilities----------------------------------------------------
   idxOk    <- which(abs(durationCDFCop)<0.003| abs(severityCDFCop)<0.003) #el y 
   P[idxOk] = NaN

   durationCDFCop[idxOk]=NaN
   severityCDFCop[idxOk]=NaN

   idxOk2 <-which(abs(P)<0.003) #el o
   durationCDFCop[idxOk2] = NaN
   severityCDFCop[idxOk2] = NaN
   P[idxOk2]   =  NaN
   P_2[idxOk]  =  NaN
   P_2[idxOk2] =  NaN

   ## Compute Return Periods (in years)  
   RT<-(mean_interar/12.)/P       #el y
   RT2<-(mean_interar/12.)/P_2    #el ó

   ## Individual Return Periods CALCULO 1
   RT_duration<-(mean_interar/12.)/(1-durationCDF)    #ajustada
   RT_severity<-(mean_interar/12.)/(1-severityCDF)    #ajustada

	   ## Conditional Return Periods
   RT_condSeverity<-RT*(1/(1-severityCDFCop))
   RT_condDuration<-RT*(1/(1-durationCDFCop))

   #data frame para plot
   df_RTcop  <- data.frame(dur_cop = as.vector(z$x), sev_cop = as.vector(z$y), RT)
   df_RTcop <- na.omit(df_RTcop)

   #dataframe for plot curves  
   df_RT2cop  <- data.frame(dur_cop = as.vector(z$x), sev_cop = as.vector(z$y), RT2)
   df_RT2cop <- na.omit(df_RT2cop)

   duraAR <- as.vector(z$x)
   seveAR <- as.vector(z$y)
   condDur <- RT_condDuration
   condSev <- RT_condSeverity
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   RT_plot <- contourplot2(df_RTcop, xlab = 'Duración (meses)', ylab='Severidad',lwd=0.9,col="black",region= FALSE) #36, 500 
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   RT2_plot <- contourplot2(df_RT2cop, xlab = 'Duración (meses)', ylab='Severidad' , cex=0.3 , lwd=0.9,
            col="black",region= FALSE) #cuts = 12

   pdf(paste(ruta_salida,name_est ,'_RTS.pdf', sep = ''), width= 10.0, height= 5.0)
   plot(ggarrange(RT_plot,RT2_plot, ncol=2, nrow=1))
   dev.off()
   #, labels = c(name_est)
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT FINAL DE COPULAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   # Graficar densidad conjunta de (X,Y) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   #for (i in 1:length(xx)){
     #for (j in 1:length(yy)){
     # zz[i, j] <- pCopula(c(xx[i], yy[j]),copula = claytonCopula(param=paramFitC))
     # }
     #}

    #zz_RT  <- matrix(RT, nrow = length(x), ncol = length(y))
    #zz_RT2 <- matrix(RT2, nrow = length(x), ncol = length(y))
    #xli <- c(min(df_RTcop$dur_cop), max(df_RTcop$dur_cop))
    #yli <- c(min(df_RTcop$sev_cop), max(df_RTcop$sev_cop))
    #contour(x, y, zz, xlim =xli,ylim=yli ,zlim = c(0, 20),nlevels = 20, main = "Curvas de nivel", xlab = "X", ylab = "Y")
    #pdf(paste(ruta_salida,name_est ,'_final_contour.pdf', sep = ''),width = 10.0, height = 5.0 )
    #par(mfrow = c(1, 2))
    #contour(x, y, zz_RT ,xlim= c(1.5, 9.0),ylim = c(3.0,9.0), main = "Curvas de nivel",  xlab = 'Duración (meses)', ylab = "Severidad", nlevels = 20 ,xaxs='i',yaxs='i')
    #contour(y, x, zz_RT2,zlim = c(0, 50),nlevels = 9, main = "Curvas de nivel",  xlab = 'Duración (meses)', ylab = "Severidad", xaxs='i',yaxs='i')
    #persp(x, y, zz, main = "Densidad de (X,Y)", xlab = "X", ylab = "Y",     zlab = "f(x,y)", theta = -35, phi = 30, col = "yellow")
    #dev.off()

   #####################################################################################################
   ########################### NOW ESTIMATE JOINT PROBABILITIES PARA EMPIRICOS #########################
   #####################################################################################################

   #ecuaciones para el calculo de cop return periods=========================
   P2 = 1 + empirical.cop - ZZ2[,1]- ZZ2[,2] ## Including marginal probabilities

   P_22 = 1 - empirical.cop                  ## Without marginal probabilities

   ## Compute probabilities
   idxOk2<- which(abs(durationCDF)<0.001| abs(severityCDF) < 0.001) #el y
   P2[idxOk2] = NaN

   durationCDF[idxOk2]=NaN
   severityCDF[idxOk2]=NaN

   idxOk22<-which(abs(P2) < 0.001) #el o
   durationCDFCop[idxOk22] = NaN
   severityCDFCop[idxOk22] = NaN
   P2[idxOk22] = NaN
   P_22[idxOk2] = NaN
   P_22[idxOk22] = NaN

   ## Compute Return Periods (in years)  
   RT<-(mean_interar/12.)/P2       #el y
   RT2<-(mean_interar/12.)/P_22    #el ó

   ## Individual Return Periods CALCULO 1
   RT_duration<-(mean_interar/12.)/(1-durationCDF)    #ajustada
   RT_severity<-(mean_interar/12.)/(1-severityCDF)    #ajustada

   ## Conditional Return Periods
   RT_condSeverity<-RT*(1/(1-severityCDF))
   RT_condDuration<-RT*(1/(1-durationCDF))

   #AQUI SI QUIERO EL PLOT AGO SI NO ON
   nill <- length(RT)
   null <- length(intera_time)

   if(nill == null){
     intera_timeFUll <- c(intera_time)
   }else if(nill > null){
     intera_timeFUll <- c(intera_time, NA)
   }

   dg<- data.frame(duration, severity,durationCDF, severityCDF, RT_duration, RT_severity,P2, RT , P_22 , RT2 , RT_condDuration, RT_condSeverity,  intera_timeFUll  )
   write.table(dg, file = paste(ruta_salida,name_est, '_FINAL_table.csv', sep = ''), row.names = FALSE,  col.names= TRUE, sep ='\t')

   #################################################################APPEND IN A DATAFRAME AND SAVE individual################################################
   #df_rt <- data.frame(duration, RT_duration,  severity , RT_severity)
   #estadistica basica
   #dur_unique           <- data.frame(dur =  df_rt$duration , dur_RT  = df_rt$RT_duration)
   #dur_unique['dur_RT'] <- round(dur_unique$dur_RT,  1)
   #dur_unique         <- unique(dur_unique)
   #dur_unique         <- dur_unique[order(dur_unique$dur),]
   #write.table(dur_unique, file = paste(ruta_salida,name_est, '_dur_uniq.csv', sep = ''), row.names = FALSE,  col.names= TRUE, sep ='\t')
   ##### 
   #sev_unique         <- data.frame(sev =  df_rt$severity , sev_RT  = df_rt$RT_severity)
   #sev_unique['sev_RT'] <- round(sev_unique$sev_RT,  1)
   #sev_unique['sev'] <- round(sev_unique$sev   ,  1)
   #sev_unique         <- unique(sev_unique)
   #sev_unique         <- sev_unique[order(sev_unique$sev),]
   #write.table(sev_unique, file = paste(ruta_salida, name_est , '_sev_uniq.csv', sep = ''), row.names = FALSE,  col.names= TRUE, sep ='\t')

   #rm(dur_unique)
   #rm(sev_unique)



   #ESCRIBIR ERRORES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   dur_out <- capture.output(fit.distDUR)
   sev_out <- capture.output(fit.distSEV)
   cop_out <- capture.output(summary(cop.fit))

   cat(paste("Resultados de ajuste Duracion-------------", idajusteD, sep =' '), dur_out , ' ' ,paste('Resultados de ajuste severidad---------------', idajusteS, sep = ' '),
             sev_out ,     file=paste(ruta_salida,name_est, '_summary.txt', sep = ''), sep="\n", append=FALSE, fill = TRUE)
   cat(' ' ,paste('Resultados de copula--------------:',name_cop, sep =''), cop_out , file=paste(ruta_salida,name_est, '_summary.txt', sep = ''), sep="\n", append=TRUE, fill = TRUE)


   #alamacenando resultados finales
   c_ests      <- append(c_ests     , name_est)
   c_corrk     <- append(c_corrk    , cor_ken)
   c_corrp     <- append(c_corrp    , cor_per)
   c_corrsp    <- append(c_corrsp   , cor_sper)
   c_marg_D    <- append(c_marg_D   , idajusteD)
   c_marg_S    <- append(c_marg_S   , idajusteS)
   c_copula    <- append(c_copula   , name_cop)
   c_pearson   <- append(c_pearson  , corCopulas)  
   c_fitline   <- append(c_fitline  , coef_fit ) 
   c_rmsecop   <- append(c_rmsecop , rmse_cop)
   c_msevRad   <- append(c_msevRad  , m)
   c_msevGrade <- append(c_msevGrade,round(m_grade,3))

}

#IN DATA FRAME AND PRINT RESULTS------
df <- data.frame(c_ests,c_msevRad,c_msevGrade,c_marg_D, c_marg_S, c_corrk, c_corrp, c_corrsp ,c_copula, c_pearson , c_fitline, c_rmsecop )
write.table(df,file  = paste(ruta_salida0,region,'_info.scv', sep = ''), sep = '\t', row.names = FALSE, col.names = TRUE)

print('###############################################')
print('#              DONE RIGHT ThinkPad            #')
print('###############################################')

#http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots  

#contour(x = xgrid, y = ygrid, z = mtrx3d, zlim = c(14,16), nlevels = 4,ylab = 'Horsepower', labcex = 1.2, )


##################################FINAL PLOT################################################
#x_dur <-c(seq(0.01,10.0,length.out = 200)) #0.05
#y_sev <-c(seq(0.01,10.0,length.out = 200))
#simPredMat <- matrix(RT, 200, 200)   #length(simPred)=1000000
#simDF <- data.frame(simDur = x, simMag = y, simPredMat)

#library(evd)


#plot.new()
#pdf('final_estoyFELIZ.pdf')
#contour(simDF$simDur, simDF$simMag, simPredMat, xaxs = 'i', yaxs = 'i', labcex = 0.6, lwd = 1, col = "black", method = "flattest", vfont = c("sans serif", "plain"))
#dev.off()

