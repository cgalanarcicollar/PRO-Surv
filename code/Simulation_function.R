simulation.function <-  function(nSim, nRand,nMeas, nTrial, beta, unorm,vnorm,phi, sigmau,sigmav, weibull,alpha,cens,dat){
 
   proportion <- function(b, se, true, level = .95, df = Inf){
    # Estimate, standard error, confidence level, true parameter, and df
    qtile <- level + (1 - level)/2 # Compute the proper quantile
    lower.bound <- b - qt(qtile, df = df)*se # Lower bound
    upper.bound <- b + qt(qtile, df = df)*se # Upper bound
    # Is the true parameter in the confidence interval? (yes = 1)
    in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
    cp <- mean(in.ci) # The coverage probability
    return(list(proportion = cp, # Return results
                in.ci = in.ci,
                ci = cbind(lower.bound, upper.bound)))
  }
  options(warn=2)
  
  # long response
  y.all<-c()
  
  #times
  e.time.all<-c()
  cens.all<-c()
  status.all<-c()
  obstime.all<-c()
  
  #ALPHA
  TScoeff.all<-c()
  JMcoeff.all<-c()
  TVCMcoeff.all<-c()
  
  TSse.all<-c()
  TVCMse.all<-c()
  JMse.all<-c()
  
  significance.TS<-0
  significance.TVCM<-0
  significance.JM<-0
  no.conv.JM<-0
  
  #LONGITUDINAL BB
  coef.bb.all<-c()
  se.bb.all <- c()
  rand.bb.all<-c()
  sigmau.bb.all<-c()
  se.bb.sigmau<-c()
  sigmav.bb.all<-c()
  se.bb.sigmav<-c()
  psi.bb.all<-c()
  se.psi<-c()
  
  #LONGITUDINAL NORMAL
  coef.norm.all<-c()
  se.norm.all<-c()
  sigmau.norm.all<-c()
  se.norm.sigmau<-c()
  sigmav.norm.all<-c()
  se.norm.sigmav<-c()
  epsilon.norm.all<-c()
  se.epsilon<-c()
  
  # SET 
  nObs<-nMeas*nRand #Numero de mediciones totales
  dat.long<- tidyr::pivot_longer(dat, cols=2:5, names_to = "meas",values_to = "time" )
  X <- model.matrix(~time,data=dat.long) #matrices del modelo longitudinal
  Z.i <- model.matrix(~z-1,data=dat.long)
  Z.t<- model.matrix(~z:time-1,data=dat.long)
  Z<-cbind(Z.i,Z.t)
  p<-1/(1+exp(-(X%*%beta+Z%*%c(unorm,vnorm))))
  
  n<-0
  ##BASELINE SURVIVAL ELEMENTS
  for (i in 1:nSim){
    no.conv <- 1
    
    while(no.conv==1){
      
      
      no.conv<-0
      
      
      # GENERATE LONGITUDINAL RESPONSES
      y<-PROreg::rBB(nObs,nTrial,p,phi)
      dat.y<-cbind(dat.long,y)
      
      # GENERATE TIME-TO-EVENT
      
      e.time<-c()
      for(j in 1:nRand){
        h<-function(s){
          h <- weibull[1]*weibull[2]*s^(weibull[2]-1)*exp(alpha*(nTrial/(1+exp(-(beta[1]+beta[2]*s+unorm[j]+vnorm[j]*s)))))
          return(h)}
        haz<-function(t){
          exp(-integrate(h,dat$x0[j],t)$value)-runif(1)}
        e.time<-c(e.time,uniroot(haz, c(0,300))$root)
      }
      
      c.times <- e.time
      censor <- sample(1:nRand,floor(cens*nRand)) #los individuos censurados si o si
      c.times[censor] <- runif(floor(cens*nRand),dat[censor,]$x0,pmin(e.time,dat[,ncol(dat)])[censor])  #El tiempo de censura de estos individuos
      
      obstimereal<-pmin(e.time,c.times,dat[,ncol(dat)])
      statusreal<-ifelse(obstimereal==e.time,1,0)
      
  
      resultado<-dplyr::bind_rows(lapply(1:nRand, function(i) dat.y[which(dat.y$z==i&dat.y$time<=obstimereal[i]),]))    #datos observados
      colnames(resultado)<-c("z","meas","x","y")
  
      
      dat.jm<-data.frame(as.factor(1:nRand),obstimereal,statusreal)
      
      colnames(dat.jm)<-c("id","obstime","status")
      
      
      #### JOINT MODELLING ASUMING NORMALITY ####
      
      fitlme<-try(nlme::lme(y~x, random=~x|z,data=resultado))
      coxFit <- try(survival::coxph(survival::Surv(dat.jm$obstime,dat.jm$status) ~1, data = dat.jm, cluster = id,x = TRUE))
      if (class(coxFit)[1]=="try-error"){
        no.conv <- 1
      }
      if (class(fitlme)[1]=="try-error"){
        no.conv <- 1
      }
      JM<-try(JM::jointModel(fitlme,coxFit,timeVar = "x", method = "piecewise-PH-aGH"))
      sum.JM<-try(summary(JM))
      if (class(JM)[1]=="try-error"){
        no.conv <- 1
        no.conv.JM<-no.conv.JM+1
        next
      } else if (JM$conv=="TRUE"){
        no.conv <- 1
        no.conv.JM<-no.conv.JM+1
        next
      } else if (class(sum.JM)[1]=="try-error"){
        no.conv <-1
        no.conv.JM<-no.conv.JM+1
        next
      }
      
      
          ############# TWO STAGE ###################
      
      #STEP 1
      
      rZ.i <- model.matrix(~z-1,data=resultado)
      rZ.t<-model.matrix(~z:x-1,data = resultado)
      rZ<-cbind(rZ.i,rZ.t)
      
      BB<- try(PROreg::BBmm(fixed.formula = y~x,Z=rZ,nRandComp =c(nRand,nRand),m=nTrial,data=resultado,maxiter = 100,show = TRUE))
      sum.BB<-try(summary(BB))
      if (class(BB)[1]=="try-error"){
        no.conv <- 1
        next
      } else if (BB$conv=="no"){
        no.conv <- 1
        next
      }else if (class(sum.BB)[1]=="try-error"){
        no.conv <- 1
        next
      }
      
      #STEP 2
        X.<-cbind(1,obstimereal)
        Z..i<-model.matrix(~-1+as.factor(1:nRand))
        Z..t<-model.matrix(~-1+as.factor(1:nRand):obstimereal)
        Z..<- cbind(Z..i,Z..t)
        etaY<-try(X.%*%BB$fixed.coef+Z..%*%BB$random.coef)
        exY<-nTrial*exp(etaY)/(1+exp(etaY))
       
        dat.jm$exY<-exY
        TS<-try(survival::coxph(survival::Surv(obstime,status)~exY, data = dat.jm))
        sum.TS<-summary(TS)
        

        
      ###LONGITUDINAL COX DATA ####
      wide <- tidyr::spread(resultado[,-4], key = meas, value = x)
      dat.lc<- survival::tmerge(data1=wide, data2=wide, id=z, death = event(obstimereal+0.0001,statusreal))
      for (i in 0:(max(count(resultado$z)$freq)-1)) {
        nam <- as.symbol(paste("x", i,sep = ""))
        dat.lc<- survival::tmerge(dat.lc, wide, id=z, measure = event(wide[[nam]]))
      }
      dat.lc<-dat.lc[-which(dat.lc$tstart==0),]
      dat.lc$y<-resultado$y
      
      # TIME DEPENDENT COX
      TVCM<-try(survival::coxph(survival::Surv(tstart,tstop,death)~y,data = dat.lc,id=z))
      sum.TVCM<-summary(TVCM)
      if (class(TVCM)[1]=="try-error"){
        no.conv <- 1
      }
      
    }
    
   
    
    # Generated times
    
    cens.all<-cbind(cens.all,cens) 
    obstime.all<-cbind(obstime.all,obstimereal)  
    status.all<-cbind(status.all,statusreal) 
    e.time.all<-cbind(e.time.all,e.time) 
    y.all <- cbind(y.all,y) 
    
    #recogemos valores
    coef.bb.all <- cbind(coef.bb.all,BB$fixed.coef)
    se.bb.all <- cbind(se.bb.all,sqrt(diag(BB$fixed.vcov))) 
    
    rand.bb.all<-cbind(rand.bb.all,BB$random.coef)
    sigmau.bb.all <- cbind(sigmau.bb.all,BB$sigma.coef[1])
    sigmav.bb.all <- cbind(sigmav.bb.all,BB$sigma.coef[2])
    se.bb.sigmau <- cbind(se.bb.sigmau,sum.BB$sigma[1,2])
    se.bb.sigmav <- cbind(se.bb.sigmav,sum.BB$sigma[2,2])
    
    psi.bb.all<-cbind(psi.bb.all, BB$psi.coef)
    se.psi<-cbind(se.psi,sum.BB$psi.table[1,2])
    
    coef.norm.all<-cbind(coef.norm.all,JM$coefficients$betas)
    se.norm.all<-cbind(se.norm.all, sum.JM$`CoefTable-Long`[,2])
    
    sigmau.norm.all<-c(sigmau.norm.all,sqrt(JM$coefficients$D[1,1]))
    sigmav.norm.all<-c(sigmav.norm.all,sqrt(JM$coefficients$D[2,2]))
    
    epsilon.norm.all<-c(epsilon.norm.all,JM$coefficients$sigma)
    
    JMcoeff.all<-cbind(JMcoeff.all,JM$coefficients$alpha)
    JMse.all<-cbind(JMse.all, sum.JM$`CoefTable-Event`[1,2])
    
    TScoeff.all<-c(TScoeff.all,TS$coefficients) 
    TSse.all<-c(TSse.all,sqrt(diag(TS$var))) 
    
    TVCMcoeff.all<-c(TVCMcoeff.all,TVCM$coefficients) 
    TVCMse.all<-c(TVCMse.all,sqrt(diag(TVCM$var))) 
    
    if (sum.TS$coefficients[,5]<0.05){
      significance.TS <- significance.TS+1
    }
    
    if (sum.JM$`CoefTable-Event`[1,4]<0.05){
      significance.JM <- significance.JM+1
    }
    
    if (sum.TVCM$coefficients[,5]<0.05){
      significance.TVCM <- significance.TVCM+1
    }
    n<-n+1
    print(paste0("Number of iteration ",n))
  }
  
  #### OUTPUTS
  #LONGITUDINAL RESULTS
  #fixed effects
  BB.mv.beta<- apply(coef.bb.all,1,mean)
  BB.V.beta <- apply(coef.bb.all,1,var)
  BB.sds.beta <- sqrt(BB.V.beta)
  BB.BIAS.beta0 <- mean(coef.bb.all[1,])-beta[1]
  BB.BIAS.beta1 <- mean(coef.bb.all[2,])-beta[2]
  BB.BIAS.beta<-c(BB.BIAS.beta0,BB.BIAS.beta1)
  BB.EMS.beta0 <- (BB.BIAS.beta0)^2+BB.V.beta[1]
  BB.EMS.beta1 <- (BB.BIAS.beta1)^2+BB.V.beta[2]
  BB.EMS.beta<-c(BB.EMS.beta0,BB.EMS.beta1)
  BB.ASE.beta0<-mean(se.bb.all[1,])
  BB.ASE.beta1<-mean(se.bb.all[2,])
  BB.ASE.beta<-c(BB.ASE.beta0,BB.ASE.beta1)
  BB.long.table<-cbind(BB.mv.beta,BB.sds.beta,BB.BIAS.beta,BB.EMS.beta,BB.ASE.beta)
  colnames(BB.long.table)<-c("est","sds","BIAS","EMS","ASE")
  rownames(BB.long.table)<-c("beta0","beta1")
  prop.bb <- proportion(coef.bb.all,se.bb.all,beta)
  prop.table <- cbind(t(prop.bb$in.ci))
  prop.both <- cbind(prop.bb$proportion)
  prop <- apply(prop.table,2,sum)/nSim
  JM.mv.beta<- apply(coef.norm.all,1,mean)
  JM.V.beta <- apply(coef.norm.all,1,var)
  JM.sds.beta <- sqrt(JM.V.beta)
  JM.ASE.beta0<-mean(se.norm.all[1,])
  JM.ASE.beta1<-mean(se.norm.all[2,])
  JM.ASE.beta<-c(JM.ASE.beta0,JM.ASE.beta1)
  JM.long.table<-cbind(JM.mv.beta,JM.sds.beta,JM.ASE.beta)
  colnames(JM.long.table)<-c("est","sds","ASE")
  rownames(JM.long.table)<-c("beta0","beta1")
  
  #random and variable effects
  BB.sigmau.est<-mean(sigmau.bb.all)
  BB.V.sigmau<- apply(sigmau.bb.all,1,var)
  BB.sds.sigmau<- sqrt(BB.V.sigmau)
  BB.BIAS.sigmau<- BB.sigmau.est-sigmau
  BB.EMS.sigmau<- BB.BIAS.sigmau^2+BB.V.sigmau
  BB.ASE.sigmau<-mean(se.bb.sigmau)
  
  BB.sigmav.est<-mean(sigmav.bb.all)
  BB.V.sigmav<- apply(sigmav.bb.all,1,var)
  BB.sds.sigmav<- sqrt(BB.V.sigmav)
  BB.BIAS.sigmav<- BB.sigmav.est-sigmav
  
  BB.EMS.sigmav<- BB.BIAS.sigmav^2+BB.V.sigmav
  BB.ASE.sigmav<-mean(se.bb.sigmav)
  
  
  BB.sigma<-cbind(rbind(BB.sigmau.est,BB.sigmav.est),rbind(BB.sds.sigmau,BB.sds.sigmav),rbind(BB.BIAS.sigmau,BB.BIAS.sigmav),rbind(BB.EMS.sigmau,BB.EMS.sigmav),rbind(BB.ASE.sigmau,BB.ASE.sigmav))
  
  psi.est<-mean(psi.bb.all)
  V.psi<-apply(psi.bb.all,1,var)
  sds.psi<-sqrt(V.psi)
  BIAS.psi<-psi.est-log(phi)
  EMS.psi<-BIAS.psi^2+V.psi
  ASE.psi<-mean(se.psi)
  psi<-cbind(psi.est,sds.psi,BIAS.psi,EMS.psi,ASE.psi)
  
  BB.var.table<-rbind(BB.sigma,psi)
  
  colnames(BB.var.table)<-c("est","SDS","BIAS","EMS","ASE")
  rownames(BB.var.table)<-c("sigmau","sigmav","psi")
  
  JM.sigmau.est<-mean(sigmau.norm.all)
  JM.V.sigmau<-var(sigmau.norm.all)
  JM.sds.sigmau<- sqrt(JM.V.sigmau)

  JM.sigmav.est<-mean(sigmav.norm.all)
  JM.V.sigmav<-var(sigmav.norm.all)
  JM.sds.sigmav<- sqrt(JM.V.sigmav)
  
  JM.sigma<-cbind(rbind(JM.sigmau.est,JM.sigmav.est),rbind(JM.sds.sigmau,JM.sds.sigmav))
  
  epsilon.est<-mean(epsilon.norm.all)
  epsilon.sds<-sqrt(var(epsilon.norm.all))
  epsilon<-cbind(epsilon.est,epsilon.sds)
  JM.var.table<-rbind(JM.sigma,epsilon)
  colnames(JM.var.table)<-c("est","SDS")
  rownames(JM.var.table)<-c("sigmau","sigmav","epsilon")
  
  #SURVIVAL RESULTS
  TS.mv.alpha<-mean(TScoeff.all)
  TVCM.mv.alpha<-mean(TVCMcoeff.all)
  JM.mv.alpha<-mean(JMcoeff.all)
  
  V.TS<- var(TScoeff.all)
  V.TVCM<- var(TVCMcoeff.all)
  V.JM<- apply(JMcoeff.all,1,var)
  
  BIAS.TS <- TS.mv.alpha-alpha
  BIAS.TVCM <- TVCM.mv.alpha-alpha
  BIAS.JM<- JM.mv.alpha-alpha
  
  EMS.TS<- (BIAS.TS)^2+V.TS
  EMS.TVCM <- (BIAS.TVCM)^2+V.TVCM
  EMS.JM <- (BIAS.JM)^2+V.JM
  
  means <- cbind(TS.mv.alpha,TVCM.mv.alpha,JM.mv.alpha)
  sds <- cbind(sqrt(V.TS),sqrt((V.TVCM)),sqrt(V.JM))
  BIAS<-cbind(BIAS.TS,BIAS.TVCM,BIAS.JM)
  EMS<-cbind(EMS.TS,EMS.TVCM,EMS.JM)
  mean.sd <- rbind(means,sds,BIAS,EMS)
  rownames(mean.sd) <- c("alpha","sd","BIAS","EMS")
  colnames(mean.sd)<-c("TS","TVCM","JM")
  prop.TS.alpha <- proportion(TScoeff.all,TSse.all,alpha)
  prop.TVCM.alpha<- proportion(TVCMcoeff.all,TVCMse.all,alpha)
  prop.JM.alpha<- proportion(JMcoeff.all,JMse.all,alpha)
  
  prop.alpha.table <- cbind(prop.TS.alpha$in.ci,prop.TVCM.alpha$in.ci,c(prop.JM.alpha$in.ci))
  colnames(prop.alpha.table)<-c("TS","TVCM","JM")
  
  prop.alpha.both <- cbind(prop.TS.alpha$proportion,prop.TVCM.alpha$proportion,prop.JM.alpha$proportion)
  
  prop.alpha <- apply(prop.alpha.table,2,sum)/nSim
  
  out <- list(Survival.effect.table=t(mean.sd),BB.Longitudinal.fixed=BB.long.table,BB.Longitudinal.var=BB.var.table,prop.in.ic.beta=prop,JM.Longitudinal.fixed=JM.long.table,JM.Longitudinal.var=JM.var.table,prop.in.ic.alpha=prop.alpha,y.all=y.all,e.time.all=e.time.all,obstime.all=obstime.all,status.all=status.all,coef.bb.all=coef.bb.all,rand.bb.all=rand.bb.all,significance.TS=significance.TS,significance.TVCM=significance.TVCM,significance.JM=significance.JM,no.conv.JM=no.conv.JM ,TS.estimation= TScoeff.all,TVCM.estimation=TVCMcoeff.all,JM.estimation=c(JMcoeff.all))
  out
}
