setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")
library(spectrolab)
library(patchwork)
library(pls)
library(ggplot2)
library(reshape2)

#########################################
## this is just a demonstration of using
## brightness-normalized reflectance in PLSR

#########################################
## read data

intact_spec_agg_train<-readRDS("SavedResults/intact_spec_agg_train.rds")
intact_spec_agg_test<-readRDS("SavedResults/intact_spec_agg_test.rds")

intact_spec_agg_train_bn<-t(apply(as.matrix(intact_spec_agg_train),1,function(x) x/sqrt(sum(x^2))))
intact_spec_agg_train_bn<-spectra(intact_spec_agg_train_bn,bands=bands(intact_spec_agg_train),names=names(intact_spec_agg_train))
meta(intact_spec_agg_train_bn)<-meta(intact_spec_agg_train)

intact_spec_agg_test_bn<-t(apply(as.matrix(intact_spec_agg_test),1,function(x) x/sqrt(sum(x^2))))
intact_spec_agg_test_bn<-spectra(intact_spec_agg_test_bn,bands=bands(intact_spec_agg_test),names=names(intact_spec_agg_test))
meta(intact_spec_agg_test_bn)<-meta(intact_spec_agg_test)

source("Scripts/senesced-trait-models/useful_functions.R")

#########################################
## train models

sol_intact<-plsr(meta(intact_spec_agg_train_bn)$solubles~as.matrix(intact_spec_agg_train_bn),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
hemi_intact<-plsr(meta(intact_spec_agg_train_bn)$hemicellulose~as.matrix(intact_spec_agg_train_bn),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
recalc_intact<-plsr(meta(intact_spec_agg_train_bn)$recalcitrant~as.matrix(intact_spec_agg_train_bn),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
Nmass_intact<-plsr(meta(intact_spec_agg_train_bn)$Nmass~as.matrix(intact_spec_agg_train_bn),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
Cmass_intact<-plsr(meta(intact_spec_agg_train_bn)$Cmass~as.matrix(intact_spec_agg_train_bn),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
LMA_intact<-plsr(meta(intact_spec_agg_train_bn)$LMA~as.matrix(intact_spec_agg_train_bn),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
Narea_intact<-plsr(meta(intact_spec_agg_train_bn)$Narea~as.matrix(intact_spec_agg_train_bn),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
Carea_intact<-plsr(meta(intact_spec_agg_train_bn)$Carea~as.matrix(intact_spec_agg_train_bn),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_sol_intact <- selectNcomp(sol_intact, method = "onesigma", plot = FALSE)
sol_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$solubles))
sol_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[sol_intact_valid],
                            Species=meta(intact_spec_agg_train_bn)$sp[sol_intact_valid],
                            Run=meta(intact_spec_agg_train_bn)$FiberRun[sol_intact_valid],
                            Measured=meta(intact_spec_agg_train_bn)$solubles[sol_intact_valid],
                            val_pred=sol_intact$validation$pred[,,ncomp_sol_intact])
ggplot(sol_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % solubles from intact-leaf spectra")+guides(color=F)

ncomp_hemi_intact <- selectNcomp(hemi_intact, method = "onesigma", plot = FALSE)
hemi_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$hemicellulose))
hemi_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[hemi_intact_valid],
                             Species=meta(intact_spec_agg_train_bn)$sp[hemi_intact_valid],
                             Run=meta(intact_spec_agg_train_bn)$FiberRun[hemi_intact_valid],
                             Measured=meta(intact_spec_agg_train_bn)$hemicellulose[hemi_intact_valid],
                             val_pred=hemi_intact$validation$pred[,,ncomp_hemi_intact])
ggplot(hemi_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(5,27),ylim=c(5,27))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % hemicellulose from intact-leaf spectra")+guides(color=F)

ncomp_recalc_intact <- selectNcomp(recalc_intact, method = "onesigma", plot = FALSE)
recalc_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$recalcitrant))
recalc_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[recalc_intact_valid],
                               Species=meta(intact_spec_agg_train_bn)$sp[recalc_intact_valid],
                               Run=meta(intact_spec_agg_train_bn)$FiberRun[recalc_intact_valid],
                               Measured=meta(intact_spec_agg_train_bn)$recalcitrant[recalc_intact_valid],
                               val_pred=recalc_intact$validation$pred[,,ncomp_recalc_intact])
ggplot(recalc_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % recalcitrants from intact-leaf spectra")+guides(color=F)

ncomp_Cmass_intact <- selectNcomp(Cmass_intact, method = "onesigma", plot = FALSE)
Cmass_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$Cmass))
Cmass_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[Cmass_intact_valid],
                              Species=meta(intact_spec_agg_train_bn)$sp[Cmass_intact_valid],
                              Run=meta(intact_spec_agg_train_bn)$FiberRun[Cmass_intact_valid],
                              Measured=meta(intact_spec_agg_train_bn)$Cmass[Cmass_intact_valid],
                              val_pred=Cmass_intact$validation$pred[,,ncomp_Cmass_intact])
ggplot(Cmass_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %C from intact-leaf spectra")

ncomp_Nmass_intact <- selectNcomp(Nmass_intact, method = "onesigma", plot = FALSE)
Nmass_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$Nmass))
Nmass_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[Nmass_intact_valid],
                              Species=meta(intact_spec_agg_train_bn)$sp[Nmass_intact_valid],
                              Run=meta(intact_spec_agg_train_bn)$FiberRun[Nmass_intact_valid],
                              Measured=meta(intact_spec_agg_train_bn)$Nmass[Nmass_intact_valid],
                              val_pred=Nmass_intact$validation$pred[,,ncomp_Nmass_intact])
ggplot(Nmass_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %N from intact-leaf spectra")

ncomp_LMA_intact <- selectNcomp(LMA_intact, method = "onesigma", plot = FALSE)
LMA_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$LMA))
LMA_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[LMA_intact_valid],
                            Species=meta(intact_spec_agg_train_bn)$sp[LMA_intact_valid],
                            Run=meta(intact_spec_agg_train_bn)$FiberRun[LMA_intact_valid],
                            Measured=meta(intact_spec_agg_train_bn)$LMA[LMA_intact_valid],
                            val_pred=LMA_intact$validation$pred[,,ncomp_LMA_intact])
ggplot(LMA_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting LMA from intact-leaf spectra")+guides(color=F)

ncomp_Carea_intact <- selectNcomp(Carea_intact, method = "onesigma", plot = FALSE)
Carea_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$Carea))
Carea_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[Carea_intact_valid],
                              Species=meta(intact_spec_agg_train_bn)$sp[Carea_intact_valid],
                              Run=meta(intact_spec_agg_train_bn)$FiberRun[Carea_intact_valid],
                              Measured=meta(intact_spec_agg_train_bn)$Carea[Carea_intact_valid],
                              val_pred=Carea_intact$validation$pred[,,ncomp_Carea_intact])
ggplot(Carea_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,200),ylim=c(0,200))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting C per area from intact-leaf spectra")

ncomp_Narea_intact <- selectNcomp(Narea_intact, method = "onesigma", plot = FALSE)
Narea_intact_valid <- which(!is.na(meta(intact_spec_agg_train_bn)$Narea))
Narea_intact_pred<-data.frame(ID=meta(intact_spec_agg_train_bn)$ID[Narea_intact_valid],
                              Species=meta(intact_spec_agg_train_bn)$sp[Narea_intact_valid],
                              Run=meta(intact_spec_agg_train_bn)$FiberRun[Narea_intact_valid],
                              Measured=meta(intact_spec_agg_train_bn)$Narea[Narea_intact_valid],
                              val_pred=Narea_intact$validation$pred[,,ncomp_Narea_intact])
ggplot(Narea_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured (g N/cm2)",y="Predicted (g N/cm2)")+
  ggtitle("Predicting N per area from intact-leaf spectra")

#######################################
## jackknife tests + prediction of validation data

sol_jack_coefs<-list()
hemi_jack_coefs<-list()
recalc_jack_coefs<-list()
Cmass_jack_coefs<-list()
Nmass_jack_coefs<-list()
Carea_jack_coefs<-list()
Narea_jack_coefs<-list()
LMA_jack_coefs<-list()

sol_jack_stats<-list()
hemi_jack_stats<-list()
recalc_jack_stats<-list()
Cmass_jack_stats<-list()
Nmass_jack_stats<-list()
Carea_jack_stats<-list()
Narea_jack_stats<-list()
LMA_jack_stats<-list()
nreps<-200

for(i in 1:nreps){
  print(i)
  
  n_cal_spec<-nrow(intact_spec_agg_train_bn)
  train_jack<-sample(1:n_cal_spec,floor(0.7*n_cal_spec))
  test_jack<-setdiff(1:n_cal_spec,train_jack)
  
  calib_jack<-intact_spec_agg_train_bn[train_jack]
  val_jack<-intact_spec_agg_train_bn[test_jack]
  
  sol_intact_jack<-plsr(meta(calib_jack)$solubles~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  hemi_intact_jack<-plsr(meta(calib_jack)$hemicellulose~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  recalc_intact_jack<-plsr(meta(calib_jack)$recalcitrant~as.matrix(calib_jack),
                           ncomp=30,method = "oscorespls",validation="none")
  Cmass_intact_jack<-plsr(meta(calib_jack)$Cmass~as.matrix(calib_jack),
                          ncomp=30,method = "oscorespls",validation="none")
  Nmass_intact_jack<-plsr(meta(calib_jack)$Nmass~as.matrix(calib_jack),
                          ncomp=30,method = "oscorespls",validation="none")
  Carea_intact_jack<-plsr(meta(calib_jack)$Carea~as.matrix(calib_jack),
                          ncomp=30,method = "oscorespls",validation="none")
  Narea_intact_jack<-plsr(meta(calib_jack)$Narea~as.matrix(calib_jack),
                          ncomp=30,method = "oscorespls",validation="none")
  LMA_intact_jack<-plsr(meta(calib_jack)$LMA~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  
  sol_jack_val_pred<-as.vector(predict(sol_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_sol_intact)[,,1])
  sol_jack_val_fit<-lm(sol_jack_val_pred~meta(val_jack)$solubles)
  sol_jack_stats[[i]]<-c(R2=summary(sol_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$solubles,sol_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$solubles,sol_jack_val_pred,0.025,0.975),
                         bias=mean(sol_jack_val_pred,na.rm=T)-mean(meta(val_jack)$solubles,na.rm=T))
  
  hemi_jack_val_pred<-as.vector(predict(hemi_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_hemi_intact)[,,1])
  hemi_jack_val_fit<-lm(hemi_jack_val_pred~meta(val_jack)$hemicellulose)
  hemi_jack_stats[[i]]<-c(R2=summary(hemi_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$hemicellulose,hemi_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$hemicellulose,hemi_jack_val_pred,0.025,0.975),
                          bias=mean(hemi_jack_val_pred,na.rm=T)-mean(meta(val_jack)$hemicellulose,na.rm=T))
  
  recalc_jack_val_pred<-as.vector(predict(recalc_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_recalc_intact)[,,1])
  recalc_jack_val_fit<-lm(recalc_jack_val_pred~meta(val_jack)$recalcitrant)
  recalc_jack_stats[[i]]<-c(R2=summary(recalc_jack_val_fit)$r.squared,
                            RMSE=RMSD(meta(val_jack)$recalcitrant,recalc_jack_val_pred),
                            perRMSE=percentRMSD(meta(val_jack)$recalcitrant,recalc_jack_val_pred,0.025,0.975),
                            bias=mean(recalc_jack_val_pred,na.rm=T)-mean(meta(val_jack)$recalcitrant,na.rm=T))
  
  Cmass_jack_val_pred<-as.vector(predict(Cmass_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Cmass_intact)[,,1])
  Cmass_jack_val_fit<-lm(Cmass_jack_val_pred~meta(val_jack)$Cmass)
  Cmass_jack_stats[[i]]<-c(R2=summary(Cmass_jack_val_fit)$r.squared,
                           RMSE=RMSD(meta(val_jack)$Cmass,Cmass_jack_val_pred),
                           perRMSE=percentRMSD(meta(val_jack)$Cmass,Cmass_jack_val_pred,0.025,0.975),
                           bias=mean(Cmass_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Cmass,na.rm=T))
  
  Nmass_jack_val_pred<-as.vector(predict(Nmass_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Nmass_intact)[,,1])
  Nmass_jack_val_fit<-lm(Nmass_jack_val_pred~meta(val_jack)$Nmass)
  Nmass_jack_stats[[i]]<-c(R2=summary(Nmass_jack_val_fit)$r.squared,
                           RMSE=RMSD(meta(val_jack)$Nmass,Nmass_jack_val_pred),
                           perRMSE=percentRMSD(meta(val_jack)$Nmass,Nmass_jack_val_pred,0.025,0.975),
                           bias=mean(Nmass_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Nmass,na.rm=T))
  
  Carea_jack_val_pred<-as.vector(predict(Carea_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Carea_intact)[,,1])
  Carea_jack_val_fit<-lm(Carea_jack_val_pred~meta(val_jack)$Carea)
  Carea_jack_stats[[i]]<-c(R2=summary(Carea_jack_val_fit)$r.squared,
                           RMSE=RMSD(meta(val_jack)$Carea,Carea_jack_val_pred),
                           perRMSE=percentRMSD(meta(val_jack)$Carea,Carea_jack_val_pred,0.025,0.975),
                           bias=mean(Carea_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Carea,na.rm=T))
  
  Narea_jack_val_pred<-as.vector(predict(Narea_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Narea_intact)[,,1])
  Narea_jack_val_fit<-lm(Narea_jack_val_pred~meta(val_jack)$Narea)
  Narea_jack_stats[[i]]<-c(R2=summary(Narea_jack_val_fit)$r.squared,
                           RMSE=RMSD(meta(val_jack)$Narea,Narea_jack_val_pred),
                           perRMSE=percentRMSD(meta(val_jack)$Narea,Narea_jack_val_pred,0.025,0.975),
                           bias=mean(Narea_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Narea,na.rm=T))
  
  LMA_jack_val_pred<-as.vector(predict(LMA_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_LMA_intact)[,,1])
  LMA_jack_val_fit<-lm(LMA_jack_val_pred~meta(val_jack)$LMA)
  LMA_jack_stats[[i]]<-c(R2=summary(LMA_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$LMA,LMA_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$LMA,LMA_jack_val_pred,0.025,0.975),
                         bias=mean(LMA_jack_val_pred,na.rm=T)-mean(meta(val_jack)$LMA,na.rm=T))
  
  sol_jack_coefs[[i]]<-as.vector(coef(sol_intact_jack,ncomp=ncomp_sol_intact,intercept=TRUE))
  hemi_jack_coefs[[i]]<-as.vector(coef(hemi_intact_jack,ncomp=ncomp_hemi_intact,intercept=TRUE))
  recalc_jack_coefs[[i]]<-as.vector(coef(recalc_intact_jack,ncomp=ncomp_recalc_intact,intercept=TRUE))
  Cmass_jack_coefs[[i]]<-as.vector(coef(Cmass_intact_jack,ncomp=ncomp_Cmass_intact,intercept=TRUE))
  Nmass_jack_coefs[[i]]<-as.vector(coef(Nmass_intact_jack,ncomp=ncomp_Nmass_intact,intercept=TRUE))
  Carea_jack_coefs[[i]]<-as.vector(coef(Carea_intact_jack,ncomp=ncomp_Carea_intact,intercept=TRUE))
  Narea_jack_coefs[[i]]<-as.vector(coef(Narea_intact_jack,ncomp=ncomp_Narea_intact,intercept=TRUE))
  LMA_jack_coefs[[i]]<-as.vector(coef(LMA_intact_jack,ncomp=ncomp_LMA_intact,intercept=TRUE))
  
}

sol_jack_pred<-apply.coefs(sol_jack_coefs,as.matrix(intact_spec_agg_test_bn))
sol_jack_stat<-t(apply(sol_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
sol_jack_df<-data.frame(pred_mean=sol_jack_stat[,1],
                        pred_low=sol_jack_stat[,2],
                        pred_high=sol_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test_bn)$solubles,
                        ncomp=ncomp_sol_intact,
                        Species=meta(intact_spec_agg_test_bn)$sp,
                        ID=meta(intact_spec_agg_test_bn)$ID)

hemi_jack_pred<-apply.coefs(hemi_jack_coefs,as.matrix(intact_spec_agg_test_bn))
hemi_jack_stat<-t(apply(hemi_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemi_jack_df<-data.frame(pred_mean=hemi_jack_stat[,1],
                         pred_low=hemi_jack_stat[,2],
                         pred_high=hemi_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test_bn)$hemicellulose,
                         ncomp=ncomp_hemi_intact,
                         Species=meta(intact_spec_agg_test_bn)$sp,
                         ID=meta(intact_spec_agg_test_bn)$ID)

recalc_jack_pred<-apply.coefs(recalc_jack_coefs,as.matrix(intact_spec_agg_test_bn))
recalc_jack_stat<-t(apply(recalc_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
recalc_jack_df<-data.frame(pred_mean=recalc_jack_stat[,1],
                           pred_low=recalc_jack_stat[,2],
                           pred_high=recalc_jack_stat[,3],
                           Measured=meta(intact_spec_agg_test_bn)$recalcitrant,
                           ncomp=ncomp_recalc_intact,
                           Species=meta(intact_spec_agg_test_bn)$sp,
                           ID=meta(intact_spec_agg_test_bn)$ID)

Cmass_jack_pred<-apply.coefs(Cmass_jack_coefs,as.matrix(intact_spec_agg_test_bn))
Cmass_jack_stat<-t(apply(Cmass_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass_jack_df<-data.frame(pred_mean=Cmass_jack_stat[,1],
                          pred_low=Cmass_jack_stat[,2],
                          pred_high=Cmass_jack_stat[,3],
                          Measured=meta(intact_spec_agg_test_bn)$Cmass,
                          ncomp=ncomp_Cmass_intact,
                          Species=meta(intact_spec_agg_test_bn)$sp,
                          ID=meta(intact_spec_agg_test_bn)$ID)

Nmass_jack_pred<-apply.coefs(Nmass_jack_coefs,as.matrix(intact_spec_agg_test_bn))
Nmass_jack_stat<-t(apply(Nmass_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass_jack_df<-data.frame(pred_mean=Nmass_jack_stat[,1],
                          pred_low=Nmass_jack_stat[,2],
                          pred_high=Nmass_jack_stat[,3],
                          Measured=meta(intact_spec_agg_test_bn)$Nmass,
                          ncomp=ncomp_Nmass_intact,
                          Species=meta(intact_spec_agg_test_bn)$sp,
                          ID=meta(intact_spec_agg_test_bn)$ID)

Carea_jack_pred<-apply.coefs(Carea_jack_coefs,as.matrix(intact_spec_agg_test_bn))
Carea_jack_stat<-t(apply(Carea_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Carea_jack_df<-data.frame(pred_mean=Carea_jack_stat[,1],
                          pred_low=Carea_jack_stat[,2],
                          pred_high=Carea_jack_stat[,3],
                          Measured=meta(intact_spec_agg_test_bn)$Carea,
                          ncomp=ncomp_Carea_intact,
                          Species=meta(intact_spec_agg_test_bn)$sp,
                          ID=meta(intact_spec_agg_test_bn)$ID)

Narea_jack_pred<-apply.coefs(Narea_jack_coefs,as.matrix(intact_spec_agg_test_bn))
Narea_jack_stat<-t(apply(Narea_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Narea_jack_df<-data.frame(pred_mean=Narea_jack_stat[,1],
                          pred_low=Narea_jack_stat[,2],
                          pred_high=Narea_jack_stat[,3],
                          Measured=meta(intact_spec_agg_test_bn)$Narea,
                          ncomp=ncomp_Narea_intact,
                          Species=meta(intact_spec_agg_test_bn)$sp,
                          ID=meta(intact_spec_agg_test_bn)$ID)

LMA_jack_pred<-apply.coefs(LMA_jack_coefs,as.matrix(intact_spec_agg_test_bn))
LMA_jack_stat<-t(apply(LMA_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df<-data.frame(pred_mean=LMA_jack_stat[,1],
                        pred_low=LMA_jack_stat[,2],
                        pred_high=LMA_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test_bn)$LMA,
                        ncomp=ncomp_LMA_intact,
                        Species=meta(intact_spec_agg_test_bn)$sp,
                        ID=meta(intact_spec_agg_test_bn)$ID)

intact_jack_df_list_bn<-list(LMA=LMA_jack_df,
                          Cmass=Cmass_jack_df,
                          Nmass=Nmass_jack_df,
                          sol=sol_jack_df,
                          hemi=hemi_jack_df,
                          recalc=recalc_jack_df,
                          Carea=Carea_jack_df,
                          Narea=Narea_jack_df)
summ<-data.frame(ncomp=unlist(lapply(intact_jack_df_list_bn,function(x) x$ncomp[1])),
                 r2=round(unlist(lapply(intact_jack_df_list_bn,function(x) summary(lm(Measured~pred_mean,data=x))$r.squared)),3),
                 rmse=signif(unlist(lapply(intact_jack_df_list_bn,function(x) RMSD(x$Measured,x$pred_mean))),3),
                 perrmse=signif(unlist(lapply(intact_jack_df_list_bn,function(x) percentRMSD(x$Measured,x$pred_mean,0.025,0.975)))*100,3))
