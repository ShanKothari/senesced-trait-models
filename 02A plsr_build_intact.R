## PLSR of fiber data
setwd("C:/Users/kotha020/Dropbox/TraitModels2018/SenescencePaper/")
library(spectrolab)
library(pls)
library(ggplot2)
library(vegan)
library(rdist)
library(reshape2)

#########################################
## read data

intact_spec_agg_train<-readRDS("SavedResults/intact_spec_agg_train.rds")
intact_spec_agg_test<-readRDS("SavedResults/intact_spec_agg_test.rds")

source("Scripts/senesced-trait-models/useful_functions.R")

#########################################
## train models

NDF_intact<-plsr(meta(intact_spec_agg_train)$NDF~as.matrix(intact_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ADF_intact<-plsr(meta(intact_spec_agg_train)$ADF~as.matrix(intact_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
perN_intact<-plsr(meta(intact_spec_agg_train)$perN~as.matrix(intact_spec_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
perC_intact<-plsr(meta(intact_spec_agg_train)$perC~as.matrix(intact_spec_agg_train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
LMA_intact<-plsr(meta(intact_spec_agg_train)$LMA~as.matrix(intact_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
perN_area_intact<-plsr(meta(intact_spec_agg_train)$perN_area~as.matrix(intact_spec_agg_train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
perC_area_intact<-plsr(meta(intact_spec_agg_train)$perC_area~as.matrix(intact_spec_agg_train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_NDF_intact <- selectNcomp(NDF_intact, method = "onesigma", plot = FALSE)
NDF_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$NDF))
NDF_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[NDF_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[NDF_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[NDF_intact_valid],
                            Measured=meta(intact_spec_agg_train)$NDF[NDF_intact_valid],
                            val_pred=NDF_intact$validation$pred[,,ncomp_NDF_intact])
ggplot(NDF_intact_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %NDF from intact-leaf spectra")+guides(color=F)

ncomp_ADF_intact <- selectNcomp(ADF_intact, method = "onesigma", plot = FALSE)
ADF_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$ADF))
ADF_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[ADF_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[ADF_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[ADF_intact_valid],
                            Measured=meta(intact_spec_agg_train)$ADF[ADF_intact_valid],
                            val_pred=ADF_intact$validation$pred[,,ncomp_ADF_intact])
ggplot(ADF_intact_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %ADF from intact-leaf spectra")+guides(color=F)

ncomp_perC_intact <- selectNcomp(perC_intact, method = "onesigma", plot = FALSE)
perC_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$perC))
perC_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[perC_intact_valid],
                             Species=meta(intact_spec_agg_train)$sp[perC_intact_valid],
                             Run=meta(intact_spec_agg_train)$FiberRun[perC_intact_valid],
                             Measured=meta(intact_spec_agg_train)$perC[perC_intact_valid],
                             val_pred=perC_intact$validation$pred[,,ncomp_perC_intact])
ggplot(perC_intact_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %C from intact-leaf spectra")

ncomp_perN_intact <- selectNcomp(perN_intact, method = "onesigma", plot = FALSE)
perN_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$perN))
perN_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[perN_intact_valid],
                             Species=meta(intact_spec_agg_train)$sp[perN_intact_valid],
                             Run=meta(intact_spec_agg_train)$FiberRun[perN_intact_valid],
                             Measured=meta(intact_spec_agg_train)$perN[perN_intact_valid],
                             val_pred=perN_intact$validation$pred[,,ncomp_perN_intact])
ggplot(perN_intact_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %N from intact-leaf spectra")

ncomp_LMA_intact <- selectNcomp(LMA_intact, method = "onesigma", plot = FALSE)
LMA_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$LMA))
LMA_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[LMA_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[LMA_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[LMA_intact_valid],
                            Measured=meta(intact_spec_agg_train)$LMA[LMA_intact_valid],
                            val_pred=LMA_intact$validation$pred[,,ncomp_LMA_intact])
ggplot(LMA_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.03),ylim=c(0,0.03))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting LMA from intact-leaf spectra")+guides(color=F)

ncomp_perC_area_intact <- selectNcomp(perC_area_intact, method = "onesigma", plot = FALSE)
perC_area_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$perC_area))
perC_area_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[perC_area_intact_valid],
                                  Species=meta(intact_spec_agg_train)$sp[perC_area_intact_valid],
                                  Run=meta(intact_spec_agg_train)$FiberRun[perC_area_intact_valid],
                                  Measured=meta(intact_spec_agg_train)$perC_area[perC_area_intact_valid],
                                  val_pred=perC_area_intact$validation$pred[,,ncomp_perC_area_intact])
ggplot(perC_area_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting C per area from intact-leaf spectra")

ncomp_perN_area_intact <- selectNcomp(perN_area_intact, method = "onesigma", plot = FALSE)
perN_area_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$perN_area))
perN_area_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[perN_area_intact_valid],
                                  Species=meta(intact_spec_agg_train)$sp[perN_area_intact_valid],
                                  Run=meta(intact_spec_agg_train)$FiberRun[perN_area_intact_valid],
                                  Measured=meta(intact_spec_agg_train)$perN_area[perN_area_intact_valid],
                                  val_pred=perN_area_intact$validation$pred[,,ncomp_perN_area_intact])
ggplot(perN_area_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.00025),ylim=c(0,0.00025))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured (g N/cm2)",y="Predicted (g N/cm2)")+
  ggtitle("Predicting N per area from intact-leaf spectra")

##########################################
## VIP plots

source("SenescencePaper/VIP.R")
VIP_intact<-data.frame(NDF=VIP(NDF_intact)[ncomp_NDF_intact,],
                       ADF=VIP(ADF_intact)[ncomp_ADF_intact,],
                       perC=VIP(perC_intact)[ncomp_perC_intact,],
                       perN=VIP(perN_intact)[ncomp_perN_intact,],
                       LMA=VIP(LMA_intact)[ncomp_LMA_intact,],
                       wavelength=400:2400)

#######################################
## jackknife tests + prediction of validation data

NDF_jack_coefs<-list()
ADF_jack_coefs<-list()
perC_jack_coefs<-list()
perN_jack_coefs<-list()
perC_area_jack_coefs<-list()
perN_area_jack_coefs<-list()
LMA_jack_coefs<-list()

NDF_jack_stats<-list()
ADF_jack_stats<-list()
perC_jack_stats<-list()
perN_jack_stats<-list()
perC_area_jack_stats<-list()
perN_area_jack_stats<-list()
LMA_jack_stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec<-nrow(intact_spec_agg_train)
  train_jack<-sample(1:n_cal_spec,floor(0.7*n_cal_spec))
  test_jack<-setdiff(1:n_cal_spec,train_jack)
  
  calib_jack<-intact_spec_agg_train[train_jack]
  val_jack<-intact_spec_agg_train[test_jack]
  
  NDF_intact_jack<-plsr(meta(calib_jack)$NDF~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  ADF_intact_jack<-plsr(meta(calib_jack)$ADF~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  perC_intact_jack<-plsr(meta(calib_jack)$perC~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  perN_intact_jack<-plsr(meta(calib_jack)$perN~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  perC_area_intact_jack<-plsr(meta(calib_jack)$perC_area~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  perN_area_intact_jack<-plsr(meta(calib_jack)$perN_area~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  LMA_intact_jack<-plsr(meta(calib_jack)$LMA~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  
  NDF_jack_val_pred<-as.vector(predict(NDF_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_NDF_intact)[,,1])
  NDF_jack_val_fit<-lm(NDF_jack_val_pred~meta(val_jack)$NDF)
  NDF_jack_stats[[i]]<-c(R2=summary(NDF_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$NDF,NDF_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$NDF,NDF_jack_val_pred,0.025,0.975),
                         bias=mean(NDF_jack_val_pred,na.rm=T)-mean(meta(val_jack)$NDF,na.rm=T))
  
  ADF_jack_val_pred<-as.vector(predict(ADF_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_ADF_intact)[,,1])
  ADF_jack_val_fit<-lm(ADF_jack_val_pred~meta(val_jack)$ADF)
  ADF_jack_stats[[i]]<-c(R2=summary(ADF_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$ADF,ADF_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$ADF,ADF_jack_val_pred,0.025,0.975),
                         bias=mean(ADF_jack_val_pred,na.rm=T)-mean(meta(val_jack)$ADF,na.rm=T))
  
  perC_jack_val_pred<-as.vector(predict(perC_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perC_intact)[,,1])
  perC_jack_val_fit<-lm(perC_jack_val_pred~meta(val_jack)$perC)
  perC_jack_stats[[i]]<-c(R2=summary(perC_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$perC,perC_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$perC,perC_jack_val_pred,0.025,0.975),
                          bias=mean(perC_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perC,na.rm=T))
  
  perN_jack_val_pred<-as.vector(predict(perN_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perN_intact)[,,1])
  perN_jack_val_fit<-lm(perN_jack_val_pred~meta(val_jack)$perN)
  perN_jack_stats[[i]]<-c(R2=summary(perN_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$perN,perN_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$perN,perN_jack_val_pred,0.025,0.975),
                          bias=mean(perN_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perN,na.rm=T))
  
  perC_area_jack_val_pred<-as.vector(predict(perC_area_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perC_area_intact)[,,1])
  perC_area_jack_val_fit<-lm(perC_area_jack_val_pred~meta(val_jack)$perC_area)
  perC_area_jack_stats[[i]]<-c(R2=summary(perC_area_jack_val_fit)$r.squared,
                               RMSE=RMSD(meta(val_jack)$perC_area,perC_area_jack_val_pred),
                               perRMSE=percentRMSD(meta(val_jack)$perC_area,perC_area_jack_val_pred,0.025,0.975),
                          bias=mean(perC_area_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perC_area,na.rm=T))
  
  perN_area_jack_val_pred<-as.vector(predict(perN_area_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perN_area_intact)[,,1])
  perN_area_jack_val_fit<-lm(perN_area_jack_val_pred~meta(val_jack)$perN_area)
  perN_area_jack_stats[[i]]<-c(R2=summary(perN_area_jack_val_fit)$r.squared,
                               RMSE=RMSD(meta(val_jack)$perN_area,perN_area_jack_val_pred),
                               perRMSE=percentRMSD(meta(val_jack)$perN_area,perN_area_jack_val_pred,0.025,0.975),
                          bias=mean(perN_area_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perN_area,na.rm=T))
  
  LMA_jack_val_pred<-as.vector(predict(LMA_intact_jack,newdata=as.matrix(val_jack),ncomp=ncomp_LMA_intact)[,,1])
  LMA_jack_val_fit<-lm(LMA_jack_val_pred~meta(val_jack)$LMA)
  LMA_jack_stats[[i]]<-c(R2=summary(LMA_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$LMA,LMA_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$LMA,LMA_jack_val_pred,0.025,0.975),
                         bias=mean(LMA_jack_val_pred,na.rm=T)-mean(meta(val_jack)$LMA,na.rm=T))
  
  NDF_jack_coefs[[i]]<-as.vector(coef(NDF_intact_jack,ncomp=ncomp_NDF_intact,intercept=TRUE))
  ADF_jack_coefs[[i]]<-as.vector(coef(ADF_intact_jack,ncomp=ncomp_ADF_intact,intercept=TRUE))
  perC_jack_coefs[[i]]<-as.vector(coef(perC_intact_jack,ncomp=ncomp_perC_intact,intercept=TRUE))
  perN_jack_coefs[[i]]<-as.vector(coef(perN_intact_jack,ncomp=ncomp_perN_intact,intercept=TRUE))
  perC_area_jack_coefs[[i]]<-as.vector(coef(perC_area_intact_jack,ncomp=ncomp_perC_area_intact,intercept=TRUE))
  perN_area_jack_coefs[[i]]<-as.vector(coef(perN_area_intact_jack,ncomp=ncomp_perN_area_intact,intercept=TRUE))
  LMA_jack_coefs[[i]]<-as.vector(coef(LMA_intact_jack,ncomp=ncomp_LMA_intact,intercept=TRUE))
  
}

NDF_jack_pred<-apply.coefs(NDF_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
NDF_jack_stat<-t(apply(NDF_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
NDF_jack_df<-data.frame(pred_mean=NDF_jack_stat[,1],
                        pred_low=NDF_jack_stat[,2],
                        pred_high=NDF_jack_stat[3],
                        Measured=meta(intact_spec_agg_test)$NDF,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

NDF_intact_val_plot<-ggplot(NDF_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured (%)",x="Predicted (%)")+
  ggtitle("Predicting NDF from intact-leaf spectra")+guides(color=F)

ADF_jack_pred<-apply.coefs(ADF_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
ADF_jack_stat<-t(apply(ADF_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
ADF_jack_df<-data.frame(pred_mean=ADF_jack_stat[,1],
                        pred_low=ADF_jack_stat[,2],
                        pred_high=ADF_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$ADF,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

ADF_intact_val_plot<-ggplot(ADF_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured (%)",x="Predicted (%)")+
  ggtitle("Predicting ADF from intact-leaf spectra")+guides(color=F)

perC_jack_pred<-apply.coefs(perC_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
perC_jack_stat<-t(apply(perC_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df<-data.frame(pred_mean=perC_jack_stat[,1],
                         pred_low=perC_jack_stat[,2],
                         pred_high=perC_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perC,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perC_intact_val_plot<-ggplot(perC_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured (%)",x="Predicted (%)")+
  ggtitle(expression("Predicting C"[mass]*" from intact-leaf spectra"))+guides(color=F)

perN_jack_pred<-apply.coefs(perN_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
perN_jack_stat<-t(apply(perN_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df<-data.frame(pred_mean=perN_jack_stat[,1],
                         pred_low=perN_jack_stat[,2],
                         pred_high=perN_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perN,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perN_intact_val_plot<-ggplot(perN_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured (%)",x="Predicted (%)")+
  ggtitle(expression("Predicting N"[mass]*" from intact-leaf spectra"))+guides(color=F)

perC_area_jack_pred<-apply.coefs(perC_area_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
perC_area_jack_stat<-t(apply(perC_area_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_area_jack_df<-data.frame(pred_mean=perC_area_jack_stat[,1],
                         pred_low=perC_area_jack_stat[,2],
                         pred_high=perC_area_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perC_area,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perC_area_intact_val_plot<-ggplot(perC_area_jack_df,
                                  aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,200),ylim=c(0,200))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured (g/m"^2*")"),x=expression("Predicted (g/m"^2*")"))+
  ggtitle(expression("Predicting C"[area]*" from intact-leaf spectra"))

perN_area_jack_pred<-apply.coefs(perN_area_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
perN_area_jack_stat<-t(apply(perN_area_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_area_jack_df<-data.frame(pred_mean=perN_area_jack_stat[,1],
                         pred_low=perN_area_jack_stat[,2],
                         pred_high=perN_area_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perN_area,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perN_area_intact_val_plot<-ggplot(perN_area_jack_df,
                                  aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured (g/m"^2*")"),x=expression("Predicted (g/m"^2*")"))+
  ggtitle(expression("Predicting N"[area]*" from intact-leaf spectra"))+guides(color=F)

LMA_jack_pred<-apply.coefs(LMA_jack_coefs_pressed,as.matrix(intact_spec_agg_test))
LMA_jack_stat<-t(apply(LMA_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df<-data.frame(pred_mean=LMA_jack_stat[,1],
                        pred_low=LMA_jack_stat[,2],
                        pred_high=LMA_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$LMA,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

LMA_intact_val_plot<-ggplot(LMA_jack_df,
                            aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured (g/m"^2*")"),x=expression("Predicted (g/m"^2*")"))+
  ggtitle("Predicting LMA from intact-leaf spectra")

R2.df<-data.frame(NDF=unlist(lapply(NDF_jack_stats,function(x) x[["R2"]])),
                  ADF=unlist(lapply(ADF_jack_stats,function(x) x[["R2"]])),
                  perC=unlist(lapply(perC_jack_stats,function(x) x[["R2"]])),
                  perN=unlist(lapply(perN_jack_stats,function(x) x[["R2"]])),
                  perC_area=unlist(lapply(perC_area_jack_stats,function(x) x[["R2"]])),
                  perN_area=unlist(lapply(perN_area_jack_stats,function(x) x[["R2"]])),
                  LMA=unlist(lapply(LMA_jack_stats,function(x) x[["R2"]])))

R2.long<-melt(R2.df)
intact_val_R2<-ggplot(R2.long,aes(y=value,x=variable))+
  geom_violin()+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="R2",x="Trait")+
  ggtitle("Intact-leaf spectra")

perRMSE.df<-data.frame(NDF=unlist(lapply(NDF_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       ADF=unlist(lapply(ADF_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       perC=unlist(lapply(perC_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       perN=unlist(lapply(perN_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       perC_area=unlist(lapply(perC_area_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       perN_area=unlist(lapply(perN_area_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                       LMA=unlist(lapply(LMA_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))))

perRMSE.long<-melt(perRMSE.df)
intact_val_perRMSE<-ggplot(perRMSE.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20))+
  labs(y="%RMSE",x="Trait")
