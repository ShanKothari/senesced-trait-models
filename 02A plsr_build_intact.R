setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")
library(spectrolab)
library(patchwork)
library(pls)
library(ggplot2)
library(reshape2)

#########################################
## read data

intact_spec_agg_train<-readRDS("SavedResults/intact_spec_agg_train.rds")
intact_spec_agg_test<-readRDS("SavedResults/intact_spec_agg_test.rds")

source("Scripts/senesced-trait-models/useful_functions.R")

#########################################
## train models

sol_intact<-plsr(meta(intact_spec_agg_train)$solubles~as.matrix(intact_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
hemi_intact<-plsr(meta(intact_spec_agg_train)$hemicellulose~as.matrix(intact_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
recalc_intact<-plsr(meta(intact_spec_agg_train)$recalcitrant~as.matrix(intact_spec_agg_train),
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

ncomp_sol_intact <- selectNcomp(sol_intact, method = "onesigma", plot = FALSE)
sol_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$solubles))
sol_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[sol_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[sol_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[sol_intact_valid],
                            Measured=meta(intact_spec_agg_train)$solubles[sol_intact_valid],
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
hemi_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$hemicellulose))
hemi_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[hemi_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[hemi_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[hemi_intact_valid],
                            Measured=meta(intact_spec_agg_train)$hemicellulose[hemi_intact_valid],
                            val_pred=hemi_intact$validation$pred[,,ncomp_hemi_intact])
ggplot(hemi_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % hemicellulose from intact-leaf spectra")+guides(color=F)

ncomp_recalc_intact <- selectNcomp(recalc_intact, method = "onesigma", plot = FALSE)
recalc_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$recalcitrant))
recalc_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[recalc_intact_valid],
                            Species=meta(intact_spec_agg_train)$sp[recalc_intact_valid],
                            Run=meta(intact_spec_agg_train)$FiberRun[recalc_intact_valid],
                            Measured=meta(intact_spec_agg_train)$recalcitrant[recalc_intact_valid],
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

ncomp_perC_intact <- selectNcomp(perC_intact, method = "onesigma", plot = FALSE)
perC_intact_valid <- which(!is.na(meta(intact_spec_agg_train)$perC))
perC_intact_pred<-data.frame(ID=meta(intact_spec_agg_train)$ID[perC_intact_valid],
                             Species=meta(intact_spec_agg_train)$sp[perC_intact_valid],
                             Run=meta(intact_spec_agg_train)$FiberRun[perC_intact_valid],
                             Measured=meta(intact_spec_agg_train)$perC[perC_intact_valid],
                             val_pred=perC_intact$validation$pred[,,ncomp_perC_intact])
ggplot(perC_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
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
ggplot(perN_intact_pred,aes(x=Measured,y=val_pred,color=Species))+
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
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
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
  coord_cartesian(xlim=c(0,200),ylim=c(0,200))+
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
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured (g N/cm2)",y="Predicted (g N/cm2)")+
  ggtitle("Predicting N per area from intact-leaf spectra")

##########################################
## VIP plots

source("VIP.R")
VIP_intact<-data.frame(sol=VIP(sol_intact)[ncomp_sol_intact,],
                       hemi=VIP(hemi_intact)[ncomp_hemi_intact,],
                       recalc=VIP(recalc_intact)[ncomp_recalc_intact,],
                       perC=VIP(perC_intact)[ncomp_perC_intact,],
                       perN=VIP(perN_intact)[ncomp_perN_intact,],
                       LMA=VIP(LMA_intact)[ncomp_LMA_intact,],
                       wavelength=400:2400)
saveRDS(VIP_intact,"SavedResults/VIP_intact.rds")

#######################################
## jackknife tests + prediction of validation data

sol_jack_coefs<-list()
hemi_jack_coefs<-list()
recalc_jack_coefs<-list()
perC_jack_coefs<-list()
perN_jack_coefs<-list()
perC_area_jack_coefs<-list()
perN_area_jack_coefs<-list()
LMA_jack_coefs<-list()

sol_jack_stats<-list()
hemi_jack_stats<-list()
recalc_jack_stats<-list()
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
  
  sol_intact_jack<-plsr(meta(calib_jack)$solubles~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  hemi_intact_jack<-plsr(meta(calib_jack)$hemicellulose~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  recalc_intact_jack<-plsr(meta(calib_jack)$recalcitrant~as.matrix(calib_jack),
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
  
  sol_jack_coefs[[i]]<-as.vector(coef(sol_intact_jack,ncomp=ncomp_sol_intact,intercept=TRUE))
  hemi_jack_coefs[[i]]<-as.vector(coef(hemi_intact_jack,ncomp=ncomp_hemi_intact,intercept=TRUE))
  recalc_jack_coefs[[i]]<-as.vector(coef(recalc_intact_jack,ncomp=ncomp_recalc_intact,intercept=TRUE))
  perC_jack_coefs[[i]]<-as.vector(coef(perC_intact_jack,ncomp=ncomp_perC_intact,intercept=TRUE))
  perN_jack_coefs[[i]]<-as.vector(coef(perN_intact_jack,ncomp=ncomp_perN_intact,intercept=TRUE))
  perC_area_jack_coefs[[i]]<-as.vector(coef(perC_area_intact_jack,ncomp=ncomp_perC_area_intact,intercept=TRUE))
  perN_area_jack_coefs[[i]]<-as.vector(coef(perN_area_intact_jack,ncomp=ncomp_perN_area_intact,intercept=TRUE))
  LMA_jack_coefs[[i]]<-as.vector(coef(LMA_intact_jack,ncomp=ncomp_LMA_intact,intercept=TRUE))
  
}

sol_jack_pred<-apply.coefs(sol_jack_coefs,as.matrix(intact_spec_agg_test))
sol_jack_stat<-t(apply(sol_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
sol_jack_df<-data.frame(pred_mean=sol_jack_stat[,1],
                        pred_low=sol_jack_stat[,2],
                        pred_high=sol_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$solubles,
                        ncomp=ncomp_sol_intact,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

hemi_jack_pred<-apply.coefs(hemi_jack_coefs,as.matrix(intact_spec_agg_test))
hemi_jack_stat<-t(apply(hemi_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemi_jack_df<-data.frame(pred_mean=hemi_jack_stat[,1],
                        pred_low=hemi_jack_stat[,2],
                        pred_high=hemi_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$hemicellulose,
                        ncomp=ncomp_hemi_intact,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

recalc_jack_pred<-apply.coefs(recalc_jack_coefs,as.matrix(intact_spec_agg_test))
recalc_jack_stat<-t(apply(recalc_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
recalc_jack_df<-data.frame(pred_mean=recalc_jack_stat[,1],
                        pred_low=recalc_jack_stat[,2],
                        pred_high=recalc_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$recalcitrant,
                        ncomp=ncomp_recalc_intact,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

perC_jack_pred<-apply.coefs(perC_jack_coefs,as.matrix(intact_spec_agg_test))
perC_jack_stat<-t(apply(perC_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df<-data.frame(pred_mean=perC_jack_stat[,1],
                         pred_low=perC_jack_stat[,2],
                         pred_high=perC_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perC,
                         ncomp=ncomp_perC_intact,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perN_jack_pred<-apply.coefs(perN_jack_coefs,as.matrix(intact_spec_agg_test))
perN_jack_stat<-t(apply(perN_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df<-data.frame(pred_mean=perN_jack_stat[,1],
                         pred_low=perN_jack_stat[,2],
                         pred_high=perN_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perN,
                         ncomp=ncomp_perN_intact,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perC_area_jack_pred<-apply.coefs(perC_area_jack_coefs,as.matrix(intact_spec_agg_test))
perC_area_jack_stat<-t(apply(perC_area_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_area_jack_df<-data.frame(pred_mean=perC_area_jack_stat[,1],
                         pred_low=perC_area_jack_stat[,2],
                         pred_high=perC_area_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perC_area,
                         ncomp=ncomp_perC_area_intact,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

perN_area_jack_pred<-apply.coefs(perN_area_jack_coefs,as.matrix(intact_spec_agg_test))
perN_area_jack_stat<-t(apply(perN_area_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_area_jack_df<-data.frame(pred_mean=perN_area_jack_stat[,1],
                         pred_low=perN_area_jack_stat[,2],
                         pred_high=perN_area_jack_stat[,3],
                         Measured=meta(intact_spec_agg_test)$perN_area,
                         ncomp=ncomp_perN_area_intact,
                         Species=meta(intact_spec_agg_test)$sp,
                         ID=meta(intact_spec_agg_test)$ID)

LMA_jack_pred<-apply.coefs(LMA_jack_coefs,as.matrix(intact_spec_agg_test))
LMA_jack_stat<-t(apply(LMA_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df<-data.frame(pred_mean=LMA_jack_stat[,1],
                        pred_low=LMA_jack_stat[,2],
                        pred_high=LMA_jack_stat[,3],
                        Measured=meta(intact_spec_agg_test)$LMA,
                        ncomp=ncomp_LMA_intact,
                        Species=meta(intact_spec_agg_test)$sp,
                        ID=meta(intact_spec_agg_test)$ID)

###################################
## save jackknife output

intact_jack_coef_list<-list(LMA=LMA_jack_coefs,
                            perC=perC_jack_coefs,
                            perN=perN_jack_coefs,
                            sol=sol_jack_coefs,
                            hemi=hemi_jack_coefs,
                            recalc=recalc_jack_coefs,
                            perC_area=perC_area_jack_coefs,
                            perN_area=perN_area_jack_coefs)
saveRDS(intact_jack_coef_list,"SavedResults/intact_jack_coefs_list.rds")

intact_jack_df_list<-list(LMA=LMA_jack_df,
                          perC=perC_jack_df,
                          perN=perN_jack_df,
                          sol=sol_jack_df,
                          hemi=hemi_jack_df,
                          recalc=recalc_jack_df,
                          perC_area=perC_area_jack_df,
                          perN_area=perN_area_jack_df)
saveRDS(intact_jack_df_list,"SavedResults/intact_jack_df_list.rds")

############################################
## violin plots

R2.df<-data.frame(sol=unlist(lapply(sol_jack_stats,function(x) x[["R2"]])),
                  hemi=unlist(lapply(hemi_jack_stats,function(x) x[["R2"]])),
                  recalc=unlist(lapply(recalc_jack_stats,function(x) x[["R2"]])),
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
  labs(y=expression(italic("R"^2)),x="Trait")+
  ggtitle("Intact-leaf spectra")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df<-data.frame(sol=unlist(lapply(sol_jack_stats,function(x) 100*x[["perRMSE"]])),
                       hemi=unlist(lapply(hemi_jack_stats,function(x) 100*x[["perRMSE"]])),
                       recalc=unlist(lapply(recalc_jack_stats,function(x) 100*x[["perRMSE"]])),
                       perC=unlist(lapply(perC_jack_stats,function(x) 100*x[["perRMSE"]])),
                       perN=unlist(lapply(perN_jack_stats,function(x) 100*x[["perRMSE"]])),
                       perC_area=unlist(lapply(perC_area_jack_stats,function(x) 100*x[["perRMSE"]])),
                       perN_area=unlist(lapply(perN_area_jack_stats,function(x) 100*x[["perRMSE"]])),
                       LMA=unlist(lapply(LMA_jack_stats,function(x) 100*x[["perRMSE"]])))

perRMSE.long<-melt(perRMSE.df)
intact_val_perRMSE<-ggplot(perRMSE.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE",x="Trait")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,max(perRMSE.long$value)*1.1))

pdf("Manuscript/FigS1.pdf",height=8,width=8)
(intact_val_R2/intact_val_perRMSE)
dev.off()
