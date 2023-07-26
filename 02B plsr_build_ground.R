setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")
library(spectrolab)
library(pls)
library(ggplot2)
library(reshape2)

#########################################
## read data

ground_spec_agg_train<-readRDS("SavedResults/ground_spec_agg_train.rds")
ground_spec_agg_test<-readRDS("SavedResults/ground_spec_agg_test.rds")

source("Scripts/senesced-trait-models/useful_functions.R")

############################################
## build initial PLSR models for calibration data

sol_ground<-plsr(meta(ground_spec_agg_train)$solubles~as.matrix(ground_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
hemi_ground<-plsr(meta(ground_spec_agg_train)$hemicellulose~as.matrix(ground_spec_agg_train),
            ncomp=30,method = "oscorespls",validation="CV",segments=10)
recalc_ground<-plsr(meta(ground_spec_agg_train)$recalcitrant~as.matrix(ground_spec_agg_train),
            ncomp=30,method = "oscorespls",validation="CV",segments=10)
Cmass_ground<-plsr(meta(ground_spec_agg_train)$Cmass~as.matrix(ground_spec_agg_train),
           ncomp=30,method = "oscorespls",validation="CV",segments=10)
Nmass_ground<-plsr(meta(ground_spec_agg_train)$Nmass~as.matrix(ground_spec_agg_train),
           ncomp=30,method = "oscorespls",validation="CV",segments=10)
LMA_ground<-plsr(meta(ground_spec_agg_train)$LMA~as.matrix(ground_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

# recalc_N_ratio_ground<-plsr(meta(ground_spec_agg_train)$recalcitrant/meta(ground_spec_agg_train)$Nmass~as.matrix(ground_spec_agg_train),
#            ncomp=30,validation="LOO")

############################################
## calibration data model fits, figures
## select the number of components for later analyses

ncomp_sol_ground <- selectNcomp(sol_ground, method = "onesigma", plot = FALSE)
sol_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$solubles))
sol_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[sol_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[sol_ground_valid],
                             Run=meta(ground_spec_agg_train)$FiberRun[sol_ground_valid],
                             Measured=meta(ground_spec_agg_train)$solubles[sol_ground_valid],
                             val_pred=sol_ground$validation$pred[,,ncomp_sol_ground])
ggplot(sol_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % solubles from ground-leaf spectra")

ncomp_hemi_ground <- selectNcomp(hemi_ground, method = "onesigma", plot = FALSE)
hemi_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$hemicellulose))
hemi_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[hemi_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[hemi_ground_valid],
                             Run=meta(ground_spec_agg_train)$FiberRun[hemi_ground_valid],
                             Measured=meta(ground_spec_agg_train)$hemicellulose[hemi_ground_valid],
                             val_pred=hemi_ground$validation$pred[,,ncomp_hemi_ground])
ggplot(hemi_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % hemicellulose from ground-leaf spectra")

ncomp_recalc_ground <- selectNcomp(recalc_ground, method = "onesigma", plot = FALSE)
recalc_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$recalcitrant))
recalc_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[recalc_ground_valid],
                            Species=meta(ground_spec_agg_train)$sp[recalc_ground_valid],
                            Run=meta(ground_spec_agg_train)$FiberRun[recalc_ground_valid],
                            Measured=meta(ground_spec_agg_train)$recalcitrant[recalc_ground_valid],
                            val_pred=recalc_ground$validation$pred[,,ncomp_recalc_ground])
ggplot(recalc_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting % recalcitrant from ground-leaf spectra")+guides(color=F)

ncomp_Cmass_ground <- selectNcomp(Cmass_ground, method = "onesigma", plot = FALSE)
Cmass_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$Cmass))
Cmass_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[Cmass_ground_valid],
                              Species=meta(ground_spec_agg_train)$sp[Cmass_ground_valid],
                              Run=meta(ground_spec_agg_train)$EARun[Cmass_ground_valid],
                              Measured=meta(ground_spec_agg_train)$Cmass[Cmass_ground_valid],
                              val_pred=Cmass_ground$validation$pred[,,ncomp_Cmass_ground])
ggplot(Cmass_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %C from ground-leaf spectra") #+guides(color=F)

ncomp_Nmass_ground <- selectNcomp(Nmass_ground, method = "onesigma", plot = FALSE)
Nmass_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$Nmass))
Nmass_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[Nmass_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[Nmass_ground_valid],
                             Run=meta(ground_spec_agg_train)$EARun[Nmass_ground_valid],
                             Measured=meta(ground_spec_agg_train)$Nmass[Nmass_ground_valid],
                             val_pred=Nmass_ground$validation$pred[,,ncomp_Nmass_ground])
ggplot(Nmass_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %N from ground-leaf spectra") # +guides(color=F)

ncomp_LMA_ground <- selectNcomp(LMA_ground, method = "onesigma", plot = FALSE)
LMA_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$LMA))
LMA_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[LMA_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[LMA_ground_valid],
                             Run=meta(ground_spec_agg_train)$EARun[LMA_ground_valid],
                             Measured=meta(ground_spec_agg_train)$LMA[LMA_ground_valid],
                             val_pred=LMA_ground$validation$pred[,,ncomp_LMA_ground])
ggplot(LMA_ground_pred,aes(x=Measured,y=val_pred,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting LMA from ground-leaf spectra")+guides(color=F)

##########################################
## VIP plots

source("VIP.R")

VIP_ground<-data.frame(sol=VIP(sol_ground)[ncomp_sol_ground,],
                       hemi=VIP(hemi_ground)[ncomp_hemi_ground,],
                       recalc=VIP(recalc_ground)[ncomp_recalc_ground,],
                       Cmass=VIP(Cmass_ground)[ncomp_Cmass_ground,],
                       Nmass=VIP(Nmass_ground)[ncomp_Nmass_ground,],
                       LMA=VIP(LMA_ground)[ncomp_LMA_ground,],
                       wavelength=400:2400)
saveRDS(VIP_ground,"SavedResults/VIP_ground.rds")

#######################################
## jackknife tests + prediction of validation data

sol_jack_coefs<-list()
hemi_jack_coefs<-list()
recalc_jack_coefs<-list()
Cmass_jack_coefs<-list()
Nmass_jack_coefs<-list()
LMA_jack_coefs<-list()

sol_jack_stats<-list()
hemi_jack_stats<-list()
recalc_jack_stats<-list()
Cmass_jack_stats<-list()
Nmass_jack_stats<-list()
LMA_jack_stats<-list()
nreps<-200

for(i in 1:nreps){
  print(i)
  
  n_cal_spec<-nrow(ground_spec_agg_train)
  train_jack<-sample(1:n_cal_spec,floor(0.7*n_cal_spec))
  test_jack<-setdiff(1:n_cal_spec,train_jack)
  
  calib_jack<-ground_spec_agg_train[train_jack]
  val_jack<-ground_spec_agg_train[test_jack]

  sol_ground_jack<-plsr(meta(calib_jack)$solubles~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  hemi_ground_jack<-plsr(meta(calib_jack)$hemicellulose~as.matrix(calib_jack),
                   ncomp=30,method = "oscorespls",validation="none")
  recalc_ground_jack<-plsr(meta(calib_jack)$recalcitrant~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  Cmass_ground_jack<-plsr(meta(calib_jack)$Cmass~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  Nmass_ground_jack<-plsr(meta(calib_jack)$Nmass~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  LMA_ground_jack<-plsr(meta(calib_jack)$LMA~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")

  sol_jack_val_pred<-as.vector(predict(sol_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_sol_ground)[,,1])
  sol_jack_val_fit<-lm(sol_jack_val_pred~meta(val_jack)$solubles)
  sol_jack_stats[[i]]<-c(R2=summary(sol_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$solubles,sol_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$solubles,sol_jack_val_pred,0.025,0.975),
                          bias=mean(sol_jack_val_pred,na.rm=T)-mean(meta(val_jack)$solubles,na.rm=T))
  
  hemi_jack_val_pred<-as.vector(predict(hemi_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_hemi_ground)[,,1])
  hemi_jack_val_fit<-lm(hemi_jack_val_pred~meta(val_jack)$hemicellulose)
  hemi_jack_stats[[i]]<-c(R2=summary(hemi_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$hemicellulose,hemi_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$hemicellulose,hemi_jack_val_pred,0.025,0.975),
                         bias=mean(hemi_jack_val_pred,na.rm=T)-mean(meta(val_jack)$hemicellulose,na.rm=T))
  
  recalc_jack_val_pred<-as.vector(predict(recalc_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_recalc_ground)[,,1])
  recalc_jack_val_fit<-lm(recalc_jack_val_pred~meta(val_jack)$recalcitrant)
  recalc_jack_stats[[i]]<-c(R2=summary(recalc_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$recalcitrant,recalc_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$recalcitrant,recalc_jack_val_pred,0.025,0.975),
                         bias=mean(recalc_jack_val_pred,na.rm=T)-mean(meta(val_jack)$recalcitrant,na.rm=T))
  
  Cmass_jack_val_pred<-as.vector(predict(Cmass_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Cmass_ground)[,,1])
  Cmass_jack_val_fit<-lm(Cmass_jack_val_pred~meta(val_jack)$Cmass)
  Cmass_jack_stats[[i]]<-c(R2=summary(Cmass_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$Cmass,Cmass_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$Cmass,Cmass_jack_val_pred,0.025,0.975),
                          bias=mean(Cmass_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Cmass,na.rm=T))
  
  Nmass_jack_val_pred<-as.vector(predict(Nmass_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_Nmass_ground)[,,1])
  Nmass_jack_val_fit<-lm(Nmass_jack_val_pred~meta(val_jack)$Nmass)
  Nmass_jack_stats[[i]]<-c(R2=summary(Nmass_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$Nmass,Nmass_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$Nmass,Nmass_jack_val_pred,0.025,0.975),
                          bias=mean(Nmass_jack_val_pred,na.rm=T)-mean(meta(val_jack)$Nmass,na.rm=T))
  
  LMA_jack_val_pred<-as.vector(predict(LMA_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_LMA_ground)[,,1])
  LMA_jack_val_fit<-lm(LMA_jack_val_pred~meta(val_jack)$LMA)
  LMA_jack_stats[[i]]<-c(R2=summary(LMA_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$LMA,LMA_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$LMA,LMA_jack_val_pred,0.025,0.975),
                         bias=mean(LMA_jack_val_pred,na.rm=T)-mean(meta(val_jack)$LMA,na.rm=T))
  
  sol_jack_coefs[[i]]<-as.vector(coef(sol_ground_jack,ncomp=ncomp_sol_ground,intercept=TRUE))
  hemi_jack_coefs[[i]]<-as.vector(coef(hemi_ground_jack,ncomp=ncomp_hemi_ground,intercept=TRUE))
  recalc_jack_coefs[[i]]<-as.vector(coef(recalc_ground_jack,ncomp=ncomp_recalc_ground,intercept=TRUE))
  Cmass_jack_coefs[[i]]<-as.vector(coef(Cmass_ground_jack,ncomp=ncomp_Cmass_ground,intercept=TRUE))
  Nmass_jack_coefs[[i]]<-as.vector(coef(Nmass_ground_jack,ncomp=ncomp_Nmass_ground,intercept=TRUE))
  LMA_jack_coefs[[i]]<-as.vector(coef(LMA_ground_jack,ncomp=ncomp_LMA_ground,intercept=TRUE))
  
}

sol_jack_pred<-apply.coefs(sol_jack_coefs,as.matrix(ground_spec_agg_test))
sol_jack_stat<-t(apply(sol_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
sol_jack_df<-data.frame(pred_mean=sol_jack_stat[,1],
                         pred_low=sol_jack_stat[,2],
                         pred_high=sol_jack_stat[,3],
                         Measured=meta(ground_spec_agg_test)$solubles,
                         ncomp=ncomp_sol_ground,
                         Species=meta(ground_spec_agg_test)$sp,
                         ID=meta(ground_spec_agg_test)$ID)

hemi_jack_pred<-apply.coefs(hemi_jack_coefs,as.matrix(ground_spec_agg_test))
hemi_jack_stat<-t(apply(hemi_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemi_jack_df<-data.frame(pred_mean=hemi_jack_stat[,1],
                        pred_low=hemi_jack_stat[,2],
                        pred_high=hemi_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$hemicellulose,
                        ncomp=ncomp_hemi_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

recalc_jack_pred<-apply.coefs(recalc_jack_coefs,as.matrix(ground_spec_agg_test))
recalc_jack_stat<-t(apply(recalc_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
recalc_jack_df<-data.frame(pred_mean=recalc_jack_stat[,1],
                        pred_low=recalc_jack_stat[,2],
                        pred_high=recalc_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$recalcitrant,
                        ncomp=ncomp_recalc_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

Cmass_jack_pred<-apply.coefs(Cmass_jack_coefs,as.matrix(ground_spec_agg_test))
Cmass_jack_stat<-t(apply(Cmass_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass_jack_df<-data.frame(pred_mean=Cmass_jack_stat[,1],
                         pred_low=Cmass_jack_stat[,2],
                         pred_high=Cmass_jack_stat[,3],
                         Measured=meta(ground_spec_agg_test)$Cmass,
                         ncomp=ncomp_Cmass_ground,
                         Species=meta(ground_spec_agg_test)$sp,
                         ID=meta(ground_spec_agg_test)$ID)

Nmass_jack_pred<-apply.coefs(Nmass_jack_coefs,as.matrix(ground_spec_agg_test))
Nmass_jack_stat<-t(apply(Nmass_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass_jack_df<-data.frame(pred_mean=Nmass_jack_stat[,1],
                        pred_low=Nmass_jack_stat[,2],
                        pred_high=Nmass_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$Nmass,
                        ncomp=ncomp_Nmass_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

LMA_jack_pred<-apply.coefs(LMA_jack_coefs,as.matrix(ground_spec_agg_test))
LMA_jack_stat<-t(apply(LMA_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df<-data.frame(pred_mean=LMA_jack_stat[,1],
                        pred_low=LMA_jack_stat[,2],
                        pred_high=LMA_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$LMA,
                        ncomp=ncomp_LMA_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

###################################
## save jackknife output

ground_jack_coef_list<-list(LMA=LMA_jack_coefs,
                            Cmass=Cmass_jack_coefs,
                            Nmass=Nmass_jack_coefs,
                            sol=sol_jack_coefs,
                            hemi=hemi_jack_coefs,
                            recalc=recalc_jack_coefs)
saveRDS(ground_jack_coef_list,"SavedResults/ground_jack_coefs_list.rds")

ground_jack_df_list<-list(LMA=LMA_jack_df,
                          Cmass=Cmass_jack_df,
                          Nmass=Nmass_jack_df,
                          sol=sol_jack_df,
                          hemi=hemi_jack_df,
                          recalc=recalc_jack_df)
saveRDS(ground_jack_df_list,"SavedResults/ground_jack_df_list.rds")

summ<-data.frame(ncomp=unlist(lapply(ground_jack_df_list,function(x) x$ncomp[1])),
                 r2=round(unlist(lapply(ground_jack_df_list,function(x) summary(lm(Measured~pred_mean,data=x))$r.squared)),3),
                 rmse=signif(unlist(lapply(ground_jack_df_list,function(x) RMSD(x$Measured,x$pred_mean))),3),
                 perrmse=signif(unlist(lapply(ground_jack_df_list,function(x) percentRMSD(x$Measured,x$pred_mean,0.025,0.975)))*100,3))
write.csv(summ,"SavedResults/stat_summary.csv")

############################################
## violin plots

R2.df<-data.frame(sol=unlist(lapply(sol_jack_stats,function(x) x[["R2"]])),
                  hemi=unlist(lapply(hemi_jack_stats,function(x) x[["R2"]])),
                  recalc=unlist(lapply(recalc_jack_stats,function(x) x[["R2"]])),
                  Cmass=unlist(lapply(Cmass_jack_stats,function(x) x[["R2"]])),
                  Nmass=unlist(lapply(Nmass_jack_stats,function(x) x[["R2"]])),
                  LMA=unlist(lapply(LMA_jack_stats,function(x) x[["R2"]])))

R2.long<-melt(R2.df)
ground_val_R2<-ggplot(R2.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y=expression(italic("R"^2)),x="Trait")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df<-data.frame(sol=unlist(lapply(sol_jack_stats,function(x) 100*x[["perRMSE"]])),
                       hemi=unlist(lapply(hemi_jack_stats,function(x) 100*x[["perRMSE"]])),
                       recalc=unlist(lapply(recalc_jack_stats,function(x) 100*x[["perRMSE"]])),
                       Cmass=unlist(lapply(Cmass_jack_stats,function(x) 100*x[["perRMSE"]])),
                       Nmass=unlist(lapply(Nmass_jack_stats,function(x) 100*x[["perRMSE"]])),
                       LMA=unlist(lapply(LMA_jack_stats,function(x) 100*x[["perRMSE"]])))

perRMSE.long<-melt(perRMSE.df)
ground_val_perRMSE<-ggplot(perRMSE.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE",x="Trait")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,max(perRMSE.long$value)*1.1))

pdf("Manuscript/FigS2.pdf",height=8,width=4.5)
(ground_val_R2/ground_val_perRMSE)
dev.off()

#######################################
## compare with spectra from the B4WARMED decomp experiment

source("Decomp/process_decomp.R")
decomp_chem<-read.csv("Decomp/stoich_leaves.csv")
colnames(decomp_chem)<-c("site","habitat","trt","species","Nmass","Cmass",
                         "sol","NDF","hemi","recalc","cellulose","ADL")
decomp_chem$habitat<-toupper(decomp_chem$habitat)
decomp_chem$trt<-toupper(decomp_chem$trt)
decomp_chem$species<-toupper(decomp_chem$species)
decomp_chem$full_id<-apply(decomp_chem[,1:4],1,paste,collapse="_")
match_ids_decomp<-match(decomp_chem$full_id,meta(decomp_agg)$full_id)

decomp_chem$recalc_N<-decomp_chem$recalc/decomp_chem$Nmass
decomp_chem$ADL_N<-decomp_chem$ADL/decomp_chem$Nmass

Cmass_jack_pred_decomp<-apply.coefs(Cmass_jack_coefs,as.matrix(decomp_agg[,400:2400]))
decomp_chem$predC<-rowMeans(Cmass_jack_pred_decomp)[match_ids_decomp]

Nmass_jack_pred_decomp<-apply.coefs(Nmass_jack_coefs,as.matrix(decomp_agg[,400:2400]))
decomp_chem$predN<-rowMeans(Nmass_jack_pred_decomp)[match_ids_decomp]

sol_jack_pred_decomp<-apply.coefs(sol_jack_coefs,as.matrix(decomp_agg[,400:2400]))
decomp_chem$pred_sol<-rowMeans(sol_jack_pred_decomp)[match_ids_decomp]

hemi_jack_pred_decomp<-apply.coefs(hemi_jack_coefs,as.matrix(decomp_agg[,400:2400]))
decomp_chem$pred_hemi<-rowMeans(hemi_jack_pred_decomp)[match_ids_decomp]

recalc_jack_pred_decomp<-apply.coefs(recalc_jack_coefs,as.matrix(decomp_agg[,400:2400]))
decomp_chem$pred_recalc<-rowMeans(recalc_jack_pred_decomp)[match_ids_decomp]

ggplot(decomp_chem,aes(x=predN,y=Nmass))+
  geom_point(aes(color=species,shape=trt),size=2)+
  geom_smooth(method="lm",se=F,size=2,color="black")+
  geom_smooth(method="lm",se=F,aes(color=species))+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+theme_bw()+
  theme(text = element_text(size=20))+
  labs(x="Predicted N",y="Measured N")+
  scale_color_brewer(palette="Set2")+
  coord_cartesian(xlim=c(0.2,1.6),ylim=c(0.2,1.6))

ggplot(decomp_chem,aes(x=pred_sol,y=sol))+
  geom_point(aes(color=species,shape=trt),size=2)+
  geom_smooth(method="lm",se=F,size=2,color="black")+
  geom_smooth(method="lm",se=F,aes(color=species))+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+theme_bw()+
  theme(text = element_text(size=20))+
  labs(x="Predicted solubles",y="Measured solubles")+
  scale_color_brewer(palette="Set2")+
  coord_cartesian(xlim=c(10,70),ylim=c(10,70))

ggplot(decomp_chem,aes(x=pred_recalc,y=recalc))+
  geom_point(aes(color=species,shape=trt),size=2)+
  geom_smooth(method="lm",se=F,size=2,color="black")+
  geom_smooth(method="lm",se=F,aes(color=species))+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+theme_bw()+
  theme(text = element_text(size=20))+
  labs(x="Predicted recalcitrant",y="Measured recalcitrant")+
  scale_color_brewer(palette="Set2")+
  coord_cartesian(xlim=c(18,65),ylim=c(18,65))

decay_rate<-read.csv("Decomp/decay_rate.csv")
decay_rate$species<-as.factor(decay_rate$species)
decay_rate$tissue_cond<-as.factor(decay_rate$tissue_cond)

levels(decay_rate$species)<-c("ACERU","ACESA","BETPA","PINBA",
                              "PINST","POPTR","QUEMA","QUERU")
levels(decay_rate$tissue_cond)<-c("AMB","WARM")
decay_rate$tissue_cond[is.na(decay_rate$tissue_cond)]<-"AMB"

decay_rate$full_id<-toupper(apply(decay_rate[,c(1,2,4,3)],1,paste,collapse="_"))

decomp_match<-match(decomp_chem$full_id,decay_rate$full_id)
decomp_chem$negexp_decay<-decay_rate$steady_state_Neg.exp[decomp_match]

ggplot(data=decomp_chem,aes(x=Nmass,y=negexp_decay,
                            color=trt,shape=habitat))+
  geom_point()+
  geom_smooth(method="lm")
