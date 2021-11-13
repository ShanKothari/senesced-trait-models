setwd("C:/Users/kotha020/Dropbox/TraitModels2018/SenescencePaper/")
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

NDF_ground<-plsr(meta(ground_spec_agg_train)$NDF~as.matrix(ground_spec_agg_train),
            ncomp=30,method = "oscorespls",validation="CV",segments=10)
# NDF_ground<-plsr(meta(ground_spec_agg_train)$NDF~log(1/as.matrix(ground_spec_agg_train)),
#                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ADF_ground<-plsr(meta(ground_spec_agg_train)$ADF~as.matrix(ground_spec_agg_train),
            ncomp=30,method = "oscorespls",validation="CV",segments=10)
# ADF_ground<-plsr(meta(ground_spec_agg_train)$ADF~log(1/as.matrix(ground_spec_agg_train)),
#                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
perC_ground<-plsr(meta(ground_spec_agg_train)$perC~as.matrix(ground_spec_agg_train),
           ncomp=30,method = "oscorespls",validation="CV",segments=10)
perN_ground<-plsr(meta(ground_spec_agg_train)$perN~as.matrix(ground_spec_agg_train),
           ncomp=30,method = "oscorespls",validation="CV",segments=10)
LMA_ground<-plsr(meta(ground_spec_agg_train)$LMA~as.matrix(ground_spec_agg_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

# ADFN_ground<-plsr(meta(ground_spec_agg_train)$ADF/meta(ground_spec_agg_train)$perN~as.matrix(ground_spec_agg_train),
#            ncomp=30,validation="LOO")

############################################
## calibration data model fits, figures
## select the number of components for later analyses

ncomp_NDF_ground <- selectNcomp(NDF_ground, method = "onesigma", plot = FALSE)
NDF_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$NDF))
NDF_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[NDF_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[NDF_ground_valid],
                             Run=meta(ground_spec_agg_train)$FiberRun[NDF_ground_valid],
                             Measured=meta(ground_spec_agg_train)$NDF[NDF_ground_valid],
                             val_pred=NDF_ground$validation$pred[,,ncomp_NDF_ground])
ggplot(NDF_ground_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %NDF from ground-leaf spectra")

ncomp_ADF_ground <- selectNcomp(ADF_ground, method = "onesigma", plot = FALSE)
ADF_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$ADF))
ADF_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[ADF_ground_valid],
                            Species=meta(ground_spec_agg_train)$sp[ADF_ground_valid],
                            Run=meta(ground_spec_agg_train)$FiberRun[ADF_ground_valid],
                            Measured=meta(ground_spec_agg_train)$ADF[ADF_ground_valid],
                            val_pred=ADF_ground$validation$pred[,,ncomp_ADF_ground])
ggplot(ADF_ground_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %ADF from ground-leaf spectra")+guides(color=F)

ncomp_perC_ground <- selectNcomp(perC_ground, method = "onesigma", plot = FALSE)
perC_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$perC))
perC_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[perC_ground_valid],
                              Species=meta(ground_spec_agg_train)$sp[perC_ground_valid],
                              Run=meta(ground_spec_agg_train)$EARun[perC_ground_valid],
                              Measured=meta(ground_spec_agg_train)$perC[perC_ground_valid],
                              val_pred=perC_ground$validation$pred[,,ncomp_perC_ground])
ggplot(perC_ground_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting %C from ground-leaf spectra") #+guides(color=F)

ncomp_perN_ground <- selectNcomp(perN_ground, method = "onesigma", plot = FALSE)
perN_ground_valid <- which(!is.na(meta(ground_spec_agg_train)$perN))
perN_ground_pred<-data.frame(ID=meta(ground_spec_agg_train)$ID[perN_ground_valid],
                             Species=meta(ground_spec_agg_train)$sp[perN_ground_valid],
                             Run=meta(ground_spec_agg_train)$EARun[perN_ground_valid],
                             Measured=meta(ground_spec_agg_train)$perN[perN_ground_valid],
                             val_pred=perN_ground$validation$pred[,,ncomp_perN_ground])
ggplot(perN_ground_pred,aes(x=Measured*100,y=val_pred*100,color=Species))+
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
  coord_cartesian(xlim=c(0,.03),ylim=c(0,.03))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.3))+
  labs(x="Measured",y="Predicted")+
  ggtitle("Predicting LMA from ground-leaf spectra")+guides(color=F)

##########################################
## VIP plots

source("Senesced_JCB/VIP.R")

VIP_ground<-data.frame(NDF=VIP(NDF_ground)[ncomp_NDF_ground,],
                       ADF=VIP(ADF_ground)[ncomp_ADF_ground,],
                       perC=VIP(perC_ground)[ncomp_perC_ground,],
                       perN=VIP(perN_ground)[ncomp_perN_ground,],
                       LMA=VIP(LMA_ground)[ncomp_LMA_ground,],
                       wavelength=400:2400)
saveRDS(VIP_ground,"VIP_ground.rds")

#######################################
## jackknife tests + prediction of validation data

NDF_jack_coefs<-list()
ADF_jack_coefs<-list()
perC_jack_coefs<-list()
perN_jack_coefs<-list()
LMA_jack_coefs<-list()

NDF_jack_stats<-list()
ADF_jack_stats<-list()
perC_jack_stats<-list()
perN_jack_stats<-list()
LMA_jack_stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n_cal_spec<-nrow(ground_spec_agg_train)
  train_jack<-sample(1:n_cal_spec,floor(0.7*n_cal_spec))
  test_jack<-setdiff(1:n_cal_spec,train_jack)
  
  calib_jack<-ground_spec_agg_train[train_jack]
  val_jack<-ground_spec_agg_train[test_jack]

  NDF_ground_jack<-plsr(meta(calib_jack)$NDF~as.matrix(calib_jack),
                   ncomp=30,method = "oscorespls",validation="none")
  ADF_ground_jack<-plsr(meta(calib_jack)$ADF~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  perC_ground_jack<-plsr(meta(calib_jack)$perC~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  perN_ground_jack<-plsr(meta(calib_jack)$perN~as.matrix(calib_jack),
                        ncomp=30,method = "oscorespls",validation="none")
  LMA_ground_jack<-plsr(meta(calib_jack)$LMA~as.matrix(calib_jack),
                         ncomp=30,method = "oscorespls",validation="none")
  
  NDF_jack_val_pred<-as.vector(predict(NDF_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_NDF_ground)[,,1])
  NDF_jack_val_fit<-lm(NDF_jack_val_pred~meta(val_jack)$NDF)
  NDF_jack_stats[[i]]<-c(R2=summary(NDF_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$NDF,NDF_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$NDF,NDF_jack_val_pred,0.025,0.975),
                         bias=mean(NDF_jack_val_pred,na.rm=T)-mean(meta(val_jack)$NDF,na.rm=T))
  
  ADF_jack_val_pred<-as.vector(predict(ADF_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_ADF_ground)[,,1])
  ADF_jack_val_fit<-lm(ADF_jack_val_pred~meta(val_jack)$ADF)
  ADF_jack_stats[[i]]<-c(R2=summary(ADF_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$ADF,ADF_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$ADF,ADF_jack_val_pred,0.025,0.975),
                         bias=mean(ADF_jack_val_pred,na.rm=T)-mean(meta(val_jack)$ADF,na.rm=T))
  
  perC_jack_val_pred<-as.vector(predict(perC_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perC_ground)[,,1])
  perC_jack_val_fit<-lm(perC_jack_val_pred~meta(val_jack)$perC)
  perC_jack_stats[[i]]<-c(R2=summary(perC_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$perC,perC_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$perC,perC_jack_val_pred,0.025,0.975),
                          bias=mean(perC_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perC,na.rm=T))
  
  perN_jack_val_pred<-as.vector(predict(perN_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_perN_ground)[,,1])
  perN_jack_val_fit<-lm(perN_jack_val_pred~meta(val_jack)$perN)
  perN_jack_stats[[i]]<-c(R2=summary(perN_jack_val_fit)$r.squared,
                          RMSE=RMSD(meta(val_jack)$perN,perN_jack_val_pred),
                          perRMSE=percentRMSD(meta(val_jack)$perN,perN_jack_val_pred,0.025,0.975),
                          bias=mean(perN_jack_val_pred,na.rm=T)-mean(meta(val_jack)$perN,na.rm=T))
  
  LMA_jack_val_pred<-as.vector(predict(LMA_ground_jack,newdata=as.matrix(val_jack),ncomp=ncomp_LMA_ground)[,,1])
  LMA_jack_val_fit<-lm(LMA_jack_val_pred~meta(val_jack)$LMA)
  LMA_jack_stats[[i]]<-c(R2=summary(LMA_jack_val_fit)$r.squared,
                         RMSE=RMSD(meta(val_jack)$LMA,LMA_jack_val_pred),
                         perRMSE=percentRMSD(meta(val_jack)$LMA,LMA_jack_val_pred,0.025,0.975),
                         bias=mean(LMA_jack_val_pred,na.rm=T)-mean(meta(val_jack)$LMA,na.rm=T))
  
  NDF_jack_coefs[[i]]<-as.vector(coef(NDF_ground_jack,ncomp=ncomp_NDF_ground,intercept=TRUE))
  ADF_jack_coefs[[i]]<-as.vector(coef(ADF_ground_jack,ncomp=ncomp_ADF_ground,intercept=TRUE))
  perC_jack_coefs[[i]]<-as.vector(coef(perC_ground_jack,ncomp=ncomp_perC_ground,intercept=TRUE))
  perN_jack_coefs[[i]]<-as.vector(coef(perN_ground_jack,ncomp=ncomp_perN_ground,intercept=TRUE))
  LMA_jack_coefs[[i]]<-as.vector(coef(LMA_ground_jack,ncomp=ncomp_LMA_ground,intercept=TRUE))
  
}

NDF_jack_pred<-apply.coefs(NDF_jack_coefs_pressed,as.matrix(ground_spec_agg_test))
NDF_jack_stat<-t(apply(NDF_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
NDF_jack_df<-data.frame(pred_mean=NDF_jack_stat[,1],
                        pred_low=NDF_jack_stat[,2],
                        pred_high=NDF_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$NDF,
                        ncomp=ncomp_NDF_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

NDF_ground_val_plot<-ggplot(NDF_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured NDF (%)",x="Predicted NDF (%)")+
  guides(color=F)

ADF_jack_pred<-apply.coefs(ADF_jack_coefs_pressed,as.matrix(ground_spec_agg_test))
ADF_jack_stat<-t(apply(ADF_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
ADF_jack_df<-data.frame(pred_mean=ADF_jack_stat[,1],
                        pred_low=ADF_jack_stat[,2],
                        pred_high=ADF_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$ADF,
                        ncomp=ncomp_ADF_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

ADF_ground_val_plot<-ggplot(ADF_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured ADF (%)",x="Predicted ADF (%)")+
  guides(color=F)

perC_jack_pred<-apply.coefs(perC_jack_coefs_pressed,as.matrix(ground_spec_agg_test))
perC_jack_stat<-t(apply(perC_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perC_jack_df<-data.frame(pred_mean=perC_jack_stat[,1],
                         pred_low=perC_jack_stat[,2],
                         pred_high=perC_jack_stat[,3],
                         Measured=meta(ground_spec_agg_test)$perC,
                         ncomp=ncomp_perC_ground,
                         Species=meta(ground_spec_agg_test)$sp,
                         ID=meta(ground_spec_agg_test)$ID)

perC_ground_val_plot<-ggplot(perC_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())++
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)

perN_jack_pred<-apply.coefs(perN_jack_coefs_pressed,as.matrix(ground_spec_agg_test))
perN_jack_stat<-t(apply(perN_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
perN_jack_df<-data.frame(pred_mean=perN_jack_stat[,1],
                        pred_low=perN_jack_stat[,2],
                        pred_high=perN_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$perN,
                        ncomp=ncomp_perN_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

perN_ground_val_plot<-ggplot(perN_jack_df,aes(y=Measured*100,x=pred_mean*100,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())++
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

LMA_jack_pred<-apply.coefs(LMA_jack_coefs_pressed,as.matrix(ground_spec_agg_test))
LMA_jack_stat<-t(apply(LMA_jack_pred,1,function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_jack_df<-data.frame(pred_mean=LMA_jack_stat[,1],
                        pred_low=LMA_jack_stat[,2],
                        pred_high=LMA_jack_stat[,3],
                        Measured=meta(ground_spec_agg_test)$LMA,
                        ncomp=ncomp_LMA_ground,
                        Species=meta(ground_spec_agg_test)$sp,
                        ID=meta(ground_spec_agg_test)$ID)

LMA_ground_val_plot<-ggplot(LMA_jack_df,
                            aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (g/m"^2*")"),
       x=expression("Predicted LMA (g/m"^2*")"))+
  guides(color=F)

###################################
## save jackknife output

ground_jack_coef_list<-list(LMA=LMA_jack_coefs,
                            perC=perC_jack_coefs,
                            perN=perN_jack_coefs,
                            NDF=NDF_jack_coefs,
                            ADF=ADF_jack_coefs)
saveRDS(ground_jack_coef_list,"SavedResults/ground_jack_coefs_list.rds")

ground_jack_df_list<-list(LMA=LMA_jack_df,
                          perC=perC_jack_df,
                          perN=perN_jack_df,
                          NDF=NDF_jack_df,
                          ADF=ADF_jack_df)
saveRDS(ground_jack_df_list,"SavedResults/ground_jack_df_list.rds")

############################################
## violin plots

R2.df<-data.frame(NDF=unlist(lapply(NDF_jack_stats,function(x) x[["R2"]])),
                  ADF=unlist(lapply(ADF_jack_stats,function(x) x[["R2"]])),
                  perC=unlist(lapply(perC_jack_stats,function(x) x[["R2"]])),
                  perN=unlist(lapply(perN_jack_stats,function(x) x[["R2"]])),
                  LMA=unlist(lapply(LMA_jack_stats,function(x) x[["R2"]])))

R2.long<-melt(R2.df)
ground_val_R2<-ggplot(R2.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="R2",x="Trait")+
  ggtitle("Ground-leaf spectra")

perRMSE.df<-data.frame(NDF=unlist(lapply(NDF_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                  ADF=unlist(lapply(ADF_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                  perC=unlist(lapply(perC_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                  perN=unlist(lapply(perN_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))),
                  LMA=unlist(lapply(LMA_jack_stats,function(x) 100*x[["RMSE"]]/(x[["max.val"]]-x[["min.val"]]))))

perRMSE.long<-melt(perRMSE.df)
ground_val_perRMSE<-ggplot(perRMSE.long,aes(y=value,x=variable))+
  geom_violin()+theme_bw()+
  theme(text = element_text(size=20))+
  labs(y="%RMSE",x="Trait")

#######################################
## compare with spectra from the B4WARMED decomp experiment

source("Senesced_JCB/Decomp/process_decomp.R")
decomp_chem<-read.csv("Senesced_JCB/Decomp/stoich_leaves.csv")
colnames(decomp_chem)<-c("site","habitat","trt","species","perN","perC",
                         "solubles","NDF","hemi","ADF","cellulose","ADL")
decomp_chem$habitat<-toupper(decomp_chem$habitat)
decomp_chem$trt<-toupper(decomp_chem$trt)
decomp_chem$species<-toupper(decomp_chem$species)
decomp_chem$full_id<-apply(decomp_chem[,1:4],1,paste,collapse="_")
match_ids_decomp<-match(decomp_chem$full_id,meta(decomp_agg)$full_id)

# N_decomp_predict<-predict(perN_ground,
#                           ncomp=ncomp_perN_ground,
#                           newdata=as.matrix(decomp_agg[,400:2500]))
# decomp_chem$predN<-N_decomp_predict[match_ids_decomp]

perC_jack_pred_decomp<-t(apply(as.matrix(decomp_agg[,400:2500]),1,function(spec) {
  preds<-lapply(perC_jack_coefs,function(coef) coef[1]+sum(coef[-1]*spec))
  return(unlist(preds))
}))
decomp_chem$predC<-rowMeans(perC_jack_pred_decomp)[match_ids_decomp]

perN_jack_pred_decomp<-t(apply(as.matrix(decomp_agg[,400:2500]),1,function(spec) {
  preds<-lapply(perN_jack_coefs,function(coef) coef[1]+sum(coef[-1]*spec))
  return(unlist(preds))
}))
decomp_chem$predN<-rowMeans(perN_jack_pred_decomp)[match_ids_decomp]

NDF_jack_pred_decomp<-t(apply(as.matrix(decomp_agg[,400:2500]),1,function(spec) {
  preds<-lapply(NDF_jack_coefs,function(coef) coef[1]+sum(coef[-1]*spec))
  return(unlist(preds))
}))
decomp_chem$predNDF<-rowMeans(NDF_jack_pred_decomp)[match_ids_decomp]

ADF_jack_pred_decomp<-t(apply(as.matrix(decomp_agg[,400:2500]),1,function(spec) {
  preds<-lapply(ADF_jack_coefs,function(coef) coef[1]+sum(coef[-1]*spec))
  return(unlist(preds))
}))
decomp_chem$predADF<-rowMeans(ADF_jack_pred_decomp)[match_ids_decomp]

ggplot(decomp_chem,aes(x=predN,y=perN/100))+
  geom_point(aes(color=species,shape=trt),size=2)+
  geom_smooth(method="lm",se=F,size=2,color="black")+
  geom_smooth(method="lm",se=F,aes(color=species))+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+theme_bw()+
  theme(text = element_text(size=20))+
  labs(x="Predicted N",y="Measured N")+
  scale_color_brewer(palette="Set2")+
  coord_cartesian(xlim=c(0.003,0.016),ylim=c(0.003,0.016))

decay_rate<-read.csv("Senesced_JCB/Decomp/decay_rate.csv")
levels(decay_rate$species)<-c("ACERU","ACESA","BETPA","PINBA",
                              "PINST","POPTR","QUEMA","QUERU")
decay_rate$full_id<-apply(decay_rate[,c(1,2,4,3)])
