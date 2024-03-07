setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(reshape2)

###################################
## VIP plotting

VIP_intact<-readRDS("SavedResults/VIP_intact.rds")
VIP_ground<-readRDS("SavedResults/VIP_ground.rds")

focal_palette=palette(brewer.pal(8,name="Set2"))

VIP_intact_long<-melt(VIP_intact,id.vars = "wavelength")
VIP_ground_long<-melt(VIP_ground,id.vars = "wavelength")

levels(VIP_intact_long$variable)<-c("sol","hemi","recalc","LMA",
                                    "Cmass","Nmass","Carea","Narea")
levels(VIP_ground_long$variable)<-c("sol","hemi","recalc","LMA",
                                    "Cmass","Nmass")

plot_left<-c("sol","hemi","recalc","LMA")

VIP_intact_long$side<-ifelse(VIP_intact_long$variable %in% plot_left,
                             "left","right")
VIP_ground_long$side<-ifelse(VIP_ground_long$variable %in% plot_left,
                             "left","right")

VIP_intact_plot_facet<-ggplot(VIP_intact_long,
                         aes(x=wavelength,y=value,color=variable))+
  geom_line(linewidth=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Intact")+
  facet_wrap(~side)+
  scale_color_manual(values=focal_palette,
                     labels=c("solubles","hemicellulose",
                              "recalcitrants","LMA",
                              expression(C[mass]),
                              expression(N[mass]),
                              expression(C[area]),
                              expression(N[area])))+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ylim(c(0,3))

VIP_ground_plot_facet<-ggplot(VIP_ground_long,
                              aes(x=wavelength,y=value,color=variable))+
  geom_line(linewidth=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Ground")+
  facet_wrap(~side)+
  scale_color_manual(values=focal_palette,
                     labels=c("solubles","hemicellulose",
                              "recalcitrants","LMA",
                              expression(C[mass]),
                              expression(N[mass]),
                              expression(C[area]),
                              expression(N[area])))+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ylim(c(0,3))+guides(color="none")

pdf("Manuscript/Fig5.pdf",height=8,width=10)
(VIP_intact_plot_facet)/
  (VIP_ground_plot_facet)+
  plot_layout(guides="collect") & 
  theme(legend.position = "bottom")
dev.off()

###################################
## read in data

intact_jack_coef_list<-readRDS("SavedResults/intact_jack_coefs_list.rds")
intact_jack_df_list<-readRDS("SavedResults/intact_jack_df_list.rds")

ground_jack_coef_list<-readRDS("SavedResults/ground_jack_coefs_list.rds")
ground_jack_df_list<-readRDS("SavedResults/ground_jack_df_list.rds")

####################################
## define plot boundaries

all.sol<-c(intact_jack_df_list$sol$Measured,
           intact_jack_df_list$sol$pred_mean,
           ground_jack_df_list$sol$pred_mean)
sol_upper<-max(all.sol,na.rm=T)+3
sol_lower<-min(all.sol,na.rm=T)-3

all.hemi<-c(intact_jack_df_list$hemi$Measured,
            intact_jack_df_list$hemi$pred_mean,
            ground_jack_df_list$hemi$pred_mean)
hemi_upper<-max(all.hemi,na.rm=T)+1
hemi_lower<-min(all.hemi,na.rm=T)-1

all.recalc<-c(intact_jack_df_list$recalc$Measured,
              intact_jack_df_list$recalc$pred_mean,
              ground_jack_df_list$recalc$pred_mean)
recalc_upper<-max(all.recalc,na.rm=T)+2
recalc_lower<-min(all.recalc,na.rm=T)-2

all.Nmass<-c(intact_jack_df_list$Nmass$Measured,
             intact_jack_df_list$Nmass$pred_mean,
             ground_jack_df_list$Nmass$pred_mean)
Nmass_upper<-max(all.Nmass,na.rm=T)+0.2
Nmass_lower<-min(all.Nmass,na.rm=T)-0.2

all.Cmass<-c(intact_jack_df_list$Cmass$Measured,
             intact_jack_df_list$Cmass$pred_mean,
             ground_jack_df_list$Cmass$pred_mean)
Cmass_upper<-max(all.Cmass,na.rm=T)+2
Cmass_lower<-min(all.Cmass,na.rm=T)-2

all.Carea<-c(intact_jack_df_list$Carea$Measured,
             intact_jack_df_list$Carea$pred_mean,
             ground_jack_df_list$Carea$pred_mean)
Carea_upper<-max(all.Carea,na.rm=T)+10
Carea_lower<-min(all.Carea,na.rm=T)-10

all.Narea<-c(intact_jack_df_list$Narea$Measured,
             intact_jack_df_list$Narea$pred_mean,
             ground_jack_df_list$Narea$pred_mean)
Narea_upper<-max(all.Narea,na.rm=T)+0.2
Narea_lower<-min(all.Narea,na.rm=T)-0.2

all.LMA<-c(intact_jack_df_list$LMA$Measured,
           intact_jack_df_list$LMA$pred_mean,
           ground_jack_df_list$LMA$pred_mean)
LMA_upper<-max(all.LMA,na.rm=T)+25
LMA_lower<-min(all.LMA,na.rm=T)-25

####################################
## create plots

sol_intact_val_plot<-ggplot(intact_jack_df_list$sol,
                            aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(sol_lower,sol_upper),ylim=c(sol_lower,sol_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  ggtitle("Intact")

hemi_intact_val_plot<-ggplot(intact_jack_df_list$hemi,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemi_lower,hemi_upper),ylim=c(hemi_lower,hemi_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

recalc_intact_val_plot<-ggplot(intact_jack_df_list$recalc,
                               aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(recalc_lower,recalc_upper),c(recalc_lower,recalc_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured recalcalcitrants (%)",x="Predicted recalcitrants (%)")+
  guides(color=F)

Cmass_intact_val_plot<-ggplot(intact_jack_df_list$Cmass,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cmass_lower,Cmass_upper),ylim=c(Cmass_lower,Cmass_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Intact")

Nmass_intact_val_plot<-ggplot(intact_jack_df_list$Nmass,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Nmass_lower,Nmass_upper),ylim=c(Nmass_lower,Nmass_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

Carea_intact_val_plot<-ggplot(intact_jack_df_list$Carea,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(c(Carea_lower,Carea_upper),
                  c(Carea_lower,Carea_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[area]*" (g/m"^2*")"),
       x=expression("Predicted C"[area]*" (g/m"^2*")"))

Narea_intact_val_plot<-ggplot(intact_jack_df_list$Narea,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(c(Narea_lower,Narea_upper),
                  c(Narea_lower,Narea_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[area]*" (g/m"^2*")"),
       x=expression("Predicted N"[area]*" (g/m"^2*")"))+
  guides(color=F)

LMA_intact_val_plot<-ggplot(intact_jack_df_list$LMA,
                            aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured LMA (g/m"^2*")"),
       x=expression("Predicted LMA (g/m"^2*")"))+
  guides(color=F)

sol_ground_val_plot<-ggplot(ground_jack_df_list$sol,
                            aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(sol_lower,sol_upper),ylim=c(sol_lower,sol_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  ggtitle("Ground")

hemi_ground_val_plot<-ggplot(ground_jack_df_list$hemi,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemi_lower,hemi_upper),ylim=c(hemi_lower,hemi_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

recalc_ground_val_plot<-ggplot(ground_jack_df_list$recalc,
                               aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(recalc_lower,recalc_upper),ylim=c(recalc_lower,recalc_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured recalcitrants (%)",x="Predicted recalcitrants (%)")

Cmass_ground_val_plot<-ggplot(ground_jack_df_list$Cmass,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cmass_lower,Cmass_upper),ylim=c(Cmass_lower,Cmass_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Ground")

Nmass_ground_val_plot<-ggplot(ground_jack_df_list$Nmass,
                              aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Nmass_lower,Nmass_upper),ylim=c(Nmass_lower,Nmass_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

LMA_ground_val_plot<-ggplot(ground_jack_df_list$LMA,
                            aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (g/m"^2*")"),
       x=expression("Predicted LMA (g/m"^2*")"))

################################
## arrange plots

pdf("Manuscript/Fig2.pdf",height=13,width=10)
(sol_intact_val_plot/hemi_intact_val_plot/recalc_intact_val_plot)|
  (sol_ground_val_plot/hemi_ground_val_plot/recalc_ground_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

pdf("Manuscript/Fig3.pdf",height=13,width=10)
(Cmass_intact_val_plot/Nmass_intact_val_plot/LMA_intact_val_plot)|
  (Cmass_ground_val_plot/Nmass_ground_val_plot/LMA_ground_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

pdf("Manuscript/Fig4.pdf",height=8,width=6)
(Carea_intact_val_plot/Narea_intact_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

###########################################
## output species mean measured and predicted traits

my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "Species")
}

intact_pred_summ<-lapply(intact_jack_df_list,function(trait_df){
  aggregate(.~Species,FUN=function(x) signif(mean(x),digits = 3),
            data=trait_df[,c("pred_mean","Measured","Species")])
})

intact_pred_summ_df<-Reduce(my_merge,intact_pred_summ)
colnames(intact_pred_summ_df)<-c("Species",
                                 paste(rep(names(intact_jack_df_list),each=2),
                                       rep(c("Pred","Meas"),times=8),sep="_"))

write.csv(intact_pred_summ_df,"SavedResults/intact_pred_summ_df.csv",row.names = F)
