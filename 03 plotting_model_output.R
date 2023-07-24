setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")
library(ggplot2)
library(patchwork)
library(RColorBrewer)

###################################
## VIP plotting

VIP_intact<-readRDS("SavedResults/VIP_intact.rds")
VIP_ground<-readRDS("SavedResults/VIP_ground.rds")

focal_palette=palette(brewer.pal(8,name="Set2")[c(1,3,4,5,6,8)])

VIP_intact_long<-melt(VIP_intact,id.vars = "wavelength")
VIP_ground_long<-melt(VIP_ground,id.vars = "wavelength")

levels(VIP_intact_long$variable)<-c("sol","hemi","recalc","Cmass","Nmass","LMA")
levels(VIP_ground_long$variable)<-c("sol","hemi","recalc","Cmass","Nmass","LMA")

VIP_intact_plot<-ggplot(VIP_intact_long,
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Intact")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ylim(c(0,3))+guides(color=F)

VIP_ground_plot<-ggplot(VIP_ground_long,
                        aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1.25)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="VIP",x="Wavelength (nm)",color = "Trait")+
  ggtitle("Ground")+
  scale_color_manual(values=focal_palette)+
  geom_hline(yintercept=0.8,linetype="dashed",size=2)+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ylim(c(0,3))

pdf("Manuscript/Fig5.pdf",height=8,width=7)
VIP_intact_plot/VIP_ground_plot+
  plot_layout(guides="collect") & theme(legend.position = "right")
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

all.perN<-c(intact_jack_df_list$perN$Measured,
            intact_jack_df_list$perN$pred_mean,
            ground_jack_df_list$perN$pred_mean)
perN_upper<-max(all.perN,na.rm=T)+0.2
perN_lower<-min(all.perN,na.rm=T)-0.2

all.perC<-c(intact_jack_df_list$perC$Measured,
            intact_jack_df_list$perC$pred_mean,
            ground_jack_df_list$perC$pred_mean)
perC_upper<-max(all.perC,na.rm=T)+2
perC_lower<-min(all.perC,na.rm=T)-2

all.perC_area<-c(intact_jack_df_list$perC_area$Measured,
                 intact_jack_df_list$perC_area$pred_mean,
                 ground_jack_df_list$perC_area$pred_mean)
perC_area_upper<-max(all.perC_area,na.rm=T)+10
perC_area_lower<-min(all.perC_area,na.rm=T)-10

all.perN_area<-c(intact_jack_df_list$perN_area$Measured,
                 intact_jack_df_list$perN_area$pred_mean,
                 ground_jack_df_list$perN_area$pred_mean)
perN_area_upper<-max(all.perN_area,na.rm=T)+0.2
perN_area_lower<-min(all.perN_area,na.rm=T)-0.2

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

perC_intact_val_plot<-ggplot(intact_jack_df_list$perC,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Intact")

perN_intact_val_plot<-ggplot(intact_jack_df_list$perN,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

perC_area_intact_val_plot<-ggplot(intact_jack_df_list$perC_area,
                                  aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(c(perC_area_lower,perC_area_upper),
                  c(perC_area_lower,perC_area_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[area]*" (g/m"^2*")"),
       x=expression("Predicted C"[area]*" (g/m"^2*")"))+
  ggtitle("Intact")

perN_area_intact_val_plot<-ggplot(intact_jack_df_list$perN_area,
                                  aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(c(perN_area_lower,perN_area_upper),
                  c(perN_area_lower,perN_area_upper))+
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

perC_ground_val_plot<-ggplot(ground_jack_df_list$perC,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Ground")

perN_ground_val_plot<-ggplot(ground_jack_df_list$perN,
                             aes(y=Measured,x=pred_mean,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred_low,xmax=pred_high),
                 color="darkslategrey",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
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
(perC_intact_val_plot/perN_intact_val_plot/LMA_intact_val_plot)|
  (perC_ground_val_plot/perN_ground_val_plot/LMA_ground_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

pdf("Manuscript/Fig4.pdf",height=8,width=5.5)
(perC_area_intact_val_plot/perN_area_intact_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()
