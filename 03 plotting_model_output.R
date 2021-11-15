setwd("C:/Users/kotha020/Dropbox/TraitModels2018/SenescencePaper/")
library(ggplot2)
library(patchwork)
library(RColorBrewer)

###################################
## VIP plotting

VIP_intact<-readRDS("SavedResults/VIP_intact.rds")
VIP_ground<-readRDS("SavedResults/VIP_ground.rds")

focal_palette=palette(brewer.pal(8,name="Set2")[c(3,4,5,6,8)])

VIP_intact_long<-melt(VIP_intact,id.vars = "wavelength")
VIP_ground_long<-melt(VIP_ground,id.vars = "wavelength")

VIP_intact_plot<-ggplot(VIP_intact_long[VIP_intact_long$variable %in% c("NDF","ADF","perC","perN","LMA"),],
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
## create plots

NDF_intact_val_plot<-ggplot(intact_jack_df_list$NDF,
                            aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured NDF (%)",x="Predicted NDF (%)")+
  guides(color=F)+
  ggtitle("Intact")

ADF_intact_val_plot<-ggplot(intact_jack_df_list$ADF,
                            aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured ADF (%)",x="Predicted ADF (%)")+
  guides(color=F)

perC_intact_val_plot<-ggplot(intact_jack_df_list$perC,
                             aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Intact")

perN_intact_val_plot<-ggplot(intact_jack_df_list$perN,
                             aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

perC_area_intact_val_plot<-ggplot(intact_jack_df_list$perC_area,
                                  aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,200),ylim=c(0,200))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured C"[area]*" (g/m"^2*")"),
       x=expression("Predicted C"[area]*" (g/m"^2*")"))+
  ggtitle("Intact")

perN_area_intact_val_plot<-ggplot(intact_jack_df_list$perN_area,
                                  aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured N"[area]*" (g/m"^2*")"),
       x=expression("Predicted N"[area]*" (g/m"^2*")"))+
  guides(color=F)

LMA_intact_val_plot<-ggplot(intact_jack_df_list$LMA,
                            aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured LMA (g/m"^2*")"),
       x=expression("Predicted LMA (g/m"^2*")"))+
  guides(color=F)

NDF_ground_val_plot<-ggplot(ground_jack_df_list$NDF,
                            aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,75),ylim=c(25,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured NDF (%)",x="Predicted NDF (%)")+
  guides(color=F)+
  ggtitle("Ground")

ADF_ground_val_plot<-ggplot(ground_jack_df_list$ADF,
                            aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(18,58),ylim=c(18,58))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured ADF (%)",x="Predicted ADF (%)")+
  guides(color=F)

perC_ground_val_plot<-ggplot(ground_jack_df_list$perC,
                             aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(40,63),ylim=c(40,63))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)+
  ggtitle("Ground")

perN_ground_val_plot<-ggplot(ground_jack_df_list$perN,
                             aes(y=Measured*100,x=pred_mean*100,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*100,xmin=pred_low*100,xmax=pred_high*100),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2.5),ylim=c(0,2.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

LMA_ground_val_plot<-ggplot(ground_jack_df_list$LMA,
                            aes(y=Measured*10000,x=pred_mean*10000,color=Species))+
  theme_bw()+
  geom_errorbarh(aes(y=Measured*10000,xmin=pred_low*10000,xmax=pred_high*10000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,300),ylim=c(0,300))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (g/m"^2*")"),
       x=expression("Predicted LMA (g/m"^2*")"))+
  guides(color=F)

################################
## arrange plots

pdf("Manuscript/Fig2.pdf",height=12,width=8)
(perC_intact_val_plot/perN_intact_val_plot/LMA_intact_val_plot)|
  (perC_ground_val_plot/perN_ground_val_plot/LMA_ground_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/Fig3.pdf",height=8,width=8)
(NDF_intact_val_plot/ADF_intact_val_plot)|
  (NDF_ground_val_plot/ADF_ground_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Manuscript/Fig4.pdf",height=12,width=6)
(perC_area_intact_val_plot/perN_area_intact_val_plot) +
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()
