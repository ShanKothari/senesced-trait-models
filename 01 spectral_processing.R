library(spectrolab)
library(caret)
library(patchwork)
library(vegan)
library(rdist)
setwd("C:/Users/querc/Dropbox/TraitModels2018/SenescencePaper/")

#######################################
## read intact spectra

source("Scripts/senesced-trait-models/read_spec.R")
intact_spec <- read.sed.multiple("Intact_Spectra_Raw/")
intact_spec <- intact_spec[,400:2400]
intact_spec_names<-names(intact_spec)
intact_spec_spl<-strsplit(x = intact_spec_names, split = "_")
sp.names<-c("ACNE","ACRU","BEPA","PIBA","PIRE","PIST",
            "QUAL","QUEL","QUMA","QURU","TIAM")
meta(intact_spec)$plot<-unlist(lapply(intact_spec_spl, function(x) x[[1]]))
meta(intact_spec)$ind<-unlist(lapply(intact_spec_spl, function(x) ifelse(x[[2]] %in% sp.names,x[[3]],x[[2]])))
meta(intact_spec)$sp<-unlist(lapply(intact_spec_spl, function(x) ifelse(x[[2]] %in% sp.names,x[[2]],x[[3]])))
meta(intact_spec)$num<-unlist(lapply(intact_spec_spl, function(x) x[[4]]))
meta(intact_spec)$ID<-apply(meta(intact_spec),1,function(x) paste(x,collapse = "_"))
## decimal rather than percent
intact_spec<-intact_spec/100

## sensor matching is imperfect but better than nothing
intact_spec_matched<-match_sensors(intact_spec, splice_at=c(983), fixed_sensor = 2, interpolate_wvl = 5)
## spline smoothing if desired
## here we'll only do it where the sensors overlap (970-1000 nm)
intact_spec_smoothed<-intact_spec_matched
intact_spec_smoothed[,970:1000]<-smooth_spline(intact_spec_smoothed[,970:1000])

intact_spec_agg<-aggregate(intact_spec_smoothed,by=meta(intact_spec_matched)$ID,
                           FUN=try_keep_txt(mean))

#######################################
## read ground spectra

ground_spec<-read_spectra(path = "Ground_Spectra_Raw","sed",exclude_if_matches = c("BAD","WR"))
ground_spec <- ground_spec[,400:2400]
ground_spec_names<-names(ground_spec)
ground_spec_spl<-strsplit(x = ground_spec_names, split = "_")
meta(ground_spec)$plot<-unlist(lapply(ground_spec_spl, function(x) x[[1]]))
meta(ground_spec)$ind<-unlist(lapply(ground_spec_spl, function(x) x[[2]]))
meta(ground_spec)$sp<-unlist(lapply(ground_spec_spl, function(x) x[[3]]))
meta(ground_spec)$num<-unlist(lapply(ground_spec_spl, function(x) x[[4]]))
meta(ground_spec)$ID<-apply(meta(ground_spec),1,function(x) paste(x,collapse = "_"))

ground_spec_agg<-aggregate(ground_spec,by=meta(ground_spec)$ID,
                           FUN=try_keep_txt(mean))

## compare which samples are in each
setdiff(names(ground_spec_agg),names(intact_spec_agg))
setdiff(names(intact_spec_agg),names(ground_spec_agg))

############################################
## calculate mean spectral distance

intact_spec_split<-split(x=data.frame(as.matrix(intact_spec_smoothed)),
                         f=meta(intact_spec_smoothed)$ID)

intact_mean_dist<-lapply(intact_spec_split,
                       function(x) mean(rdist(x,metric="angular")))
median(unlist(intact_mean_dist)*180/pi)

ground_spec_split<-split(x=data.frame(as.matrix(ground_spec)),
                         f=meta(ground_spec)$ID)

ground_mean_dist<-lapply(ground_spec_split,
                       function(x) mean(rdist(x,metric="angular")))
median(unlist(ground_mean_dist)*180/pi)

############################################
## attach trait data

## read fiber data
fiber<-read.csv("FiberAnalysis/FiberSummary.csv")
fiber<-fiber[-which(fiber$REDO %in% c("YES")),]
fiber$solubles<-100-fiber$NDF
fiber$hemicellulose<-fiber$NDF-fiber$ADF

meta(intact_spec_agg)$solubles<-fiber$solubles[match(meta(intact_spec_agg)$ID,fiber$SampleName)]
meta(intact_spec_agg)$hemicellulose<-fiber$hemicellulose[match(meta(intact_spec_agg)$ID,fiber$SampleName)]
meta(intact_spec_agg)$recalcitrant<-fiber$ADF[match(meta(intact_spec_agg)$ID,fiber$SampleName)]
meta(intact_spec_agg)$FiberRun<-fiber$Run[match(meta(intact_spec_agg)$ID,fiber$SampleName)]

meta(ground_spec_agg)$solubles<-fiber$solubles[match(meta(ground_spec_agg)$ID,fiber$SampleName)]
meta(ground_spec_agg)$hemicellulose<-fiber$hemicellulose[match(meta(ground_spec_agg)$ID,fiber$SampleName)]
meta(ground_spec_agg)$recalcitrant<-fiber$ADF[match(meta(ground_spec_agg)$ID,fiber$SampleName)]
meta(ground_spec_agg)$FiberRun<-fiber$Run[match(meta(ground_spec_agg)$ID,fiber$SampleName)]

## elemental data
EAdata<-read.csv("ElementalAnalysis/SummaryTables/CNSummary_all.csv")
EAdata<-EAdata[-which(EAdata$REDO %in% c("YES")),]
EAdata<-EAdata[-grep("duplicate",EAdata$SampleID),]
EAdata<-EAdata[-grep("apple",EAdata$SampleID),]

meta(intact_spec_agg)$Cmass<-EAdata$C.percent[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]
meta(intact_spec_agg)$Nmass<-EAdata$N.percent[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]
meta(intact_spec_agg)$EARun<-EAdata$Run[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]

meta(ground_spec_agg)$Cmass<-EAdata$C.percent[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]
meta(ground_spec_agg)$Nmass<-EAdata$N.percent[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]
meta(ground_spec_agg)$EARun<-EAdata$Run[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]

## calculate LMA for broadleafs in g/cm2
LMAdata_broad<-read.csv("SLA/SLA_senesced_broadleaf_2_6_2020.csv")
LMAdata_broad<-LMAdata_broad[-which(LMAdata_broad$Bad=="YES"),]
LMAdata_broad$ID<-apply(LMAdata_broad[,c(1,3,2,4)],1,function(x) paste(x,collapse="_"))
LMAdata_broad$ID<-gsub(pattern = " ",replacement="",x = LMAdata_broad$ID)
LMAdata_broad$LMA<-with(LMAdata_broad,weight/(size.per.hole.punch..cm.2.*number.of.discs))
LMAdata_broad<-LMAdata_broad[,c("Plot","Species","Location","Number","LMA","ID")]

## calculate area for needleleafs in cm2
LMAdata_needle<-read.csv("SLA/SLA_senesced_needleleaf_3_17_2020.csv")
LMAdata_needle$area<-LMAdata_needle$LengthMm/10*LMAdata_needle$WidthMm/10
LMAdata_needle$LengthMm<-NULL
LMAdata_needle$WidthMm<-NULL
LMAdata_needle$SLA<-NULL
LMAdata_needle$NeedleNumber<-NULL

## aggregate masses and areas for a sample
## then calculate LMA in g/cm2
LMAdata_needle_agg<-aggregate(.~Plot*Species*Location*IndNumber,data=LMAdata_needle,FUN=sum)
LMAdata_needle_agg$LMA<-LMAdata_needle_agg$MassG/LMAdata_needle_agg$area
LMAdata_needle_agg$MassG<-NULL
LMAdata_needle_agg$area<-NULL
LMAdata_needle_agg$ID<-apply(LMAdata_needle_agg[,1:4],1,function(x) paste(x,collapse="_"))
colnames(LMAdata_needle_agg)<-c("Plot","Species","Location","Number","LMA","ID")

LMAdata<-rbind(LMAdata_broad,LMAdata_needle_agg)

## converting g/cm2 to g/m2
meta(intact_spec_agg)$LMA<-LMAdata$LMA[match(meta(intact_spec_agg)$ID,LMAdata$ID)]*10000
meta(intact_spec_agg)$Carea<-with(meta(intact_spec_agg),Cmass*LMA)/100
meta(intact_spec_agg)$Nmass_area<-with(meta(intact_spec_agg),Nmass*LMA)/100

meta(ground_spec_agg)$LMA<-LMAdata$LMA[match(meta(ground_spec_agg)$ID,LMAdata$ID)]*10000
meta(ground_spec_agg)$Carea<-with(meta(ground_spec_agg),Cmass*LMA)/100
meta(ground_spec_agg)$Nmass_area<-with(meta(ground_spec_agg),Nmass*LMA)/100

####################################
## divide training and testing

train_sample <- createDataPartition(
  y = meta(ground_spec_agg)$sp,
  p = .75,
  list = FALSE
)
test_sample<-setdiff(1:nrow(ground_spec_agg),train_sample)

ground_spec_agg_train<-ground_spec_agg[train_sample,]
ground_spec_agg_test<-ground_spec_agg[test_sample,]

## this assigns the one intact specimen not in the ground dataset
## to the training data
intact_spec_agg_train<-intact_spec_agg[-which(names(intact_spec_agg) %in% names(ground_spec_agg_test)),]
intact_spec_agg_test<-intact_spec_agg[which(names(intact_spec_agg) %in% names(ground_spec_agg_test)),]

saveRDS(intact_spec_agg,"SavedResults/intact_spec_agg.rds")
saveRDS(intact_spec_agg_train,"SavedResults/intact_spec_agg_train.rds")
saveRDS(intact_spec_agg_test,"SavedResults/intact_spec_agg_test.rds")
write.csv(as.data.frame(intact_spec_agg)[,-1],"SavedResults/intact_spec.csv",row.names = F)

saveRDS(ground_spec_agg,"SavedResults/ground_spec_agg.rds")
saveRDS(ground_spec_agg_train,"SavedResults/ground_spec_agg_train.rds")
saveRDS(ground_spec_agg_test,"SavedResults/ground_spec_agg_test.rds")
write.csv(as.data.frame(ground_spec_agg)[,-1],"SavedResults/ground_spec.csv",row.names = F)

####################################
## plot spectra

intact_quantiles<-quantile(intact_spec_agg,probs=c(0.025,0.25,0.5,0.75,0.975))
intact_CV<-apply(intact_spec_agg,2,function(x) sd(x)/mean(x))

intact_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(intact_quantiles)[1,],
                  ymax = as.matrix(intact_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(intact_quantiles)[2,],
                  ymax = as.matrix(intact_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(intact_quantiles)[3,]),
            linewidth=1,color="black")+
  geom_line(aes(x=400:2400,y=intact_CV),
            linewidth=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x="Wavelength",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ggtitle("Intact-leaf spectra")

ground_quantiles<-quantile(ground_spec_agg,probs=c(0.025,0.25,0.5,0.75,0.975))
ground_CV<-apply(ground_spec_agg,2,function(x) sd(x)/mean(x))

ground_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ground_quantiles)[1,],
                  ymax = as.matrix(ground_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ground_quantiles)[2,],
                  ymax = as.matrix(ground_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(ground_quantiles)[3,]),
            linewidth=1,color="black")+
  geom_line(aes(x=400:2400,y=ground_CV),
            linewidth=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="Wavelength",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  ggtitle("Ground-leaf spectra")

pdf("Manuscript/Fig1.pdf",height=8,width=6)
intact_spec_plot/ground_spec_plot
dev.off()
