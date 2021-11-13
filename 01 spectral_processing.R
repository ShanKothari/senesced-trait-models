library(spectrolab)
setwd("C:/Users/kotha020/Dropbox/TraitModels2018/SenescencePaper/")

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

## sensor matching is imperfect
intact_spec_matched<-match_sensors(intact_spec, splice_at=c(983), fixed_sensor = 2, interpolate_wvl = 5,
                                   factor_range = c(0.5, 3))
## spline smoothing if desired
intact_spec_matched<-smooth(intact_spec_matched)

intact_spec_agg<-aggregate(intact_spec_matched,by=meta(intact_spec_matched)$ID,
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

ground_spec_split<-split(x=data.frame(reflectance(ground_spec)),f=meta(ground_spec)$ID)

mean_spec_dist<-lapply(ground_spec_split,
                       function(x) mean(rdist(x,metric="angular")))
mean(unlist(mean_spec_dist)*180/pi)

############################################
## attach trait data

## read fiber data
fiber<-read.csv("FiberAnalysis/FiberSummary.csv")
fiber<-fiber[-which(fiber$REDO %in% c("YES")),]

meta(intact_spec_agg)$NDF<-fiber$NDF[match(meta(intact_spec_agg)$ID,fiber$SampleName)]/100
meta(intact_spec_agg)$ADF<-fiber$ADF[match(meta(intact_spec_agg)$ID,fiber$SampleName)]/100
meta(intact_spec_agg)$FiberRun<-fiber$Run[match(meta(intact_spec_agg)$ID,fiber$SampleName)]

meta(ground_spec_agg)$NDF<-fiber$NDF[match(meta(ground_spec_agg)$ID,fiber$SampleName)]/100
meta(ground_spec_agg)$ADF<-fiber$ADF[match(meta(ground_spec_agg)$ID,fiber$SampleName)]/100
meta(ground_spec_agg)$FiberRun<-fiber$Run[match(meta(ground_spec_agg)$ID,fiber$SampleName)]

## elemental data
EAdata<-read.csv("ElementalAnalysis/SummaryTables/CNSummary_all.csv")
EAdata<-EAdata[-which(EAdata$REDO %in% c("YES")),]
EAdata<-EAdata[-grep("duplicate",EAdata$SampleID),]
EAdata<-EAdata[-grep("apple",EAdata$SampleID),]

meta(intact_spec_agg)$perC<-EAdata$C.percent[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]/100
meta(intact_spec_agg)$perN<-EAdata$N.percent[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]/100
meta(intact_spec_agg)$EARun<-EAdata$Run[match(meta(intact_spec_agg)$ID,EAdata$SampleID)]

meta(ground_spec_agg)$perC<-EAdata$C.percent[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]/100
meta(ground_spec_agg)$perN<-EAdata$N.percent[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]/100
meta(ground_spec_agg)$EARun<-EAdata$Run[match(meta(ground_spec_agg)$ID,EAdata$SampleID)]

## LMA data
LMAdata_broad<-read.csv("SLA/SLA_senesced_broadleaf_2_6_2020.csv")
LMAdata_broad<-LMAdata_broad[-which(LMAdata_broad$Bad=="YES"),]
LMAdata_broad$ID<-apply(LMAdata_broad[,c(1,3,2,4)],1,function(x) paste(x,collapse="_"))
LMAdata_broad$ID<-gsub(pattern = " ",replacement="",x = LMAdata_broad$ID)
LMAdata_broad$LMA<-with(LMAdata_broad,weight/(size.per.hole.punch..cm.2.*number.of.discs))
LMAdata_broad<-LMAdata_broad[,c("Plot","Species","Location","Number","LMA","ID")]

## calculate area
LMAdata_needle<-read.csv("SLA/SLA_senesced_needleleaf_3_17_2020.csv")
LMAdata_needle$area<-LMAdata_needle$LengthMm/10*LMAdata_needle$WidthMm/10
LMAdata_needle$LengthMm<-NULL
LMAdata_needle$WidthMm<-NULL
LMAdata_needle$SLA<-NULL
LMAdata_needle$NeedleNumber<-NULL

## aggregate and calculate LMA
LMAdata_needle_agg<-aggregate(.~Plot*Species*Location*IndNumber,data=LMAdata_needle,FUN=sum)
LMAdata_needle_agg$LMA<-LMAdata_needle_agg$MassG/LMAdata_needle_agg$area
LMAdata_needle_agg$MassG<-NULL
LMAdata_needle_agg$area<-NULL
LMAdata_needle_agg$ID<-apply(LMAdata_needle_agg[,1:4],1,function(x) paste(x,collapse="_"))
colnames(LMAdata_needle_agg)<-c("Plot","Species","Location","Number","LMA","ID")

LMAdata<-rbind(LMAdata_broad,LMAdata_needle_agg)

meta(intact_spec_agg)$LMA<-LMAdata$LMA[match(meta(intact_spec_agg)$ID,LMAdata$ID)]
meta(intact_spec_agg)$perC_area<-with(meta(intact_spec_agg),perC*LMA)
meta(intact_spec_agg)$perN_area<-with(meta(intact_spec_agg),perN*LMA)

meta(ground_spec_agg)$LMA<-LMAdata$LMA[match(meta(ground_spec_agg)$ID,LMAdata$ID)]
meta(ground_spec_agg)$perC_area<-with(meta(ground_spec_agg),perC*LMA)
meta(ground_spec_agg)$perN_area<-with(meta(ground_spec_agg),perN*LMA)

####################################
## divide training and testing

train_sample <- createDataPartition(
  y = meta(ground_spec_agg)$sp,
  p = .8,
  list = FALSE
)
test_sample<-setdiff(1:nrow(reflectance(ground_spec_agg)),train_sample)

ground_spec_agg_train<-ground_spec_agg[train_sample,]
ground_spec_agg_test<-ground_spec_agg[test_sample,]

intact_spec_agg_train<-intact_spec_agg[train_sample,]
intact_spec_agg_test<-intact_spec_agg[test_sample,]

####################################
## plot spectra

intact_quantiles<-quantile(intact_spec_agg,probs=c(0,0.025,0.5,0.975,1))
intact_CV<-apply(intact_spec_agg,2,function(x) sd(x)/mean(x))

intact_spec_plot<-ggplot()+
  geom_line(aes(x=400:2400,y=reflectance(intact_quantiles)[3,]),size=1)+
  geom_line(aes(x=400:2400,y=reflectance(intact_quantiles)[2,]),size=1,linetype="dashed")+
  geom_line(aes(x=400:2400,y=reflectance(intact_quantiles)[4,]),size=1,linetype="dashed")+
  geom_line(aes(x=400:2400,y=reflectance(intact_quantiles)[1,]),size=1,linetype="dotted")+
  geom_line(aes(x=400:2400,y=reflectance(intact_quantiles)[5,]),size=1,linetype="dotted")+
  geom_line(aes(x=400:2400,y=intact_CV),size=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="Wavelength",y="Reflectance (or CV)")+lims(y=c(0,1))+
  ggtitle("Intact-leaf spectra")

ground_quantiles<-quantile(ground_spec_agg,probs=c(0,0.025,0.5,0.975,1))
ground_CV<-apply(ground_spec_agg,2,function(x) sd(x)/mean(x))

ground_spec_plot<-ggplot()+
  geom_line(aes(x=400:2400,y=reflectance(ground_quantiles)[3,]),size=1)+
  geom_line(aes(x=400:2400,y=reflectance(ground_quantiles)[2,]),size=1,linetype="dashed")+
  geom_line(aes(x=400:2400,y=reflectance(ground_quantiles)[4,]),size=1,linetype="dashed")+
  geom_line(aes(x=400:2400,y=reflectance(ground_quantiles)[1,]),size=1,linetype="dotted")+
  geom_line(aes(x=400:2400,y=reflectance(ground_quantiles)[5,]),size=1,linetype="dotted")+
  geom_line(aes(x=400:2400,y=ground_CV),size=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(x="Wavelength",y="Reflectance (or CV)")+lims(y=c(0,1))+
  ggtitle("Ground-leaf spectra")