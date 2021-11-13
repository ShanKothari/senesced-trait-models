read.sed.multiple<-function(path){
  spec.files<-list.files(path=path,pattern=".sed")
  spec.files<-spec.files[-grep(pattern="WR|BAD|TEST|CHECK",spec.files)]
  
  read.sed<-function(path,filename){
    lines<-readLines(paste(path,filename,sep=""))
    lines<-lines[-(1:26)]
    lines.split<-strsplit(lines,split="\t")
    refl.col<-which(lines.split[[1]]=="Reflect. %")
    reflectance.char<-lapply(lines.split,function(x) x[[refl.col]])
    refl.sp<-as.numeric(gsub(pattern = " ", replacement = "", reflectance.char[-1]))
    return(refl.sp)
  }

  spec.list<-lapply(spec.files,function(x) read.sed(path,x))
  spec.matrix<-t(matrix(unlist(spec.list),
                        ncol=length(spec.list),
                        nrow=length(spec.list[[1]])))
  colnames(spec.matrix)<-350:2500
  spec.all<-spectra(spec.matrix,bands=350:2500,names=spec.files)
  
  return(spec.all)
}
