XYZ.mic<-function(f,cores=1)
{
  orig<-getwd()
setwd(f)

library(bioimagetools)
if (cores>1)
{
  library(parallel)
  options("mc.cores"=cores)
}

    files<-sample(list.files("blue"))
    dir<-"XYZmic"
    if(length(list.files(dir))==0)dir.create(dir)
    
    for (file in files)
    {
 try({
      print(file)
      img<-readTIF(paste("rgb",file,sep="/"))
      Xmic<-attr(img,"x.resolution")
      Ymic<-attr(img,"y.resolution")
      Zmic<-as.numeric(attr(img,"slices"))*as.numeric(attr(img,"spacing"))
      remove(img)
      write(c(Xmic,Ymic,Zmic),file=paste(dir,"/",file,".txt",sep=""))
   })
}
  setwd(orig)
}

