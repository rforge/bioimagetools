dapimask<-function(f,cores=1)
  {
orig<-getwd()
setwd(f)

  library(bioimagetools)
  if (cores>1)
    {
    library(parallel)
    options("mc.cores"=cores)
  }
files<-(list.files("blue"))
if(length(list.files("dapimask"))==0)dir.create("dapimask")


if(cores>1)jobs <- mclapply(files, dapimask.file, mc.preschedule=FALSE)
if(cores==1)jobs <- lapply(files, dapimask.file)
setwd(orig)
}

dapimask.file<-function(file){
  print(file)
  blau<-readTIF(paste("blue/",file,sep=""))
  mb<-apply(blau,3,mean)
  mbr<-0.3*sum(range(mb))
  mbr<-which(mbr<mb)
  small<-min(mbr):max(mbr)
  dims0<-dim(blau)
  blau<-blau[,,small]
  dims<-dim(blau)
  blau<-blau-median(blau)
  blau[blau<0]<-0
  blau<-array(blau,dims)
  blau<-blau/max(blau)
  blau<-filter(blau,"var",4,1/3)
  b<-blau>mean(blau)
  #b<-blau>quantile(blau,.8)
  b2<-array(0,dims0)
  b2[,,small]<-array(as.integer(b),dim(b))
  n<-5
  mask<-1-outside(b2,0,n)
  brush<-makeBrush(2*n-1,shape='box')
  mask<-erode(mask,brush)
  
  if(0)
  { 
    z<-25
    bb<-mask[,,z]-.5
    bb<-bb*blau[,,z]
    image(bb)
  }
  
  #mask0<-1-outside(b,0,15)
  #mask<-array(0,dims0)
  #mask[,,small]<-array(as.integer(mask0),dim(mask0))
  writeTIF(mask,paste("dapimask/",file,sep=""),bps=8)
}  
