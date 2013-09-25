find.mode<-function(x)
{
  d<-density(x)
  return(d$x[which(d$y==max(d$y))[1]])
}


split.channels.file<-function(file)
{
test<-try({
img<-readTIF(paste("rgb",file,sep="/"))

img<-img-min(img)
img<-img/(max(img))

D<-length(dim(img))
if (D==4)
{
  red<-img[,,1,]
  green<-img[,,2,]
  blue<-img[,,3,]
}
if (D==3)
{
  Z<-dim(img)[3]
  red<-img[,,seq(1,Z,by=3)]  
  green<-img[,,seq(2,Z,by=3)]  
  blue<-img[,,seq(3,Z,by=3)]  
}

red<-red-find.mode(red)
green<-green-find.mode(green)
blue<-blue-find.mode(blue)

red[red<0]<-0
green[green<0]<-0
blue[blue<0]<-0

red<-red-min(red)
green<-green-min(green)
blue<-blue-min(blue)
red<-red/max(red)
green<-green/max(green)
blue<-blue/max(blue)


writeTIF(blue,paste("blue/",file,sep=""),bps=16L)
writeTIF(green,paste("green/",file,sep=""),bps=16L)
writeTIF(red,paste("red/",file,sep=""),bps=16L)

Xmic<-attr(img,"x.resolution")
Ymic<-attr(img,"y.resolution")
Zmic<-as.numeric(attr(img,"slices"))*as.numeric(attr(img,"spacing"))
write(c(Xmic,Ymic,Zmic),file=paste("XYZmic/",file,".txt",sep=""))

remove(img,red,green,blue)
gc(verbose=FALSE)
},silent=TRUE)
if(class(test)=="try-error")cat(paste0(file,": ",attr(test,"condition"),"\n"))
else(cat(paste0(file," OK\n")))
}

split.channels<-function(f,cores=1)
{
  orig<-getwd()
  setwd(f)
  require(bioimagetools)
  if(cores>1)
  {
    require(parallel)
    options("mc.cores"=cores)
  }
  
  files<-list.files("rgb")
  cat(paste(length(files),"files.\n"))
                
  if(length(list.files("red"))==0)dir.create("red")
  if(length(list.files("blue"))==0)dir.create("blue")
  if(length(list.files("green"))==0)dir.create("green")
  if(length(list.files("XYZmic"))==0)dir.create("XYZmic")
  
  if(cores>1)jobs <- mclapply(files,split.channels.file)
  if(cores==1)jobs <- lapply(files,split.channels.file)
  setwd(orig)
}
