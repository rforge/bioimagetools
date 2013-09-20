classify<-function(f,N,cores=1)
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
if(length(list.files(paste0("class",N)))==0)dir.create(paste0("class",N))

if(cores>1)jobs <- mclapply(files, classify.file, N=N, mc.preschedule=FALSE)
if(cores==1)jobs <- lapply(files, classify.file, N=N)
setwd(orig)
}

classify.file<-function(file,N)
{
  mask<-try(readTIF(paste("dapimask/",file,sep=""),native=TRUE))
  print(summary(mask))
  blau<-readTIF(paste("blue/",file,sep=""))
  blau<-array(blau,dim(blau))
  blau<-round(blau*2^16)
  storage.mode(blau)<-"integer"
  test<-try({jobss
    img.seg<-segment(blau,N,0.1,1/3,mask=(mask==1),maxit=50,varfixed=TRUE,
                     inforce.nclust=TRUE, start="equal")
    classes<-array(as.integer(img.seg$class),dim(blau))
    remove(img.seg)
    writeTIF(classes/N,paste("class",N,"/",file,sep=""),bps=8)
  })
}
