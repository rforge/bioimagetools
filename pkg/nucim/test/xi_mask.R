library(bioimagetools)

setwd("/home/schmid/projects/marion/C2C12")
files<-sample(list.files("green"))
file<-files[1]
    print(file)
    mask<-readTIF(paste("dapimask/",file,sep=""))
    blau<-readTIF(paste("blue/",file,sep=""))
    prot<-readTIF(paste("green/",file,sep=""))
    pm<-median(prot)
    prot[mask==0]<-pm
 #   prot<-prot*2^16
#    storage.mode(prot)<-"integer"
    
    brush<-makeBrush(25,shape="gaussian",sigma=5)
    prot1<-filterImage2d(prot,brush)
    
    prot1<-prot1-min(prot1)
    prot1<-prot1/max(prot1)
    
    prot13<-prot1>.385
    
    prot4<-bwlabel3d(prot13)
    prot5<-cmoments3d(prot4,prot/(2^16))
    
    nr.xi<-2
    if(length(grep("mES",file))==0)nr.xi<-2
    if(dim(prot5)[1]==1)nr.xi<-1
    
    which<-rev(order(prot5[,5]))[1:nr.xi]
    xi<-prot5[which,2:4]
    
    xi.mask<-list()
    for (i in 1:nr.xi)
      xi.mask[[i]]<-(prot4==which[i])
    
    #   blue<-readImage(paste("rgb/",substr(file,0,18),"blue.tif",sep=""))
    
    #    test<-rgbImage(blue=blue,green=prot,red=xi.mask[[1]])
    #    writeImage(test,file=paste("xistbereich/1_",file,sep=""))
    #    if(nr.xi==2)
    #    {
    #      test<-rgbImage(blue=blue,green=prot,red=xi.mask[[2]])
    #     writeImage(test,file=paste("test_",file,sep=""))
    #    }
    
    save(xi.mask,file=paste("xist/",substr(file,0,17),".Rdata",sep=""))
  })
}
