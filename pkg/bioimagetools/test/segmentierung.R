require(bioimagetools)
require(EBImage)

test<-readImage("/home/schmid/bioimg/software/bioimagetools/test/b_090623_25_DAPI.tif")

test2<-test[,,35]
image(test2)

Z<-which(apply(test,3,mean)>.04)
mask<-maskdapi(img=test,thresh=.2,filter=3)
image(mask[,,35])


seg1<-segment(test,nclust=5,beta=.2,mask=mask)
image(seg1$class[,,35])

seg1<-segment(test,nclust=5,beta=.2,mask=mask,varfixed=FALSE)
image(seg1$class[,,35])

seg1<-segment(test2,nclust=5,beta=.2,mask=mask[,,35],inforce.nclust=TRUE)
image(seg1$class)

seg1<-segment(test2,nclust=5,beta=.2,mask=mask[,,35],varfixed=FALSE)
image(seg1$class)

seg1<-segment(test2,nclust=5,beta=.2,mask=mask[,,35],varfixed=FALSE)
image(seg1$class[,,35])


std<-standardize(test,mask=mask,N=32,sd=6)
table.n(std,32)
colors.in.classes(std,std)
  
test2<-runif(128*128,0,1)
test2<-sort(test2)
test2<-array(test2,c(128,128))
green<-array(test2*runif(128*128,0,1),c(128,128))
blue<-array(runif(128*128,0,1)^2,c(128,128))
image(blue)

std<-standardize(test2,N=32,sd=4)
t<-table.n(std,32)
barplot(t)
colors.in.classes(std,green)
colors.in.classes(std,blue,col1="blue")
