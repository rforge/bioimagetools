require(bioimagetools)
require(EBImage)

test<-readImage("/home/schmid/bioimg/software/bioimagetools/test/b_090623_25_DAPI.tif")

test2<-test[,,35]
image(test2)

seg<-segment(test2,nclust=5,beta=.2)
image(seg$class)

seg3d<-segment(test,nclust=5,beta=.2)
image(seg3d$class[,,35])

mask<-maskdapi(img=test,thresh=.3)

std<-standardize(test,mask=mask,N=32,sd=6)
table0(std)
colors.in.classes(std,std)
  