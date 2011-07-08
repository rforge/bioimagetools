\name{segment}
\alias{segment}
\title{Segmentation of 3D images using EM algorithm}
\description{
Segmentation of 3D images using EM algorithm
}
\usage{
segment(img, nclust, beta, z.scale=0,varfixed=TRUE,maxit=100, mask=array(TRUE,dim(img)), priormu=rep(NA,length(nclust)), priormusd=rep(NULL,length(nclust)), min.eps=10^{-7} ) }
\arguments{
  \item{img}{is a 3d arrary representing an image.}
  \item{nclust}{is the number of clusters/classes to be segmented.}
  \item{beta}{is a matrix of size nclust x nclust, representing the prior weight of classes neighbouring each other.}
  \item{z.scale}{ratio of voxel dimension in x/y direction and z direction. Will be multiplied on beta for neighbouring voxel in z direction.}
  \item{varfixed}{is a logical variable. If TRUE, variacne is equal in each class.}
  \item{maxit}{is the maximum number of iterations.}
  \item{mask}{is a logical array, representing the voxels to be used in the segmentation.}
  \item{priormu}{is a vector with mean of the normal prior of the expected values of all classes. Default is NA, which represents no prior assumption.}
  \item{priormusd}{is a vector with standard deviations of the normal prior of the expected values of all classes.}
  \item{min.eps}{stop criterion. Minimal change in sum of squared estimate of mean in order to stop.}
}
\value{
 A list with "class": 3d array of class per voxel; "mu" estimated means; "sigma": estimated standard deviations. 
}
\examples{
original<-array(1,c(300,300,50))
original[,,20:50]<-2
original[1:150,,]<-3
original[,200:300,]<-1

img<-array(rnorm(300*300*50,original,.2),c(300,300,50))
img<-img-min(img)
img<-img/max(img)

beta<-matrix(rep(-.25,9),nrow=3)
beta<-beta+0.75*diag(3)

seg.img<-segment(img,3,beta,z.scale=.3,varfixed=FALSE)

print(sum(seg.img$class!=original))
}
\author{Volker J. Schmid \email{volkerschmid@users.sourceforge.net}}


