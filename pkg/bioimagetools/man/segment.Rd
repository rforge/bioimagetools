\name{segment}
\alias{segment}
\title{Segmentation of 3D images using EM algorithm}
\description{
Segmentation of 3D images using EM algorithm
}
\usage{
segment(img, nclust, beta, z.scale=0,varfixed=TRUE,maxit=30, mask=array(TRUE,dim(img)), priormu=rep(NA,length(nclust)), priormusd=rep(NULL,length(nclust)) ) }
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
}
\value{
 A list with "class": 3d array of class per voxel; "mu" estimated means; "sigma": estimated standard deviations. 
}
\author{Volker J. Schmid \email{volkerschmid@users.sourceforge.net}}


