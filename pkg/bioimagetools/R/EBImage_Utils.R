# Utilities for EBImage
# 
# Author: fabians
###############################################################################

#' Labels connected objects in a binary image stack.
#' @param im a stack of binary images (or a 3d-array)   
#' @return A Grayscale Image object or an array, containing the labelled version of im.
bwlabel3d <- function(im){
	res <- array(0, dim = dim(im))
	depth <- dim(im)[3]
	
	label <- bwlabel(im[,,1])
	#make sure labels are unique for each slice
	label[label != 0] <- label[label != 0] + 10^(1+floor(log(max(label), 10)))
	updatedLabels <- activeLabels <- unique(as.vector(label[label != 0])) 
	
	cluster <- 1:max(label)
	res[,,1] <- label

	cat("labelling.")
	
	for(i in 2:depth){
		cat(".")
		newlabel <- bwlabel(im[,,i])
		not0 <- newlabel != 0
		#make sure labels are unique for each slice
		newlabel[not0] <- newlabel[not0] + i*10^(1+floor(log(max(newlabel), 10)))
		done <- !not0
		alreadyWritten <- !not0
		
		if(length(activeLabels)){
			for(l in activeLabels){
				#propagate label 'l' into this slice in stack if there is overlap
				overlap <- not0 & label == l
				if(any(overlap)){
					#overwrite all pixels with the same labels as those in the overlap with label 'l'
					overwrite <- apply(newlabel, 2, "%in%", unique(as.vector(newlabel[overlap])))
			
					# check if pixels in overlap have already been relabeled,
					# if yes the label to be assigned can be replaced with the previously assigned one 
					if(any(alreadyWritten[overwrite]>0)){
						relabel <- which(res == l, arr.ind=T)
						previousLabel <-  unique(res[,,i][overwrite])
						#necessary if overwrite covers area of more than 2 different labels
						previousLabel <- previousLabel[previousLabel!=0] 
						res[relabel] <- previousLabel
						#remove label 'l' from activeLabels 
						updatedLabels <- updatedLabels[updatedLabels!=l]
					} else {
						res[,,i][overwrite] <- l
					}	
					done[overwrite] <- TRUE
					
					alreadyWritten <- alreadyWritten + overwrite
				} else {
					#remove label 'l' from activeLabels if there is no overlap
					updatedLabels <- updatedLabels[updatedLabels!=l]
				}
			}
			
		}
		if(any(!done)){
			#add new labels for objects
			updatedLabels <-c(updatedLabels, unique(as.vector(newlabel[!done])))
			res[,,i][!done] <- newlabel[!done]
		} 
		activeLabels <- updatedLabels
		label <- res[,,i] 
	}
	
	cat(", relabelling")

	#re-label with 1 to no. of objects
  #changed V.S. 11Sep13
	labels <- sort(unique(as.vector(res)))[-1]
	n.labels <- length(labels)
  newlabel<-1
  
  for (i in 1:n.labels)
  { 
    while (sum(newlabel==labels)>0)newlabel<-newlabel+1
    if (labels[i]>n.labels){
      if(labels%%10==1)cat(".")
      res[res==labels[i]]<-newlabel
      labels[i]<-newlabel
    }
  }

  cat(", done.")
	
	return(as.Image(res))
}

#' Computes moments from image objects
#' 
#' Computes intensity-weighted centers of objects and their mass (sum of intensities)
#' 
#' @param mask a labeled stack as returned from bwlabel3d
#' @param ref the original image stack
#' @return a matrix with the moments of the objects in the stack
cmoments3d <- function(mask, ref){
	labels <- 1:max(mask)
	ret <- t(sapply(labels, function(x){
						ind <- which(mask == x, arr.ind=T)
						w <- ref[ind]
						return(c(x, apply(ind, 2, weighted.mean, w=w), sum(w)))
					}))
	colnames(ret) <- c("label","m.x","m.y","m.z","w")
	return(ret)
}	


#' Get central moments of objects in a single-channel image stack
#' 
#' Uses the methodology used for segmentation in the RBImage vignette 
#'  (threshhold->opening->fillHull) from all 3 spatial directions and 
#'  overlays these results to get a binary image which is then segmented 
#'  with bwlabel3d. Central moments are extracted with cmoments3d 
#' 
#' @param file the path of the image stack
#' @param threshold  the quantile of intensities used for thresholding if quantile=TRUE 
#' 		  or the intensity value if quantile=FALSE, defaults to the 80% quantile
#' @param threshW, threshH width and height of the moving rectangular window for threshold, defaults to 5.  
#' @param brushsize the brushsize for makeBrush for opening, defaults to 3 
#' @param quantile defaults to TRUE
#' @return a list with the original stack, the labeled stack, and the matrix of central moments of the found objects 
preprocess <- function(file, threshold=.95, threshW = 5,  threshH = 5, brushsize=3, quantile=TRUE){
	stopifnot(require(EBImage))
	
	cat("reading file....")
	#read file
	im <- readImage(file)
	
	cat("thresholding & smoothing image....")
	#threshold
#	if(quantile){
#		thresh <- quantile(im, threshold)
#		while(thresh == 0 &&  threshold < 1){
#			threshold <- threshold + .01
#			thresh <- quantile(im, threshold)
#			cat("\n cutoff was 0. increasing threshold to ", threshold, "\n")
#		}
#	} else {
#		thresh <- threshold
#	} 
	if(quantile){
		thresh <- quantile(im, threshold)
		stopifnot(thresh != 0)
		thresh <- quantile(im, threshold)
	} else {
		thresh <- threshold
	} 
	
	
	
	imThresh <- thresh(im, threshW, threshH, thresh)
	
	#smooth binary image from all directions
	mask1 <- opening(imThresh, makeBrush(brushsize, shape='disc'))
	mask1 <- fillHull(mask1)
	mask2 <- opening(aperm(imThresh, c(1,3,2)), makeBrush(brushsize, shape='disc'))
	mask2 <- fillHull(mask2)
	mask3 <- opening(aperm(imThresh, c(2,3,1)), makeBrush(brushsize, shape='disc'))
	mask3 <- fillHull(mask3)
	mask <- fillHull(mask1 + aperm(mask2, c(1,3,2)) + aperm(mask3, c(3, 1, 2)))
	stopifnot(any(as.logical(mask)))
	
	cat("segmentation....")
	#segmentation
	label <- bwlabel3d(mask)
	
	cat("get moments\n")
	#extract moments 
	mom <- cmoments3d(label, im)
	
	cat("found", max(label), "objects in image.\n")
	return(list(im=im, 
				label=label, 
				moments=mom))	
}


#' Compute cross-type nearest neighbor distances 
#' @param dist a distance matrix, the upper n1 x n1 part contains distances between objects of type 1
#' 			the lower n2 x n2 part contains distances between objects of type 2
#' @param n1, n2  numbers of objects of type 1 and 2 respectively
#' @param w optional weights of the objects (length n1+n2), defaults to equal weights
#' @return  a (n1+n2) x 2 matrix with the cross-type nearest neighbor distances and 
#' 			weights given as the sum of the weights of the involved objects
crossNN <- function(dist, n1, n2, w = rep(1, n1+n2)){
	use <- dist[-(1:n1),-((n1+1):(n1+n2))] #use only lower left block containing the cross type distances
	whereMin <- rbind(
			cbind( n1 + apply(use, 2, which.min), 1:n1), #for each type1 which is closest type2
			cbind((n1+1):(n1+n2), apply(use, 1, which.min))) #for each type2 which is closest type1
	
	cnn <- dist[whereMin]
	
	#use sum of weights of the involved objects as weight for the cross-type distance
	weights <- apply(whereMin, 1, function(x){
				sum(w[x])
			})
	return(cbind(cnn=cnn, w=w))
}

#' Permutation Test for cross-type nearest neighbor distances
#' @param dist a distance matrix, the upper n1 x n1 part contains distances between objects of type 1
#' 			the lower n2 x n2 part contains distances between objects of type 2
#' @param n1, n2  numbers of objects of type 1 and 2 respectively
#' @param w (optional) weights of the objects (length n1+n2)
#' @param B number of permutations to generate
#' @param alternative alternative hypothesis ("less" to test H0:Colocalization )
#' @param returnSample return sampled null distibution
#' @param papply which apply function to use for generating the null distribution, 
#' 		defaults to mclapply if multicore is available, else lapply
#' @param ... additional arguments for papply
#' @return a list with the p.value, the observed weighted mean of the cNN-distances, alternative and (if returnSample) the simulated null dist 
cnnTest <- function(dist, n1, n2, w = rep(1, n1+n2), 
		B = 999, alternative = "less", returnSample = TRUE,  
		papply = if (require("multicore")) mclapply else lapply, 
		...){
	
	teststat <- function(dist, n1, n2, w){
		cnn <- crossNN(dist, n1, n2, w)
		return(weighted.mean(x = cnn[,'cnn'], w = cnn[,'w']))
	}
	
	obs <- teststat(dist, n1, n2, w)
	
	permutations <- replicate(B, sample(n1+n2), simplify = FALSE)
	nulldist <- unlist(papply(permutations, function(x){
						teststat(dist[x, x], n1, n2, w[x])
					}, ...))
	p.value <- switch(alternative, 
			greater = mean(obs < nulldist),
			less  = mean(obs > nulldist),
			two.sided = 0.5 * min(mean(obs < nulldist), mean(obs > nulldist)))
	ret <- list(statistic = "weighted mean of cNN-distances", 
			p.value = p.value, estimate = obs, alternative = alternative)
	ret$sample <- nulldist
	return(ret)
}


#' Permutation Test for cross-type nearest neighbor distances
#' @param im1, im2  image stacks as returned by preprocess
#' @param hres, vres horizontal and vertical resolution of the stacks
#' @param B number of permutations to generate
#' @param alternative alternative hypothesis ("less" to test H0:Colocalization )
#' @param returnSample return sampled null distibution
#' @param ... additional arguments for papply
#' @return a list with the p.value, the observed weighted mean of the cNN-distances
testColoc <- function(im1, im2, hres = 0.1023810, vres = 0.2500000, B=999, alternative = "less", returnSample = TRUE, ...){
	#extract centers and adjust to resolution
	centers <- rbind(im1$moments[,c('m.x','m.y','m.z')], im2$moments[,c('m.x','m.y','m.z')])
	centers <- t(t(centers) * c(rep(hres,2),vres))
	
	n1 <- nrow(im1$moments)
	n2 <- nrow(im2$moments)
	w <- c(im1$moments[,'w'], im2$moments[,'w'])
	dist <- as.matrix(dist(centers)) 
	return(cnnTest(dist, n1, n2, w, B, alternative, returnSample, ...))
}

