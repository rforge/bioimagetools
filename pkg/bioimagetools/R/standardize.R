standardize<-function(img,mask=array(TRUE,dim(img)),N=32,sd=1/6)
{
int<-img[which(mask==1)]
mb<-mean(int,na.rm=TRUE)
sdb<-sd(int,na.rm=TRUE)
img<-img-mb
img<-img/sdb
img<-img/sd # Stauchung
img<-img*N
img<-img+(N/2)+.5
img<-round(img)
img[img<1]<-0
img[img>N]<-N+1
return(img)
}


colors.in.classes<-function(classes,color1,color2=NULL,mask=array(TRUE,dim(img)),N=32,sd1=2,sd2=2,col1="green",col2="red")
{
  no2<-ifelse(is.null(color2),TRUE,FALSE)
  classes<-array(classes,dim(classes))
  color1<-array(color1,dim(color1))
  if(!no2)color2<-array(color2,dim(color2))
  mask<-array(mask,dim(mask))
  
  thresh1<-mean(color1)+sd1*sd(color1)
  if(!no2)thresh2<-mean(color2)+sd2*sd(color2)
  
  t1<-table0(classes,N)
  t1<-t1/sum(t1)
  
  t2<-table0(classes[color1>thresh1],N)
  t2<-t2/sum(t2)
  
  t3<-0
  if(!no2){
    t3<-table0(classes[color2>thresh2],N)
    t3<-t3/sum(t3)
  }
  
  barplot(t1,ylim=c(0,max(c(t1,t2,t3))))
  barplot(t2,xlab=NULL,col=col1,density=40,border=col1,add=TRUE,axes=FALSE)
  if(!no2)barplot(t3,xlab=NULL,col=col2,density=40,border=col2,add=TRUE,axes=FALSE)

  ret<-thresh1
  if(!no2)ret<-c(thresh1,thresh2)
  return(ret)
}
