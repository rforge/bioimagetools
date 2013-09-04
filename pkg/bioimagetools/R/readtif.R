if(0)
  {
  readTIF<-function(file=file.choose(),native=FALSE)
{
  require(tiff)
  li<-readTIFF(file,all=TRUE,info=TRUE,as.is=TRUE,native=native)
  Z<-length(li)
  img<-array(0,c(dim(li[[1]]),Z))
  for (i in 1:Z)img[,,i]<-li[[i]]
  storage.mode(img)<-"integer"
  temp<-attributes(li[[1]])
  tmp<-gregexpr("\n",temp$description)
  if(length(tmp)>0)if (tmp[[1]][1]!=-1)
    {
      temp2<-regmatches(temp$description,tmp,invert=TRUE)[[1]]
      temp3<-temp4<-c()
      for (i in temp2)
      {
        j<-gregexpr("=",i)[[1]]
        j<-regmatches(i,j,invert=TRUE)[[1]]
        temp3<-c(temp3,j[1])
        temp4<-c(temp4,j[2])
      }
      names(temp4)<-temp3
      temp<-c(temp,temp4)
  }
  temp<-temp[!(names(temp)=="")]
  K<-as.integer(temp$channels)
  if(length(K)==0)K<-1
  if (K>1)
  {
    img<-array(img,c(dim(li[[1]]),K,Z/K))
    storage.mode(img)<-"integer"
    #for (i in 1:K)img0[,,i,]<-img[,,seq(i,Z,by=K)]
    #img<-img0
  }
  if(min(img)<0){img<-array(bitFlip(img,bitWidth=temp$bits.per.sample),dim(img))}
  img<-aperm(img,c(2,1,3:length(dim(img))))
  temp$dim<-dim(img)
  attributes(img)<-temp
  return(img)
}

writeTIF<-function(img,file,bps=NULL)
{
  require(tiff)
  if(is.null(bps))if(!is.null(attr(img,"bits.per.sample")))bps<-attr(img,"bits.per.sample")
  if(is.null(bps))bps<-8L
  imglist<-list()
  if (length(dim(img))==3)
  {
    Z<-dim(img)[3]
    for (i in 1:Z)
    imglist[[i]]<-img[,,i]/max(img[,,i])
  }
  if (length(dim(img))==4)
  {
    C<-dim(img)[3]
    Z<-dim(img)[4]
    k<-0
    maxi<-1:C
    for (j in 1:C)maxi[i]<-max(img[,,j,],na.rm=TRUE)
    for (i in 1:Z)
      for (j in 1:C)
        {
        k<-k+1
        imglist[[k]]<-img[,,j,i]/maxi[j]
        }
  }
  Z<-length(imglist)
  ati<-attributes(img)
  ati$dim<-dim(imglist[[1]])
  for (i in 1:Z)
    attributes(imglist[[i]])<-ati
  writeTIFF(what=imglist,where=file,reduce=TRUE,bits.per.sample=bps)
}
}