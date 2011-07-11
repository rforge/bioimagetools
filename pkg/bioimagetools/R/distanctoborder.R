distance2border<-function(points,img.classes,x.microns,y.microns,z.microns,class1,class2=NULL,mask=array(TRUE,dim(img.classes)),plot=FALSE,main="Minimal distance to border", xlab="Distance in Microns", xlim=c(-.3,.3),n=100,stats=TRUE,file=NULL)
{
dims<-dim(img.classes)
X<-dims[1]
Y<-dims[2]
Z<-dims[3]

points.discrete<-data.frame("x"=1+floor(X*points$X/x.microns),"y"=1+floor(Y*points$Y/y.microns),"z"=1+floor(Z*points$Z/z.microns))

cat(".")
valid<-c()

for (i in 1:dim(points.discrete)[1])
if(mask[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]==1)
if(img.classes[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]==class1)
{
valid<-rbind(valid,points[i,])
}
names(valid)<-c("x","y","z")

cat(".")

x<-rep(1:X,Y*Z)
y<-rep(rep(1:Y,each=X),Z)
z<-rep(1:Z,each=X*Y)
which<-(mask==1)

nucleus<-data.frame("x"=x[which],"y"=y[which],"z"=z[which],"class"=(img.classes)[which])
if(is.null(class2))chromatin<-nucleus[nucleus$class!=class1,1:3]
if(!is.null(class2))chromatin<-nucleus[nucleus$class==class2,1:3]

cat(".")

#Translate voxel coordinates to microns
chromatin$x<-(chromatin$x-1)/X*x.microns
chromatin$y<-(chromatin$y-1)/Y*y.microns
chromatin$z<-(chromatin$z-1)/Z*z.microns

abstand1<-apply(valid,1,bioimagetools..find.min.distance,chromatin,c(x.microns/X,y.microns/Y,z.microns/Z))

cat(".")

if (is.null(class2))
{
valid<-c()
for (i in 1:dim(points)[1])
if(mask[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]==1)
if(img.classes[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]!=class1)
{
valid<-rbind(valid,points[i,])
}

}
else
{
valid<-c()
for (i in 1:dim(points)[1])
if(mask[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]==1)
if(classes[points.discrete[i,1],points.discrete[i,2],points.discrete[i,3]]==class2)
{
valid<-rbind(valid,points[i,])
}
}
names(valid)<-c("x","y","z")
cat(".")

chromatin<-nucleus[nucleus$class==class1,1:3]

#Translate voxel coordinates to microns
chromatin$x<-(chromatin$x-1)/X*x.microns
chromatin$y<-(chromatin$y-1)/Y*y.microns
chromatin$z<-(chromatin$z-1)/Z*z.microns

cat(".")
abstand2<-apply(valid,1,bioimagetools..find.min.distance,chromatin,c(x.microns/X,y.microns/Y,z.microns/Z))
abstand<-c(abstand1,-abstand2)

cat(".\n")
if (plot)
{
if(!is.null(file))png(file)
temp<-hist(abstand[abstand<xlim[2]&abstand>xlim[1]],breaks=seq(xlim[1],xlim[2],length=n),main=main,xlab=xlab)
if(stats)text(xlim[2]*.85,max(temp$counts)-2,paste("mean: ",round(mean(1000*abstand),1),"\n median: ",round(median(1000*abstand),1),"\n st.dev.: ",round(sd(1000*abstand),2)))
box()
if(!is.null(file))dev.off()
}
else
{
return(abstand)
}

}

bioimagetools..find.min.distance<-function(point,voxels,microns)
{
x<-y<-z<-rep(0,dim(voxels)[1])
which<-(point[1]<voxels$x)
if(sample(100,1)<5)cat(".")
x[which]<-point[1]-voxels$x[which]
which<-(point[1]>(voxels$x+microns[1]))
x[which]<-voxels$x[which]+microns[1]-point[1]
which<-(point[2]<voxels$y)
y[which]<-point[2]-voxels$y[which]
which<-(point[2]>(voxels$y+microns[2]))
y[which]<-voxels$y[which]+microns[2]-point[2]
which<-(point[3]<voxels$z)
z[which]<-point[3]-voxels$z[which]
which<-(point[3]>(voxels$z+microns[3]))
z[which]<-voxels$z[which]+microns[3]-point[3]
dist<-sqrt(x^2+y^2+z^2)
found<-which(dist==min(dist))
return(dist[found][1])
}
