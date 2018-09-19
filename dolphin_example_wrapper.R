############################################################
# Dolphin Case Study Wrapper
############################################################
# 1.) Section 1: Libraries and functions
# 2.) Section 2: JAGS Model
# 3.) Section 3: Data Processing and JAGS run
# 4.) Section 4: Results post processing
# 5.) Section 5: Graphing

############################################################
# Section 1: libraries and functions
############################################################
# load needed libraries
library(MASS) 
library(clusterGeneration)
library(repmis)
library(rgl)
library(shape)
library(magick)
library(jagsUI) # to run jags


#Function to calcuate SEV
SEV.function<-function(sigma) {
  if (eigen(sigma)$values[3] > 0)
  {
    axises<-sqrt(eigen(sigma)$values)
  }
  else {
    axises<-NA
  }
  #out<-list()
  #out$a<-axises[1]
  #out$b<-axises[2]
  #out$c<-axises[3]
  #out$SEV<-prod(axises,pi,4/3)
  return(prod(axises,pi))
}
# Calculate convex hull
hull<-function (x, y) 
{
  chI <- chull(x, y)
  chI <- c(chI, chI[1])
  hullx <- x[chI]
  hully <- y[chI]
  out<- list()
  out$x<-hullx
  out$y<-hully
  out
}

# calculate standard ellipse area
SEA.function<-function(sigma) {
  axises<-sqrt(eigen(sigma)$values)
  out<-list()
  out$a<-axises[1]
  out$b<-axises[2]
  out$SEV<-prod(axises,pi)
  return(out)
}

# calculate distance between two points in 3D
distance<-function(a,b){
  l<-sqrt((a[,1]-b[,1])^2+(a[,2]-b[,2])^2)
  return(l)
}

# shortcut for quantile calculation
quants<-function(x){
  quantile(x,probs=c(0.025,0.5,0.975),na.rm=TRUE)
}

# calculate mode
mode <- function(s) {
  d <- density(s,na.rm=TRUE)
  d$x[which.max(d$y)]
}

# Functions for ellipsoid overlap
MVNdist <- function(mu,S){function(x){t(x-mu)%*%solve(S)%*%(x-mu)}}
CubeCount3D <- function(mu1,mu2,S1,S2,HowFine=1.0,ShowProgress=TRUE){
  # set-up stuff:
  OneStep <- 10^(-HowFine)
  CubeVol <- OneStep^3
  CubeCt <- 0
  A1  <- eigen(S1)$vect
  dD1 <- eigen(S1)$val  
  MVNdist2 <- MVNdist(mu2,S2)
  
  # the loops:
  zlist <- seq(from=0,to=dD1[2]^.5,by=OneStep)
  zlist <- if(length(zlist)>1){c(rev(-zlist),zlist[2:length(zlist)])}
  for(Z in zlist){
    yBd <- (dD1[2]*(1-Z^2/dD1[3]))^.5
    # if(is.nan(yBd)|is.na(yBd)){errorlog <- c(errorlog,list(yBd,Z))}
    ylist <- seq(from=0,to=yBd,by=OneStep)
    ylist <- if(length(ylist)>1){c(rev(-ylist),ylist[2:length(ylist)])}
    for(Y in ylist){
      xBd <- (dD1[1]*(1-Y^2/dD1[2]-Z^2/dD1[3]))^.5
      # if(is.nan(xBd)|is.na(xBd)){errorlog <- c(errorlog,list(xBd,Y,Z))}
      xlist <- seq(from=0,to=xBd,by=OneStep)
      xlist <- if(length(xlist)>1){c(rev(-xlist),xlist[2:length(xlist)])}
      for(X in xlist){
        if(MVNdist2(A1%*%c(X,Y,Z)+mu1)<1){CubeCt=CubeCt+1}
      }
    }
    # This will show a "progress report" in the console window as this script is running:
    # (0=start, 1=end)
    if(ShowProgress){print((.5*(Z+dD1[3]^.5))/dD1[3]^.5)}
  }
  return(CubeCt*CubeVol)
} 

############################################################
# Section 2: Model File
############################################################
sink("Dolphin_example.txt")
cat("
    model{
    #priors
    for (k in 1:nGroups){
      # Mean Vector
      mu[k,1] ~ dnorm(0,0.001)
      mu[k,2] ~ dnorm(0,0.001)
      mu[k,3] ~ dnorm(0,0.001)
      # Covariance Matrix
      prec[1:3,1:3,k] ~ dwish(S3,4)
      cov[1:3,1:3,k] <- inverse(prec[1:3,1:3,k])
    }#k
    
    #likelihood
    for (i in 1:obs){
      y[i,1:3] ~ dmnorm(mu[stock[i],],prec[,,stock[i]])
    }#i
  }
    
    ",fill = TRUE)
sink()


############################################################
# Section 3: Data Processing and JAGS Run
############################################################
# Read in data
# Other data files can be used but formated as [observations, c(group, isotope1, isotope2....)] 
dolphin<-as.data.frame(matrix(c(
               "G",	-10.533,	11.513,	7.2,
               "G",	-11.665,	13.495,	8.7,
               "G",	-12.66,	11.928,	9.2,
               "G",	-13.734,	12.814,	9.9,
               "G",	-12.61,	12.736,	10.1,
               "G",	-11.723,	12.877,	11.6,
               "G",	-11.223,	12.179,	12.7,
               "G",	-10.953,	12.787,	12.8,
               "G",	-12.7,	12.132,	15,
               "G",	-12.381,	12.537,	15.7,
               "IN",	-10.482,	11.049,	4.5,
               "IN",	-11.727,	11.956,	6.2,
               "IN",	-10.53,	10.82,	7.9,
               "IN",	-11.592,	12.56,	10.2,
               "IN",	-10.54,	11.69,	5.5,
               "IN",	-13.63,	11.794,	5.8,
               "IN",	-8.23,	10.906,	6,
               "IN",	-11.47,	11.83,	6.1,
               "IN",	-8.27,	10.561,	6.1,
               "IN",	-8.48,	10.33,	6.4,
               "IN",	-8.797,	12.949,	6.7,
               "IN",	-8.85,	13.694,	7,
               "IN",	-10,	11.3,	7.5,
               "IN",	-12.2,	11.6,	8.3,
               "IN",	-12.413,	15.838,	8.8,
               "IN",	-12.584,	14.672,	8.8,
               "IN",	-12.22,	12.57,	9.2,
               "OFF",	-11.085,	11.833,	7,
               "OFF",	-12.134,	12.764,	13.8,
               "OFF",	-10.447,	12.503,	15.7,
               "OFF",	-13.992,	14.04,	15.9,
               "OFF",	-12.46,	12.53,	15.9,
               "OFF",	-11.392,	14.743,	16.8,
               "OFF",	-11.237,	13.312,	17,
               "OFF",	-12.637,	13.277,	17.8,
               "OFF",	-12.941,	13.578,	18.4), ncol=4, byrow = TRUE))

colnames(dolphin)<-c("Stock",	"d13C",	"d15N",	"d34S")
dolphin[,2:ncol(dolphin)]<-apply(dolphin[,2:ncol(dolphin)], 2, as.numeric)
# Make objects for jags
nIsotopes<-ncol(dolphin)-1
isotope.names<-colnames(dolphin)[-1]
stock.names<-unique(dolphin$Stock)
nGroups<-length(stock.names)
obs<-nrow(dolphin)
stock<-as.numeric(as.factor(dolphin$Stock))

#package the data
jags.data=list(y=as.matrix(dolphin[,2:(nIsotopes+1)]), 
               nGroups=nGroups, obs=obs, S3=(diag(3)*3),
               stock=stock)

# Provide intial values for covariance matrix
init.prec<-array(NA, dim=c(nIsotopes,nIsotopes,nGroups))
for (i in 1:nGroups){
    temp<-colMeans(apply(dolphin[,2:(nIsotopes+1)], 2, as.numeric))
    init.prec[,,i]<-cov(cbind(rnorm(5,temp[1],5),rnorm(5,temp[2],5),rnorm(5,temp[3],5)))
}
inits=function(){list( prec=init.prec)}

# Parameters to save
params<-c("mu","cov")

# MCMC settings
nt <- 5
nb <- 10000
nc <- 3
na<-10000
ni <- 20000 

# Run Jags
out <- jags(data = jags.data,
            #inits = inits,
            parameters.to.save = params,
            model.file = "Dolphin_example.txt",
            n.chains = nc,
            n.adapt = na,
            n.iter = ni,
            n.burnin = nb,
            n.thin = nt#,parallel=TRUE
)

# Check trace
#traceplot(out)

# Summary
out

############################################################
# Section 4: Results Post Processing
############################################################
# Calculate SEV
SEV<-matrix(NA,ncol=nGroups,nrow=nrow(out$sims.list$cov))
SEV.summary<-matrix(NA,nrow=nGroups,ncol=5)
colnames(SEV.summary)<-c("n","2.5","50.0","97.5","mode")
for (k in 1:nGroups){
  SEV[,k]<-apply(out$sims.list$cov[,,,k],c(1),SEV.function)
  SEV.summary[k,]<-c(sum(dolphin$Stock==stock.names[k]),quants(SEV[,k]),mode(SEV[,k]))
}
rownames(SEV.summary)<-stock.names
colnames(SEV.summary)<-c("n","2.5","50","97.5","mode")
print(SEV.summary,digits=3)

# SEV Pairwise comparison
SEV.diff<-rep(list(matrix(NA,ncol=nGroups,nrow=nrow(SEV))),nGroups)
for (i in 1:nGroups){
  for (j in 1:nGroups){
    SEV.diff[[i]][,j]<-SEV[,i]-SEV[,j]
  }
}
diff.table<-matrix(NA,nrow=nGroups,ncol=3)

for (i in 1:nGroups){
  for (j in 1:nGroups){
    diff.table[i,j]<-sum(SEV.diff[[i]][,j]>0)/length(SEV.diff[[i]][,j])
  }
}
p.table<-matrix(NA,nrow=nGroups,ncol=nGroups)
colnames(p.table)<-paste(">",stock.names, sep="")
rownames(p.table)<-stock.names
for (i in 1:nGroups){
  for (j in 1:nGroups){
    p.table[i,j]<-sum(SEV.diff[[i]][,j]>0)/length(SEV.diff[[i]][,j])
  }
}
print(p.table)


#Pairwise comparison centroid 
mu.diff<-array(NA, dim=c(nIsotopes,nrow(out$sims.list$mu),nGroups,nGroups))
for (k in 1:nIsotopes){  
  for (i in 1:nrow(out$sims.list$mu)){
    for (j in 1:nGroups){
      for (m in 1:nGroups){
        mu.diff[k,i,j,m]<-out$sims.list$mu[i,j,k]-out$sims.list$mu[i,m,k]
      }}}}
mu.p<-rep(list(matrix(NA,ncol=nGroups,nrow=nGroups)),3)
for (k in 1:3){
  for(j in 1:nGroups){
    for(m in 1:nGroups){
      mu.p[[k]][j,m]<-sum(mu.diff[k,,j,m]>0)/nrow(out$sims.list$mu)
    }}}
for (i in 1:nGroups){
  rownames(mu.p[[i]])<-stock.names
  colnames(mu.p[[i]])<-paste(">",stock.names, sep="")
}
names(mu.p)<-isotope.names
mu.p

# Distance between centroids
null_mu<-array(out$sims.list$mu[seq(1,nrow(out$sims.list$mu),by=2),,],dim=c(nrow(out$sims.list$mu)/2,3,nGroups))
test_mu<-array(out$sims.list$mu[seq(2,nrow(out$sims.list$mu),by=2),,],dim=c(nrow(out$sims.list$mu)/2,3,nGroups))
dist.diff1<-array(NA, dim=c(nGroups,nGroups,nrow(out$sims.list$mu)/2))
for (j in 1:nGroups){
  for (m in 1:nGroups){
    dist.diff1[j,m,]<-distance(test_mu[,j,],test_mu[,m,]) - distance(test_mu[,j,],null_mu[,j,]) - distance(test_mu[,m,],null_mu[,m,])
  }}
dist.diff2<-array(NA, dim=c(nGroups,nGroups,nrow(out$sims.list$mu)/2))
for (j in 1:nGroups){
  for (m in 1:nGroups){
    dist.diff2[j,m,]<-distance(test_mu[,j,],test_mu[,m,])
  }}
dist.sum<-rep(list(matrix(NA,ncol=nGroups,nrow=nGroups)),3)
for (i in 1:3){
  for (j in 1:nGroups){
    for (m in 1:nGroups){
      dist.sum[[i]][j,m]<-quantile(dist.diff2[j,m,],probs=c(0.025,0.5,0.975)[i])
    }}}

names(dist.sum)<-c("q2.5","q50","q97.5")
for(i in 1:3){
  colnames(dist.sum[[i]])<-stock.names
  rownames(dist.sum[[i]])<-stock.names
}
pdist<-matrix(NA, nrow=nGroups,ncol=nGroups)
for(j in 1:nGroups){
  for(m in 1:nGroups){
    pdist[j,m]<-sum(dist.diff1[j,m,]>0)/(nrow(out$sims.list$mu)/2)
  }
}
colnames(pdist)<-row.names(pdist)<-stock.names
print(pdist, digits=2)
print(dist.sum, digits=2)

############################################################
# Section 5: Graphs
############################################################
# 3D Graphs
#get mean covariance matrix and mean vector
sigma<-out$mean$cov
mu<-out$mean$mu

# Call and construct plot
colorit<-c("red","yellow","blue")
open3d()
for (k in nGroups:1){
  nam <- paste("A", k, sep = "") 
  assign(nam, ellipse3d(qmesh=TRUE,sigma[,,k],centre=out$mean$mu[k,],level=0.40, subdivide = 5, trans=diag(4)))
  wire3d(get(paste("A", k, sep="")), col=colorit[k])
  nam2 <- paste("cords", k, sep = "")
  assign(nam2,t(get(paste("A", k, sep=""))$vb)[,1:3])
}
axes3d(lwd=5, col="black")
grid3d("x+",lwd=5,col="black")
grid3d("y+",lwd=5,col="black")
grid3d("z+",lwd=5,col="black")
axes3d(c('x--','x++','y--','y++','z--','z++'),expand=2)
title3d(,,'Carbon','Nitrogen','Sulfur')

# Save a movie rotation
M <- par3d("userMatrix")
if (!rgl.useNULL())
  play3d( par3dinterp(time=(0:2)*0.75,userMatrix=list(M,
                                                      rotate3d(M, pi/2, 1, 0, 0),
                                                      rotate3d(M, pi/2, 0, 1, 0) ) ), 
          duration=3 )
movie3d( spin3d(axis=c(0,1,0),rpm=3), duration=10, dir=getwd() )


# 2D Graphs (must run 3d graph first)
cords<-data.frame(matrix(NA, nrow=nrow(cords1),ncol=nGroups*3))
seq3<-seq(1,nGroups*3,by=3)
for (i  in 1:nGroups){
  cords[,seq3[i]:(seq3[i]+2)]<-get(paste("cords",i,sep=""))
}

cmin<-floor(min(cords[,seq3]))
cmax<-ceiling(max(cords[,seq3]))
nmin<-floor(min(cords[,(seq3+1)]))
nmax<-ceiling(max(cords[,(seq3+1)]))
smin<-floor(min(cords[,(seq3+2)]))-1
smax<-ceiling(max(cords[,(seq3+2)]))
carbon<-expression(paste(delta^"13"*"C (\211)"))
nitrogen<-expression(paste(delta^"15"*"N (\211)"))
sulfur<-expression(paste(delta^"34"*"S (\211)"))
rgb<-matrix(c(1,0,0,0,1,0,0,0,1),ncol=3)

#Carbon Nitrogen
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,ylim=c(nmin,nmax),xlim=c(cmin,cmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=1),col="white")
abline(h=seq(-100,100,by=1),col="white")
axis(side = 1, at = seq(-100,100,by=1), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=1),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(carbon, side=1, line=3.3, cex=2)
par(las=3)
mtext(nitrogen, side=2, line=1.8, cex=2)
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col=rgb(rgb[i,1],rgb[i,2],rgb[i,3],0.5),border=NA)
}

#Carbon Sulfur
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,ylim=c(smin,smax),xlim=c(cmin,cmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=1),col="white")
abline(h=seq(-100,100,by=2),col="white")
axis(side = 1, at = seq(-100,100,by=1), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=2),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(carbon, side=1, line=3.3, cex=2)
par(las=3)
mtext(sulfur, side=2, line=1.8, cex=2)

for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+2)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+2)]),col=rgb(rgb[i,1],rgb[i,2],rgb[i,3],0.5),border=NA)
}

# Sulfur Nitrogen
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,xlim=c(smin,smax),ylim=c(nmin,nmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=2),col="white")
abline(h=seq(-100,100,by=1),col="white")
axis(side = 1, at = seq(-100,100,by=2), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=1),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(sulfur, side=1, line=3.3, cex=2)
par(las=3)
mtext(nitrogen, side=2, line=1.8, cex=2)
for(i in 1:nGroups){
  polygon(hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col=rgb(rgb[i,1],rgb[i,2],rgb[i,3],0.5),border=NA)
}



