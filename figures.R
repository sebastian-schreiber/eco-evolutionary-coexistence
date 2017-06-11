# Load the base code
source("base-code.R") 


####################################################
# Figure 3 - Sample eco-evolutionary dynamics
# associated with bistability and intransitivities
####################################################

# set parameters, time sequence, colors, line types, etc.
parms=list(b=1,K=200,d=0.2,delta=c(0.1,0.1),a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,c=rbind(c(0.2,0.1),c(0.1,0.2),c(0.15,0.15)))
time=seq(0,200,length=200)
cols=c(1,1,"darkgray","darkgray")
ltys=c(1,2,1,2) 
initialn=c(12.497123 ,0.0,0.01,0.01,9.219349)
start=length(time)*0+1
end=length(time)

# set up layout for figure
pdf("figure/bistable.pdf",height=6,width=6)
par(cex.lab=1.25,cex.axis=1.25,bty="n",mar=c(4.5,4.5,0.5,1))
layout(rbind(c(1,1,2,2),c(3,3,4,4),c(5,5,5,6)),heights = c(1,1,1,0.5))

# create first panel
out=ode(y=initialn,times=time,func=diploid,parms=parms)
dens=cbind(out[start:end,2]*2+out[start:end,4],
           out[start:end,3]*2+out[start:end,4],
           out[start:end,5],out[start:end,6])
matplot(time[start:end]*parms$d,dens,type="l",
        col=cols,lwd=2,lty=ltys,
        xlab="prey generations",ylab="densities",bty="n")
legend("topleft","A",bty="n",cex=1.25)
text(30,10,"Eco-evolutionary\n bistability",cex=1.25)


# create second panel
initialn=c(0,12.497123,0.01,9.219349,0.01)
start=length(time)*0+1
end=length(time)
out=ode(y=initialn,times=time,func=diploid,parms=parms)
dens=cbind(out[start:end,2]*2+out[start:end,4],
           out[start:end,3]*2+out[start:end,4],
           out[start:end,5],out[start:end,6])
matplot(time[start:end]*parms$d,dens,type="l",
        col=cols,lwd=2,lty=ltys,
        xlab="prey generations",ylab="densities",bty="n")
#legend("bottomright","B",bty="n",cex=1.25)

# create third panel
parms=list(b=1,K=200,d=0.2,delta=c(0.1,0.12),a=rbind(c(2,1),c(1,2),c(1.5,1.5))/25,c=rbind(c(0.15,0.2),c(0.2,0.15),c(0.16,0.175)))
time=seq(0,500,length=400)
parms$a=rbind(c(2,1),c(1,2),c(1,1))/10
out=ode(y=initialn,times=time,func=diploid,parms=parms)
start=length(time)*0+1
end=length(time)
par(cex.lab=1.25,cex.axis=1.25,bty="n")
dens=cbind(out[start:end,2]*2+out[start:end,4],
           out[start:end,3]*2+out[start:end,4],
           out[start:end,5],out[start:end,6])
matplot(time[start:end]*parms$d,dens,type="l",
        col=cols,lwd=2,lty=ltys,
	xlab="prey generations",ylab="densities",bty="n")
legend("topleft","B",bty="n",cex=1.25)
text(70,15,"Coexistence by\n unstable intransitivity",cex=1.25)

# create fourth panel
parms=list(b=0.5,K=100,d=0.2,delta=c(0.1,0.12),a=rbind(c(2,1),c(1,2),c(1.5,1.5))/10,c=rbind(c(0.2,0.1),c(0.1,0.2),c(0.15,0.15)))
time=seq(0,2000,length=400)
initialn=c(2.5,0.1,0,2.18,0.1)
start=length(time)*0.75+1
end=length(time)
out=ode(y=initialn,times=time,func=diploid,parms=parms)
dens=cbind(out[start:end,2]*2+out[start:end,4],
           out[start:end,3]*2+out[start:end,4],
           out[start:end,5],out[start:end,6])
matplot(time[start:end]*parms$d,dens,type="l",
        col=cols,lwd=2,lty=ltys,
	xlab="prey generations",ylab="densities",bty="n")
#legend("topright","D",bty="n",cex=1.25)


# create fifth panel
parms=list(b=1,K=200,d=0.2,delta=c(0.1,0.12),a=rbind(c(2,1),c(1,2),c(1.5,1.5))/25,c=rbind(c(0.15,0.2),c(0.2,0.15),c(0.16,0.175)))
time=seq(0,3100,length=4200)
cols=c(1,1,"darkgray","darkgray")
ltys=c(1,2,1,2)
initialn=c(8.333333,0,0.1,9.479167,0.1)
start=length(time)*0+1
end=length(time)
out=ode(y=initialn,times=time,func=diploid,parms=parms)
dens=cbind(out[start:end,2]*2+out[start:end,4],
           out[start:end,3]*2+out[start:end,4],
           out[start:end,5],out[start:end,6])
matplot(time[start:end]*parms$d,dens,type="l",
        col=cols,lwd=2,lty=ltys,
        xlab="prey generations",ylab="density",bty="n",ylim=c(0,36))
legend("topleft","C",bty="n",cex=1.25)
text(520,30,"Exclusion\n by stable\n intransitivity",cex=1.25)


# create legend in fifth panel
par(mar=c(0,0,0,0))
plot(0,col="white",xaxt="n",yaxt="n",xlab="",ylab="")
legend("center",c("prey allele 1","prey allele 2","predator 1","predator 2"),col=cols,lwd=3,lty=ltys,bty="n",seg.len=4,cex=1.25)

# close pdf device

dev.off()

####################################################
# Numerical version of Figure 4.
# Final version was edited with PowerPoint 
####################################################


parms=list(b=1,K=4000,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.15,0.2),c(0.2,0.15),c(0.175,0.175)))

k=60 
cs=exp(seq(log(0.6),log(0.1),length=k))
a3s=exp(seq(log(0.1),log(0.9),length=k))
pleio=numeric(k)
dominance=numeric(k)
cond1=matrix(0,k,k)
cond2=matrix(0,k,k)
for(i in 1:k){
  temp=c(cs[i],0.2)
  temp2=c(0.2,cs[i])
  temp3=(temp+temp2)/2
  parms$c=rbind(temp,temp2,temp3)
  pleio[i]=-log(temp[1]/temp[2])
  for(j in 1:k){
    temp=c(1,2)/25
    temp2=c(2,1)/25
    temp3=c(1,1)/25*(a3s[j]*1+(1-a3s[j])*2)
    dominance[j]=log((temp[2]-temp3[2])/(temp3[2]-temp2[2]))
    parms$a=rbind(temp,temp2,temp3)
    out=coexist(parms)
    cond1[i,j]=out$ratio1
    cond2[i,j]=out$ratio3
  }
}

cols=c(1,1,"darkgray","darkgray")


par(mar=c(4.5,4.5,1,2))
contour(pleio,dominance,cond1,levels=c(1,1),
        xlab=expression(paste("ecological pleiotropy ",alpha[i])),yaxt="n",
        ylab=expression(paste("defensive dominance  ",beta[i])),
        drawlabels = FALSE,
        lwd=2)

temp=which(cond1[,1]>1)

contour(pleio[temp],dominance,cond2[temp,],levels=c(1,1),add=TRUE,
        drawlabels = FALSE,lwd=2)
abline(h=0,lty=2)
abline(v=0,lty=2)




####################################################
# Figure 5-  Four bifurcation figures
####################################################



parms=list(b=1,K=400,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.2,0.2),c(0.2,0.2),c(0.2,0.2)))

# Varying ecological pleiotropy
h=1.5
out1=blahblah(h=h)
# Varying dominance and then K
out2=blahblahh()
out3=blahblahK(0.46)

pdf("bifurcations.pdf",height=5,width=5)
par(cex.lab=1.25,cex.axis=1.25,mar=c(4,4.5,2,1))
layout(matrix(c(1,2,3),nrow = 3,ncol=1),heights = c(1.2,1.2,1.2))

x=-log(out1$c1/0.2)
y=cbind(out1$p1mins,out1$p1maxs,out1$a1mins,out1$a1maxs)
temp=max(which(out1$p1min>0.6))+3
temp2=length(out1$p1min)
matplot(x,y[,c(2,4)],type="l",col=cols[c(1,3)],lty=c(2,1),
        xlab=expression(paste("mean pleiotropy ",(alpha[1]+alpha[2])/2)),
        ylab="",ylim=c(0,1),
        lwd=2)
matplot(x[temp:temp2],y[temp:temp2,c(1,3)],add=TRUE,
        col=cols[c(1,3)],lty=c(2,1),type="l",lwd=2)
legend("bottomleft",c("prey allele 1","predator 1"),
       col=cols[c(1,3)],lwd=2,lty=c(2,1),bty="n")
abline(v=x[temp],lty=2)
abline(v=-0.12,lty=2)
text(-1,0.8,"bistability")
text(-0.4,0.5,"intransitive\nexclusion")
text(0.25,0.9,"coexistence")
legend("topright","A",bty="n",cex=1.25)

x=log(out2$hs/(1-out2$hs))
y=cbind(out2$p1mins,out2$p1maxs,out2$a1mins,out2$a1maxs)
matplot(x,y,type="l",col=cols,lty=c(2,2,1,1),
        xlab=expression(paste("mean dominance  ",(beta[1]+beta[2])/2)),
        ylab="min/max frequency",ylim=c(0,1),
        lwd=2)
abline(v=log(0.9),lty=2)
text(log(0.85),0.5,"intransitive\nexclusion.")
text(log(1.1),0.9,"coexistence")
legend("topright","B",bty="n",cex=1.25)


x=out3$Ks
y=cbind(out3$p1mins,out3$p1maxs,out3$a1mins,out3$a1maxs)
matplot(x,y,type="l",col=cols,lty=c(2,2,1,1),
         xlab=expression(paste("Prey carrying capacity ",K)),
         ylab="",ylim=c(0,1),
         lwd=2)
text(30,0.5,"coexistence")
text(80,0.5,"intransitive\n exclusion")
abline(v=45,lty=2)
legend("topright","C",bty="n",cex=1.25)
dev.off()


#####################################
# APPENDIX MUTATION FIGURE
######################################

parms=list(b=0.5,K=500,d=0.05,delta=c(0.1,0.15),
           a=rbind(c(2,1),c(1,2))/10,
           c=rbind(c(0.2,0.1),c(0.1,0.2)),
           mu=0)
parms$c=rbind(c(0.18,0.2),c(0.2,0.18))
initialn=c(5,0.1,2.125,0.1)
time=seq(0,10000,length=10000)
blerg=function(mu){
  parms$mu=mu
  out=ode(y=initialn,times=time,func=haploid,parms=parms)
  n=out[5000:10000,2:3]
  x=n[,1]/(n[,1]+n[,2])
  n=out[5000:10000,4:5]
  y=n[,1]/(n[,1]+n[,2])
  return(list(xmin=min(x),xmax=max(x),xmean=mean(x),
    ymin=min(y),ymax=max(y),ymean=mean(y)))
}

k=25
logmus=seq(-50,-1,length=k)
stuff=matrix(0,k,6)
for(i in 1:k){
  out=blerg(10^logmus[i])
  stuff[i,]=as.numeric(out)
}

par(cex.lab=1.25,cex.axis=1.25,mfrow=c(1,2))
matplot(logmus,log10(stuff[,1:3]),type="l",col=c(1,1,1,2,2,2),
        lty=c(2,2,1),lwd=2,
        xlab=expression(paste(log[10]," mutation rate")),
        ylab=expression(paste(log[10]," frequency of prey allele 1")))
legend("bottomright",c("max","mean","min"),lty=c(2,1,2),lwd=2,bty="n")
matplot(logmus,log10(stuff[,4:6]),type="l",col=c(1,1,1,2,2,2),
        lty=c(2,2,1),lwd=2,
        xlab=expression(paste(log[10]," mutation rate")),
        ylab=expression(paste(log[10]," frequency of predator 1")))


