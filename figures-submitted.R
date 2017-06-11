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
# Figure 4-  First bifurcation figure
####################################################



parms=list(b=1,K=400,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.15,0.2),c(0.2,0.15),c(0.175,0.175)))

k=60 
cs=seq(0,0.6,length=k)
a3s=seq(0,1,length=k)
cond1=matrix(0,k,k)
cond2=matrix(0,k,k)
for(i in 1:k){
  temp=c(cs[i],0.2)
  temp2=c(0.2,cs[i])
  temp3=(temp+temp2)/2
  parms$c=rbind(temp,temp2,temp3)
  for(j in 1:k){
    temp=c(1,2)/25
    temp2=c(2,1)/25
    temp3=c(1,1)/25*(a3s[j]*1+(1-a3s[j])*2)
    parms$a=rbind(temp,temp2,temp3)
    out=coexist(parms)
    cond1[i,j]=out$ratio1
    cond2[i,j]=out$ratio3
  }
}

cols=c(1,1,"darkgray","darkgray")

pdf("figure/bif.pdf",height=6,width=7)
par(cex.lab=1.25,cex.axis=1.25,mar=c(2,4.5,1,1))
layout(matrix(c(1,3,2,3),nrow = 2,ncol=2),heights = c(1,1.5))
h=1.5
out=blahblah(h=h) 
temp=max(which(out$c1/0.2<2))
x=out$c1/0.2
y=cbind(out$p1mins,out$p1maxs,out$a1mins,out$a1maxs)
matplot(x[1:temp],y[1:temp,],type="l",col=cols,lty=c(2,2,1,1),
        xlim=c(0,max(x)),xlab="",ylab="min/max frequency",ylim=c(0,1),
        lwd=2)
matplot(x[(temp+1):length(x)],y[(temp+1):length(x),],type="l",
        col=cols,lty=c(2,2,1,1),add=TRUE,lwd=2)
text(0.5,0.85,"coexistence")
text(1.5,0.5,"exclusion")
text(2.5,0.5,"bistability")
legend("bottomright","A",bty="n",cex=1.25)
abline(v=(2-h)/(h-1),lty=2)
abline(v=2,lty=2)
h=1
out=blahblah(h)
temp=max(which(out$c1/0.2<2))
x=out$c1/0.2
y=cbind(out$p1mins,out$p1maxs,out$a1mins,out$a1maxs)
matplot(x[1:temp],y[1:temp,],type="l",col=cols,lty=c(2,2,1,1),
        xlim=c(0,max(x)),xlab="",ylab="min/max frequency",ylim=c(0,1),
        lwd=2)
matplot(x[(temp+1):length(x)],y[(temp+1):length(x),],type="l",
        col=cols,lty=c(2,2,1,1),add=TRUE,lwd=2)
text(1,0.75,"coexistence")
text(2.5,0.5,"bistability")
legend("bottomright","B",bty="n",cex=1.25)
abline(v=(2-h)/(h-1),lty=2)
abline(v=2,lty=2) 
legend("bottomleft",c("prey allele 1","predator 1"),
       col=cols[c(1,3)],lwd=2,lty=c(2,1),bty="n")
par(mar=c(4.5,4.5,1,2))
contour(cs/0.2,a3s,cond1,levels=c(1,1),
        xlab=expression(paste("pleiotropic mismatch ",alpha[i])),yaxt="n",
        ylab=expression(paste("defensive dominance  ",beta[i])),
        drawlabels = FALSE,
        lwd=2)
contour(cs[1:40]/0.2,a3s,cond2[1:40,],levels=c(1,1),add=TRUE,
        drawlabels = FALSE,lwd=2)
legend("bottomright","C",bty="n",cex=1.25)
axis(side=2,at=c(0,0.25,0.5,0.75,1),labels = c("0",F,"1",F,expression(infinity)))
mtext(side = 4, at = 0.5, text="A",las=2,line=1)
mtext(side = 4, at = 0.75, text="B",las=2,line=1)

abline(h=0.5,lty=2)
abline(h=0.75,lty=2)
text(0.1/0.2,0.8,"coexistence",cex=1.25)
text(0.25/0.2,0.2,"intransitive exclusion",cex=1.25)
text(0.5/0.2,0.6,"bistability",cex=1.25)
dev.off()




####################################################
# Figure 5-  Second bifurcation figure
####################################################



parms=list(b=1,K=400,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.2,0.2),c(0.2,0.2),c(0.2,0.2)))


pdf("figure/density-bif.pdf",height=6,width=7)

par(cex.lab=1.25,cex.axis=1.25,mar=c(4.5,4.5,1,1))
layout(matrix(c(2,3,1,3),nrow = 2,ncol=2),heights = c(1,1.5))
cols=c(1,1,"darkgray","darkgray")

# second orbital bifurcatin diagram with h=0.5
out=blahblahh()
x=out$hs/(1-out$hs)
y=cbind(out$p1mins,out$p1maxs,out$a1mins,out$a1maxs)
matplot(x,y,type="l",col=cols,lty=c(2,2,1,1),
        xlab=expression(paste("defensive dominance  ",beta[i])),
        ylab="min/max frequency",ylim=c(0,1),
        lwd=2)
abline(v=0.9,lty=2)
text(0.85,0.25,"intrans.\n excl.")
text(1.1,0.25,"coexistence")
text(1.18,0.9,"B",bty="n",cex=1.25)


# first orbital bifurcation diagram with h=0.46
out2=blahblahK(0.46)
x=out2$Ks
y=cbind(out2$p1mins,out2$p1maxs,out2$a1mins,out2$a1maxs)
matplot(x,y,type="l",col=cols,lty=c(2,2,1,1),
        xlab=expression(paste("Prey carrying capacity ",K)),ylab="min/max frequency",ylim=c(0,1),
        lwd=2)
text(30,0.5,"coexistence")
text(80,0.5,"intransitive\n exclusion")
abline(v=45,lty=2)
text(95,0.9,"A",bty="n",cex=1.25)
legend("bottomright",c("prey allele 1","predator 1"),
       col=cols[c(1,3)],lwd=2,lty=c(2,1),bty="n")
# 2d bif diagram
k=100 
Ks=seq(12.51,100,length=k)
a3s=seq(0.2,0.6,length=k)
cond1=matrix(0,k,k)
cond2=matrix(0,k,k)
for(i in 1:k){
  parms$K=Ks[i]
  for(j in 1:k){
    temp=c(1,2)/25
    temp2=c(2,1)/25
    temp3=c(1,1)/25*(a3s[j]*1+(1-a3s[j])*2)
    parms$a=rbind(temp,temp2,temp3)
    out=coexist(parms)
    cond1[i,j]=out$ratio1
    cond2[i,j]=out$ratio3
  }
}
par(mar=c(4.5,4.5,1,2))
contour(Ks,a3s,cond2,levels=c(1,1),
        xlab=expression(paste("Prey carrying capacity ",K)),yaxt="n",
        ylab=expression(paste("defensive dominance  ",beta[i])),
        drawlabels = FALSE,
        lwd=2,xlim=c(0,max(Ks)))
abline(v=12.5,lty=1,lwd=2)
axis(side=2,at=c(0.2,0.46,0.5,0.6),labels = c(0.8,0.85,1,1.5))
mtext(side = 3, at = 75, text="B",las=1,line=0)
mtext(side = 4, at = 0.46, text="A",las=2,line=1)
abline(v=75,lty=2)
abline(h=0.46,lty=2)
text(60,0.35,"intransitive exclusion",cex=1.25)
text(60,0.55,"coexistence",cex=1.25)
text(0.5,0.35,"predators\n can't establish",cex=1.25,srt=90)
legend("bottomright","C",bty="n",cex=1.25)

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


