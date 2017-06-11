# Load the base code
source("base-code.R") 

parms=list(b=1,K=40,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.15,0.2),c(0.2,0.15),c(0.175,0.175)))


blahblah=function(h=1.5){
  blah=function(c1){
    parms$a=rbind(c(1,2),c(2,1),c(h,h))/25
    c2=0.2
    parms$c=rbind(c(c1,c2),c(c2,c1),c(c1+c2,c1+c2)/2)
    time=seq(0,6000,length=1000)
    initialn=c(1,0.1,0,1,0.1)
    out=ode(y=initialn,times=time,func=diploid,parms=parms)
    start=length(time)/3
    end=length(time)
    out=out[start:end,]
    dens=cbind(out[,2]*2+out[,4],
               out[,3]*2+out[,4],
               out[,5],out[,6])
    return(list(a1min=min(dens[,1]/(dens[,1]+dens[,2])),
                a1max=max(dens[,1]/(dens[,1]+dens[,2])),
                a1mean=mean(dens[,1]/(dens[,1]+dens[,2])),
                p1min=min(dens[,3]/(dens[,3]+dens[,4])),
                p1max=max(dens[,3]/(dens[,3]+dens[,4])),
                p1mean=mean(dens[,3]/(dens[,3]+dens[,4])))
    )
  }
  
  k=30
  c1s=exp(seq(log(0.6),log(0.1),length=k))
  p1mins=numeric(k)
  p1maxs=numeric(k)
  p1means=numeric(k)
  a1mins=numeric(k)
  a1maxs=numeric(k)
  a1means=numeric(k)
  for(i in 1:k){
    out=blah(c1s[i])
    a1mins[i]=out$a1min
    a1maxs[i]=out$a1max
    a1means[i]=out$a1mean
    p1mins[i]=out$p1min
    p1maxs[i]=out$p1max
    p1means[i]=out$p1mean
  }
  
  return(list(a1mins=a1mins,a1maxs=a1maxs,p1mins=p1mins,
              p1maxs=p1maxs,c1s=c1s))
}


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


h=1.5
out1=blahblah(h=h)
h=1.0
out2=blahblah(h=h) 



#load("Figure4.Rdata")
par(cex.lab=1.25,cex.axis=1.25,mar=c(2,4.5,1,1))
layout(matrix(c(1,2,3),nrow = 3,ncol=1),heights = c(1,1,2.5))
 
x=out1$c1/0.2
y=cbind(out1$p1mins,out1$p1maxs,out1$a1mins,out1$a1maxs)
matplot(-log(x),y,type="l",col=cols,lty=c(2,2,1,1),
        xlab="",ylab="min/max frequency",ylim=c(0,1),
        lwd=2)

x=out2$c1/0.2
y=cbind(out2$p1mins,out2$p1maxs,out2$a1mins,out2$a1maxs)
matplot(-log(x),y,type="l",col=cols,lty=c(2,2,1,1),
        xlab="",ylab="min/max frequency",ylim=c(0,1),
        lwd=2)

par(mar=c(4.5,4.5,1,2))
contour(pleio,dominance,cond1,levels=c(1,1),
        xlab=expression(paste("ecological pleiotropy ",log(alpha[i]))),yaxt="n",
        ylab=expression(paste("defensive dominance  ",log(beta[i]))),
        drawlabels = FALSE,
        lwd=2)

temp=which(cond1[,1]>1)

contour(pleio[temp],dominance,cond2[temp,],levels=c(1,1),add=TRUE,
        drawlabels = FALSE,lwd=2)
legend("bottomright","C",bty="n",cex=1.25)
abline(h=0,lty=2)
abline(v=0,lty=2)

#save.image("Figure4.Rdata")






parms=list(b=1,K=400,d=0.2,delta=c(0.1,0.1),
           a=rbind(c(1,2),c(2,1),c(1.5,1.5))/25,
           c=rbind(c(0.2,0.2),c(0.2,0.2),c(0.2,0.2)))


k=10 
Ks=seq(12.51,100,length=k)
a3s=seq(0.2,0.6,length=k)
cond1=matrix(0,k,k)
cond2=matrix(0,k,k)
dominance=numeric(k)
for(i in 1:k){
  parms$K=Ks[i]
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
par(mar=c(4.5,4.5,1,2))
contour(Ks,dominance,cond2,levels=c(1,1),
        xlab=expression(paste("Prey carrying capacity ",K)),
        ylab=expression(paste("defensive dominance  ",beta[i])),
        lwd=2,xlim=c(0,max(Ks)))
abline(v=12.5,lty=1,lwd=2)
abline(h=0,lty=2,lwd=1)
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