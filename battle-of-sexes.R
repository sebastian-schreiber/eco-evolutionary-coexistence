# simulator for asymmetric games via replicator equations

parms=list(A=rbind(c(0,1),c(1,0)),B=-rbind(c(0,1),c(1,0)))


f=function(t,x,parms){
  with(parms,{
    dx1=x[1]*(1-x[1])*(A[1,2]+0.2*x[1]-(A[1,2]+A[2,1])*x[2])   #use 0.2*x[1], -0.2*x[1]
    dx2=x[2]*(1-x[2])*(B[1,2]+0.2*x[2]-(B[1,2]+B[2,1])*x[1])   # use 0.2x[2], -0.2*x[2]
    return(list(c(dx1,dx2)))
    })
}


require(deSolve)

times=seq(0,50,length=200)
x0=c(0.4,0.4)
out=lsoda(y = x0,times = times,func = f,parms = parms)

plot(out[,2],out[,3],type="l",
     xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),
     bty="n",lwd=2)
lines(x=c(0,0,1,1,0),y=c(0,1,1,0,0),lwd=2)
points(x=c(0,0,1,1),y=c(0,1,0,1),pch=21,bg="black",cex=2)


# BISTABLE

parms=list(A=rbind(c(0,1),c(1,0)),B=rbind(c(0,1),c(1,0)))


f=function(t,x,parms){
  with(parms,{
    dx1=x[1]*(1-x[1])*(A[1,2]-(A[1,2]+A[2,1])*x[2])   #use 0.2*x[1], -0.2*x[1]
    dx2=x[2]*(1-x[2])*(B[1,2]-(B[1,2]+B[2,1])*x[1])   # use 0.2x[2], -0.2*x[2]
    return(list(c(dx1,dx2)))
  })
}


times=seq(0,1500,length=5000)
x0=c(0.015,0.01)
out=lsoda(y = x0,times = times,func = f,parms = parms)

plot(out[,2],out[,3],type="l",
     xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),
     bty="n",lwd=2)

x0=c(0.995,0.99)
out=lsoda(y = x0,times = times,func = f,parms = parms)
lines(out[,2],out[,3],type="l",lwd=2)

x0=c(0.99,0.995)
out=lsoda(y = x0,times = times,func = f,parms = parms)
lines(out[,2],out[,3],type="l",lwd=2)

x0=c(0.01,0.015)
out=lsoda(y = x0,times = times,func = f,parms = parms)
lines(out[,2],out[,3],type="l",lwd=2)

lines(x=c(0,1),y=c(0,1),lwd=2)
lines(x=c(0,1),y=c(1,0),lwd=2)
lines(x=c(0,0,1,1,0),y=c(0,1,1,0,0),lwd=2)
points(x=c(0,0,1,1),y=c(0,1,0,1),pch=21,bg="black",cex=2)
points(x=c(0.5),y=c(0.5),pch=21,bg="black",cex=2)