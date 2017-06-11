# This R file contains the code used to create the figures
# for the following paper:
# Title: "Evolution as a coexistence mechanism: Does genetic architecture matter?" 
# Authors: Sebastian J. Schreiber, Swati Patel, and Casey ter Horst
# Journal: American Naturalist


# load the ODE solving package

library(deSolve)


#################################################
# Functions for the exploitative competition ODEs
#################################################

# Inputs for all of these functions are times t, 
# vector of densities n, and parameters parms
# Output for all of these functions are the vector dn/dt
# for the community

# Function for the haploid exploitative competition model
# n1 and n2 are the densities of the two prey genotypes
# p1 and p2 are the densities of the two predator species
# b is the intrinsic birth rate of the prey
# K the number of habitat sites for the prey
# d the per-capita death rate of the prey
# mu the mutation probability (symmetric)
# a[i,j] is the per-capita attack rate on prey i by predator j
# c[i,j] is the conversion efficience of predator j on prey i
# delta[i] is the per-capita mortality rate of predator i

haploid=function(t,n,parms){
	with(parms,{
	  n=pmax(n,0)
		p=n[3:4]
		dn1=(1-mu)*n[1]*b*(1-(n[1]+n[2])/K)+mu*n[2]*b*(1-(n[1]+n[2])/K)-d*n[1]-a[1,1]*n[1]*p[1]-a[1,2]*n[1]*p[2]
		dn2=(1-mu)*n[2]*b*(1-(n[1]+n[2])/K)+mu*n[1]*b*(1-(n[1]+n[2])/K)-d*n[2]-a[2,1]*n[2]*p[1]-a[2,2]*n[2]*p[2]
		dp1=p[1]*(c[1,1]*a[1,1]*n[1]+c[2,1]*a[2,1]*n[2]-delta[1])
		dp2=p[2]*(c[1,2]*a[1,2]*n[1]+c[2,2]*a[2,2]*n[2]-delta[2])
		dn=c(dn1,dn2,dp1,dp2)
		list(dn)
		})
}

# Function for the diploid exploitative competition model
# n1, n2, n3 are the densities of A1A1, A1A2, and A2A2 prey genotypes
# x[i] are the frequencies of the prey genotypes
# parameters as in the haploid model


diploid=function(t,n,parms){
	with(parms,{
		n=pmax(n,0)
		p=n[4:5]
		N=sum(n[1:3])
		x=n[1:3]/N
		dn1=N*b*(1-N/K)*(x[1]^2+x[1]*x[3]+x[3]^2/4)-d*n[1]-a[1,1]*n[1]*p[1]-a[1,2]*n[1]*p[2]
		dn2=N*b*(1-N/K)*(x[2]^2+x[2]*x[3]+x[3]^2/4)-d*n[2]-a[2,1]*n[2]*p[1]-a[2,2]*n[2]*p[2]
		dn3=N*b*(1-N/K)*(x[1]*x[3]+x[2]*x[3]+x[1]*x[2]*2+x[3]^2/2)-d*n[3]-a[3,1]*n[3]*p[1]-a[3,2]*n[3]*p[2]
		dp1=p[1]*(c[1,1]*a[1,1]*n[1]+c[2,1]*a[2,1]*n[2]+c[3,1]*a[3,1]*n[3]-delta[1])
		dp2=p[2]*(c[1,2]*a[1,2]*n[1]+c[2,2]*a[2,2]*n[2]+c[3,2]*a[3,2]*n[3]-delta[2])
		dn=c(dn1,dn2,dn3,dp1,dp2)
		list(dn)
		})
}

#################################################
# Two auxiliary functions
#################################################

# function to compute alpha (ecological pleiotropy)
# and beta (defensive dominance) for the diploid model
# as described in the manuscript. 

alphabeta=function(parms){
  with(parms,{
    alphas=log(c(c[2,1]/c[1,1],c[1,2]/c[2,2]))
    betas=log(c((a[2,1]-a[3,1])/(a[3,1]-a[1,1]),
                (a[1,2]-a[3,2])/(a[3,2]-a[2,2])))
    return(list(alphas=alphas,betas=betas))
  })
}

# function that computes the ratios 
# n11^1/n11^2 and n22^2/n22^1 of 
# break even densities both need to be >1 for coexistence

coexist=function(parms){
  with(parms,{
    n111=delta[1]/(c[1,1]*a[1,1])
    p111=b*(1-n111/K)/a[1,1]
    n112=delta[1]/(c[1,2]*a[1,2])
    p112=b*(1-n112/K)/a[1,2]
    n222=delta[2]/(c[2,2]*a[2,2])
    p222=b*(1-n222/K)/a[2,2]
    n221=delta[2]/(c[2,1]*a[2,1])
    p221=b*(1-n221/K)/a[2,1]
    del22=abs((a[2,1]-a[3,1])/(a[2,2]-a[3,2]))
    del11=abs((a[1,2]-a[3,2])/(a[1,1]-a[3,1]))
    return(list(ratio1=n111/n112,ratio2=n222/n221,
  ratio3=del11*del22*(p221/p222)*(p112/p111)*(n111/n112)*(n222/n221)))
  })
}

#################################################
# Functions for creating the bifurcation diagrams
#################################################


# Function for creating the ecological pleiotropy
# orbital bifurcation diagram (Figure 5A)
# Input: h which determines the attack rate on the heterozygote genotypes
# Output: range of c[1,1]=c[2,2] values for the bifurcation diagram,
#         min/max/mean frequencies of prey A1A1 genotype and predator 1
#         attack rates 

blahblah=function(h=1.5){
  blah=function(c1){
    parms$a=rbind(c(1,2),c(2,1),c(h,h))/25
    c2=0.2
    parms$c=rbind(c(c1,c2),c(c2,c1),c(c1+c2,c1+c2)/2)
    time=seq(0,6000,length=12000)
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
  
  k=60
  c1s=exp(seq(log(0.2*3),log(0.2/2),length=k))
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
              p1maxs=p1maxs,c1s=c1s,a=parms$a,c=parms$c))
}

# Function for creating the prey carrying capacity
# orbital bifurcation diagram (Figure 5C)
# Input: h which determines the conversion efficiency on heterozygote
#        prey genotypes (default value is what is used in the figure)
# Output: range of K values for the bifurcation diagram,
#         min/max/mean frequencies of prey A1A1 genotype 

blahblahK=function(h=1/2){
  blah=function(K){
    temp=c(1,2)/25
    temp2=c(2,1)/25
    temp3=c(1,1)/25*(h*1+(1-h)*2)
    parms$a=rbind(temp,temp2,temp3)
    parms$K=K
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
  
  k=25
  Ks=seq(12.505,100,length=k)
  p1mins=numeric(k)
  p1maxs=numeric(k)
  p1means=numeric(k)
  a1mins=numeric(k)
  a1maxs=numeric(k)
  a1means=numeric(k)
  for(i in 1:k){
    out=blah(Ks[i])
    a1mins[i]=out$a1min
    a1maxs[i]=out$a1max
    a1means[i]=out$a1mean
    p1mins[i]=out$p1min
    p1maxs[i]=out$p1max
    p1means[i]=out$p1mean
  }
  
  return(list(a1mins=a1mins,a1maxs=a1maxs,p1mins=p1mins,
              p1maxs=p1maxs,Ks=Ks))
}


# Function for creating the defensive dominance 
# orbital bifurcation diagram (Figure 5B)
# Input: K prey carrying capacity (default value is what is used in the figure)
# Output: range of weights h used to determine the 
#         the attack rate on the heterozygous prey
#         for the bifurcation diagram,
#         min/max/mean frequencies of prey A1A1 genotype 

blahblahh=function(K=75){
  blah=function(h){
    temp=c(1,2)/25
    temp2=c(2,1)/25
    temp3=c(1,1)/25*(h*1+(1-h)*2)
    parms$a=rbind(temp,temp2,temp3)
    parms$K=K
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
  
  k=25
  convertbeta2h=function(beta)beta/(1+beta)
  hs=seq(convertbeta2h(0.8),convertbeta2h(1.2),length=k)
  p1mins=numeric(k)
  p1maxs=numeric(k)
  p1means=numeric(k)
  a1mins=numeric(k)
  a1maxs=numeric(k)
  a1means=numeric(k)
  for(i in 1:k){
    out=blah(hs[i])
    a1mins[i]=out$a1min
    a1maxs[i]=out$a1max
    a1means[i]=out$a1mean
    p1mins[i]=out$p1min
    p1maxs[i]=out$p1max
    p1means[i]=out$p1mean
  }
  
  return(list(a1mins=a1mins,a1maxs=a1maxs,p1mins=p1mins,
              p1maxs=p1maxs,hs=hs))
}



