
rm(list=ls())

library(fda.usc)
library(PASWR)
library(PairedData)



# Under the null hypothesis

delta=1.20                        # Modify delta under the alternative
pvalue=NULL
tsim=1000                         #number of repetitions
nv=120                            # Sample size (number of curves generated) at each simulation




par(mfrow=c(2,2))
for (m in 1:tsim)
     {
       N=1000                      # Number of simulations
       nt= 100                     # Number of time periods
       t=seq(0, 10, length=nt)     # Time interval
                   
--------------------------------------
# Generating a Brownian motion $\nu(t)$
#--------------------------------------

#z=rnorm (nt, 0, 1)

length(nt)  
v=rep (0, (nt+1))
for (i in 2:(nt+1))
     {
      v[i]=v[i-1]+rnorm(1,0,0.5)
     }
v=v[-1]
plot(t,v, type="l", ylab=expression(v(t)), xlab="t") 



#--------------------------------------
# Generating curves treatment 1
#--------------------------------------

x1=matrix(0, nrow=length(t), ncol=nv)
mu=sin(2*pi*t)
for (j in 1:nv)
     {
      error= runif(nt, -1, 1)
    # error= rbeta(nt, 2, 2)
      x1[,j]=mu+error     
      }
matplot(t, x1, col=2, lty=1,type="l", ylab=expression('X'[i1](t)), xlab=expression(t), ylim=c(-2,3))


#--------------------------------------
# Generating curves treatment 2
#--------------------------------------

x2=matrix(0, nrow=length(t), ncol=nv)
mu=sin(2*pi*t)
for (j in 1:nv)
     {
      error= runif(nt, -1, 1)
    # error= rbeta(nt, 2, 2)
      x2[,j]=mu+error
      }
matplot(t, x2, col=3, lty=1,type="l", ylab=expression('X'[i2](t)), xlab=expression(t),ylim=c(-2,3))


#--------------------------------------
# Generating curves treatment 3
#--------------------------------------

x3=matrix(0, nrow=length(t), ncol=nv)
mu=sin(2*pi*t)
for (j in 1:nv)
     {
      error= runif(nt, -1, 1)
    # error=rbeta(nt, 2, 2)
      x3[,j]=mu+error+delta
      }
matplot(t, x3, col=4, lty=1,type="l", ylab=expression('X'[i3](t)), xlab=expression(t),ylim=c(-2,3))


# Random projections treatment 1

dim(x1)
length(v)
fd1<-x1*v
fdataobj1<-fdata(t(fd1),t)
rp1<-int.simpson(fdataobj1)


# Random projections treatment 2

dim(x2)
length(v)
fd2<-x2*v
fdataobj2<-fdata(t(fd2),t)
rp2<-int.simpson(fdataobj2)


# Random projections treatment 3

dim(x3)
length(v)
fd3<-x3*v
fdataobj3<-fdata(t(fd3),t)
rp3<-int.simpson(fdataobj3)

# Kruskal-Wallis

y=c(rp1,rp2,rp3)
trat=as.factor(c(rep(1,nv),rep(2,nv),rep(3,nv)))
data=as.data.frame(cbind(trat, y))
kw=kruskal.test(y ~ trat, data = data)
pvalue[m]= kw[3][[1]]

}

pvalue
power=sum(ifelse(pvalue<0.05,1,0))/tsim
power



#-------------------
# Results 
#-------------------

n=10
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,  0.22, 0.24, 0.28, 0.32, 0.38, 0.44, 0.60, 0.70,  0.90, 1.20) 
power1= c(0.046, 0.071, 0.147, 0.31, 0.49, 0.68, 0.75, 0.79, 0.85, 0.86,  0.89,  0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97,  0.98, 0.99)
plot(delta, power1,  xlab=expression(delta), ylab="Power", type="l", lwd=2, col=1)
abline(h=0.05, lty=2)


n=30
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,  0.22, 0.24, 0.28, 0.32,  0.38, 0.44, 0.60, 0.70, 0.90, 1.20) 
power2= c(0.04,  0.13,  0.45,  0.72, 0.81, 0.85, 0.89, 0.91, 0.92, 0.93,  0.94,  0.95, 0.96, 0.97, 0.973, 0.976, 0.979, 0.985, 0.993, 0.995, 0.999)
lines(delta, power2, type="l", lwd=2, col=2)


#n=50
#delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,   0.20,   0.22,   0.24,   0.28,  0.32,  0.38,   0.44,   0.60,  0.70,  0.90,  1.20) 
#power3= c(0.04,  0.20,  0.63,  0.80, 0.86, 0.88, 0.90, 0.91, 0.93, 0.94,  0.953,   0.956,  0.958,  0.96, 0.963, 0.965,   0.97,  0.975, 0.983,  0.99,  0.997)
#lines(delta, power3, type="l", lwd=2, col=4)


n=80
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,   0.22,  0.24,    0.28,  0.32,   0.38,   0.44, 0.60, 0.70, 0.90,   1.20) 
power4= c(0.04,  0.31,  0.74,  0.87, 0.91, 0.93, 0.93, 0.95, 0.96, 0.965, 0.968, 0.971, 0.974,   0.976,  0.979,  0.981, 0.983, 0.986, 0.989, 0.992,  0.995)
lines(delta, power4, type="l", lwd=2, col=4)


n=120
delta = c(0.00, 0.02,   0.04, 0.06, 0.08, 0.10, 0.12, 0.14,  0.16,  0.18,  0.20,  0.22,   0.24,  0.28, 20.32,  0.38,  0.44,   0.60,   0.70,   0.90,   1.20) 
power5= c(0.05, 0.46,   0.83, 0.89, 0.92, 0.93, 0.95, 0.96,  0.962, 0.967, 0.97,  0.972,  0.976, 0.98, 0.982,  0.985, 0.987,  0.988,  0.992, 0.992,  0.995)                  
lines(delta, power5, type="l", lwd=2, col=3)



#-------------------
# Plot paper
#-------------------

n=10
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,  0.22, 0.24, 0.28, 0.32, 0.38, 0.44, 0.60, 0.70) 
power1= c(0.046, 0.071, 0.147, 0.31, 0.49, 0.68, 0.75, 0.79, 0.85, 0.86,  0.89,  0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97)
plot(delta, power1,  xlab=expression(delta), ylab="Power", type="l", lwd=1, col=4)
abline(h=0.05, lty=3, lwd=2)


n=30
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,  0.22, 0.24,  0.28, 0.32,   0.38, 0.44,  0.60, 0.70) 
power2= c(0.04,  0.13,  0.45,  0.72, 0.81, 0.85, 0.89, 0.91, 0.92, 0.93,  0.94,  0.95, 0.96,  0.97, 0.97,   0.98, 0.98,  0.99, 0.99)
lines(delta, power2, type="l", lwd=1, col=3)


n=80
delta = c(0.00,  0.02,  0.04,  0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,  0.20,   0.22,  0.24, 0.28,  0.32,  0.38, 0.44, 0.60, 0.70) 
power4= c(0.04,  0.31,  0.74,  0.87, 0.91, 0.93, 0.94, 0.95, 0.96, 0.97,  0.97,   0.97,  0.97, 0.98,  0.98,  0.98, 0.98, 0.99, 0.99)
lines(delta, power4, type="l", lwd=1, col=2)


n=120
delta = c(0.00, 0.02,   0.04, 0.06, 0.08, 0.10, 0.12, 0.14,  0.16,  0.18,  0.20,  0.22,   0.24,  0.28, 0.32,  0.38,  0.44,   0.60,   0.70) 
power5= c(0.05, 0.46,   0.83, 0.89, 0.92, 0.93, 0.95, 0.96,  0.962, 0.967, 0.97,  0.972,  0.976, 0.98, 0.98,  0.985, 0.987,  0.988,  0.99)                  
lines(delta, power5, type="l", lwd=1, col=9)





