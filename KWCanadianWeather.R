
#--------------------------------------------------------------------
# Kruskal-Wallis test for functional data based on random projections
# Application to Canadian Weather Dataset
#--------------------------------------------------------------------


rm(list=ls())
library(fda.usc)
library(fda)
library(BSDA)
library( fdANOVA )
library(car)
library(dplyr)
library(devtools)
library(ggpubr)
# library(normtest) 
library(nortest) 
library(moments) 


#------------------------
# Canadian Weather Data
#-----------------------

data(CanadianWeather)
names(CanadianWeather)
temp = CanadianWeather$dailyAv[,,1]
daytime = (1:365)
matplot(daytime,temp,ylab="Temperature Degrees C", xlab="Day", type='l', main ="Canadian Temperature Curves")


# We can also plot by region; Atlantic, Pacific, Central and North. 
matplot(daytime, temp,  ylab="Temperature Degrees C", xlab="Day",  
        type='l', col=as.factor(CanadianWeather$region))


#---------------------------------------------------------------------
# Smoothing with a Fourier basis with harmonic acceleration operator
#---------------------------------------------------------------------
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)
daySmooth <- with(CanadianWeather, smooth.basis(day.5,
       dailyAv[,,"Temperature.C"],
       daybasis65, fdnames=list("Day", "Station", "Degrees C")) )

daySmooth2 <- smooth.fdPar(daySmooth$fd, lambda=1e5)
plot(daySmooth)
plot(daySmooth2)
class(daySmooth2)


#---------------------------
# names for (climate zones)
#-----------------------------

zonenames <- c("Atlantic", "Pacific ", "Continental", "Arctic")
temp.fdATL <- daySmooth2[1:15]
temp.fdCON <- daySmooth2[16:24]
temp.fdPAC <- daySmooth2[25:31]
temp.fdART <- daySmooth2[32:35]

#-------------------------------------------------------
# ANOVA for Functional Data with Temperature Dataset
# Note: This test is based on the normality assumption
# See Ramsay and Silverman book 
#-------------------------------------------------------

fdG1 = rep('temp.fdATL',15)
fdG2 = rep ('temp.fdCON',9)
fdG3 = rep('temp.fdPAC',7)
fdG4 = rep ('temp.fdART',4)
period <- as.factor(c(fdG1,fdG2,fdG3,fdG4))
plot(temp.fdATL, col="1")
class(daySmooth2)
fdata11 <- fdata(daySmooth2)
class(fdata11)
plot(fdata11)
res=fanova.onefactor(fdata11,period,nboot=50,plot=TRUE)
res




#--------------------------------
# Curves for each region
#--------------------------------

nt=365   
par(mfrow=c(2,2))
temp.fdata1 <- fdata(temp.fdATL, argvals = seq(1,365,length=nt))
temp.fdata2 <- fdata(temp.fdCON ,argvals = seq(1,365,length=nt))
temp.fdata3 <- fdata(temp.fdPAC,argvals = seq(1,365,length=nt))
temp.fdata4 <- fdata(temp.fdART,argvals = seq(1,365,length=nt))
plot(temp.fdata1, col=3, main=" Atlantic",ylim= c(-30,22), ylab="Temperature C", xlab=" (a) Day", lty=1)
plot(temp.fdata2, col=1, main=" Continental",ylim=c(-30,22), ylab="Temperature C", xlab=" (b) Day", lty=1)
plot(temp.fdata3, col="grey", main="Pacific", ylim= c(-30,22),ylab="Temperature C", xlab=" (c) Day")
plot(temp.fdata4, col="blue", main="Arctic",ylim= c(-35,22), ylab="Temperature C", xlab=" (d) Day")



#------------------------------------------------------
# Generating a Brownian Motion and Random Projections
#-----------------------------------------------------

                      
t=seq(0,365, length= nt)        
v=rep (0, nt)
for (i in 2:nt+1)
     {
      v[i]=v[i-1]+rnorm(1,0,0.1)
     }
v=v[-1]
plot(t,v, type="l", ylab=expression(v(t)), xlab="t")


# Plot curves $X(t)v(t)$ 

par(mfrow=c(2,2))
plot(temp.fdata1*v, main=" X1(t)*v(t)", ylab="Temp. C", xlab="Day")
plot(temp.fdata2*v, main=" X2(t)*v(t)", ylab="Temp. C", xlab="Day")
plot(temp.fdata3*v, main=" X3(t)*v(t)", ylab="Temp. C", xlab="Day")
plot(temp.fdata4*v, main=" X4(t)*v(t)", ylab="Temp. C", xlab="Day")


#-------------------------------------------------
# Integrating X(t)*v(t) for each group separately
#--------------------------------------------------

x1 <- int.simpson(temp.fdata1*v)
x2 <- int.simpson(temp.fdata2*v)  
x3 <- int.simpson(temp.fdata3*v) 
x4 <- int.simpson(temp.fdata4*v)
length(x1)
length(x2)
length(x3)
length(x4)

#-------------------------------------------------
# testing for normality with random projections
#------------------------------------------------

ad.test(x1)
ad.test(x2)
ad.test(x3)      # The test cannot be used in this case 
ad.test(x4)      # The test cannot be used in this case
lillie.test(x1)
lillie.test(x2)
lillie.test(x3)
lillie.test(x4)  # The test cannot be used in this case

# Note: The assumption of normality is not valid 
# The random projections (rp from groups 1, 2 and 3) do not follows
# a normal distribution 

#--------------------------------------------
# Kruskal-Wallis test with random projections
#--------------------------------------------

kruskal.test(list(x1,x2,x3,x4))

# Another way to write the test

projections <- c(x1,x2,x3,x4)
g <- factor(rep(1:4, c(15, 9,7,4)),labels = c("Atlantic", "Pacific ", "Continental", "Artic"))
kruskal.test(projections, g)
leveneTest(projections,g)
pairwise.wilcox.test(projections, g, p.adjust.method = "bonf")
pairwise.wilcox.test(projections, g, p.adjust.method = "BH")
pairwise.wilcox.test(projections, g, p.adjust.method = "holm")
   
   
   
#-----------------------------------------------------------------------------------------------------
# Testing the null hypothesis 1000 times (generating at each iteration a different random projection)
#-----------------------------------------------------------------------------------------------------

rep=1000
pvalue=NULL
for (j in 1:rep)
     {

     nt=365                         
     t=seq(0,365, length= nt)        
     v=rep (0, nt)
            for (i in 2:nt+1)
                {
                 v[i]=v[i-1]+rnorm(1,0,0.1)
                }
     v=v[-1]
# random_projections <- int.simpson(temp.fdata*v,  method = NULL)


#-------------------------------------------------
# Integrating X(t)*v(t) for each group separately
#--------------------------------------------------

     x1 <- int.simpson(temp.fdata1*v)
     x2 <- int.simpson(temp.fdata2*v)  
     x3 <- int.simpson(temp.fdata3*v) 
     x4 <- int.simpson(temp.fdata4*v)


#--------------------------------------------
# Kruskal-Wallis test with random projections
#--------------------------------------------

    pvalue[j]=kruskal.test(list(x1,x2,x3,x4))[[3]]
     }

boxplot(pvalue, pch=16, col=4)
alpha=0.05
p_value=sum(ifelse(pvalue<alpha,0,1))/rep
p_value 
