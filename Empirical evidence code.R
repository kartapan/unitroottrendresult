rm(list=ls())
setwd("C:/Users/Tapan Kar/Dropbox/Thesis/RSC-A")
df= read.csv("RBI 2000-2023 Interest rate data.csv")
df= df$Interest.rate
N= length(df)
yt= df
library(moments)
kurtosis(yt)
skewness(yt)

DF=c(rep(0,N))
wt= c(rep(0,N))
ut= c(rep(0,N))
st= c(rep(0,N))
ht= c(rep(0,N))
error= c(rep(0,N))
sumt= c(rep(0,N))
sumsqt= c(rep(0,N))
link=  c(rep(0,N))
inf= c(rep(0,N))
cal= c(rep(0,N))
varerr= c(rep(0,N))
varl= c(rep(0,N))
cor= c(rep(0,N)) 
w= c(rep(0,N)) 
v= c(rep(0,N)) 
num= c(rep(0,N)) 
den=  c(rep(0,N)) 
num1= c(rep(0,N)) 
num2= c(rep(0,N-1))
den2= c(rep(0,N-1)) 
den1 = c(rep(0,N)) 
pt=c(rep(0,N-1))
ft=c(rep(0,N))
f1t=c(rep(0,N-1))
pIV=c(rep(0,N-1))

TRDF=c(rep(0,N))
TRgt= c(rep(0,N))
TRvt= c(rep(0,N))
TRht= c(rep(0,N))
TRerror= c(rep(0,N))
TRsumt= c(rep(0,N))
TRsumsqt= c(rep(0,N))
TRlink=  c(rep(0,N))
TRinf= c(rep(0,N))
TRcal= c(rep(0,N))
TRvarerr= c(rep(0,N))
TRvarl= c(rep(0,N))
TRcor= c(rep(0,N)) 
v= c(rep(0,N)) 
TRpt=c(rep(0,N-1))
TRft=c(rep(0,N))
TRf1t=c(rep(0,N-1))
TRpIV=c(rep(0,N-1))

TRnum= c(rep(0,N)) 
TRden=  c(rep(0,N)) 
TRnum1= c(rep(0,N)) 
TRnum2= c(rep(0,N-1))
TRden2= c(rep(0,N-1)) 
TRden1 = c(rep(0,N))


yt= df

# computing t-drift test


DF[1]=0
for(i in 2:N){
  DF[i]= yt[i]-yt[i-1]
}
ybar= mean(yt)
pt=DF[-1]
ft[1]=0
for(i in 2:N){
  ft[i]=yt[i-1]
}
f1t=ft[-1]
reg=  lm(pt~f1t)
resid= reg$residuals

f= function(x){ 
  del= x[1]
  v = x[2]
  lamda = x[3]
  lo = lamda*(exp((del^-2)/2))*sinh(v/del)
  
  
  
  
  
  for(i in 1:N-1){
    den2[i]= (((lamda)^2)+(((resid[i]-lo)^2)))
    num2[i]= (resid[i]-lo)
    pIV[i] = (1/2)*log(den2[i])+((1/2)*((v+(del*log((num2[i]+sqrt(den2[i]))/lamda)))^2))-log(del)+1/2*log(2*pi)
  }   
  mle = sum(pIV[1:N-1])  
  return(mle)
}  
A = constrOptim(c(1,0,1),f,grad= NULL,method =  "Nelder-Mead",ui= rbind(c(1,0,0),c(0,0,1)),ci=c(0,0))  
D = A$par
del1 = D[1]
v1 = D[2]
lamda1= D[3]
lo1= lamda1*(exp((del1^(-2))/2))*sinh(v1/del1)
ybar= mean(yt)
ht[1]=0
num1[1]=0
den1[1]=0
fr= function(b){
  
  for(i in 2:N){
    num1[i]= (yt[i]-yt[i-1]-lo1-b)
    den1[i]= (((lamda1)^2)+((yt[i]-yt[i-1]-lo1-b)^2))
    ht[i]= (1/2)*log(den1[i])+((1/2)*((v1+(del1*log((num1[i]+sqrt(den1[i]))/lamda1)))^2))
  }  
  mle= sum(ht[2:N]) 
  
  return(mle)
  
}
ans1= optimize(fr,interval = c(-100,100))
beta= ans1$minimum
error[1]=0
for(i in 2:N){
  error[i]= yt[i]-yt[i-1]-beta
}
meanerr= sum(error[1:N])/N
link[1]=0
num[1]=0
den[1]=0
for(i in 2:N){ 
  num[i]= (error[i]-lo1)
  den[i]= (((lamda1)^2)+((error[i]-lo1)^2))
  
  link[i]= (num[i]/(den[i]))+(((v1+del1*log((num[i]+(sqrt(den[i])))/lamda1)))*((del1)/(sqrt(den[i]))))
}
meanlink= sum(link[1:N])/N 
cal[1]=0 
for(i in 2:N){
  cal[i]= (yt[i-1]-ybar)*(link[i])
}
exp= sum(cal[1:N])
sumsqt[1]=0
for(i in 2:N){
  sumsqt[i]= ((yt[i-1]-ybar)^2)
}
exp1= sum(sumsqt[1:N])
test= (exp)/(sqrt(exp1))


varerr[1]=0
for(i in 2:N){
  varerr[i]= (error[i]-meanerr)^2
}
var= sum(varerr[1:N])/N
sderr= sqrt(var)
varl[1]=0
for(i in 2:N){
  varl[i]= (link[i]-meanlink)^2
}
varlink= sum(varl[1:N])/N
sdlink= sqrt(varlink)
cor[1]=0
for(i in 2:N){
  cor[i]= (error[i]-meanerr)*(link[i]-meanlink)
}
cov= sum(cor[1:N])/N
co= cov/(sderr*sdlink)
ratio= sdlink/sderr 

# critical value for t-drift test
cri= (-2.853)+(1- co)*0.780+ ((1- co)^2)*0.548+ ((1- co)^3)*(-0.118) 

# test statistic t-drift test
test1= (test)/sdlink 

# scaling constant estimation

lamda_drift= 1.066	+0.040*(1-co)-0.032*(1-co)^2+	0.047*(1-co)^3

# U strategy value
u_test= lamda_drift*cri

# test1< u_test, so U strategy is able to reject

# computing t- trend test

TRDF[1]=0
for(i in 2:N){
  TRDF[i]= yt[i]-yt[i-1]
} 
TRft[1]=0
for(i in 2:N){
  TRft[i]=yt[i-1]
}
TRf1t=TRft[-1]

TRgt[1]=0
for(i in 2:N){
  TRgt[i]= (i-1)
}
TRg1t=TRgt[-1]
ybar= mean(yt)
TRpt=TRDF[-1]
TRvt[1]=0
for(i in 2:N){
  TRvt[i]= i
}
TRv1t=TRvt[-1]


TRreg=  lm(TRpt~TRg1t+TRf1t)
TRresid= TRreg$residuals
TRf= function(x){ 
  del= x[1]
  v = x[2]
  lamda = x[3]
  lo = lamda*(exp((del^-2)/2))*sinh(v/del)
  
  
  
  for(i in 1:N-1){
    TRden2[i]= (((lamda)^2)+(((TRresid[i]-lo)^2)))
    TRnum2[i]= (TRresid[i]-lo)
    TRpIV[i] = (1/2)*log(TRden2[i])+((1/2)*((v+(del*log((TRnum2[i]+sqrt(TRden2[i]))/lamda)))^2))-log(del)+1/2*log(2*pi)
  }   
  mle = sum(TRpIV[1:N-1])  
  return(mle)
}  
TRA = constrOptim(c(1,0,1),TRf,grad= NULL,method =  "Nelder-Mead",ui= rbind(c(1,0,0),c(0,0,1)),ci=c(0,0))  
TRD = TRA$par
TRdel1 = TRD[1]
TRv1 = TRD[2]
TRlamda1= TRD[3]
TRlo1= TRlamda1*(exp((TRdel1^(-2))/2))*sinh(TRv1/TRdel1)
ybar= mean(yt)
TRht[1]=0
TRnum1[1]=0
TRden1[1]=0
fr= function(q){
  b= q[1]
  c= q[2]
  
  for(i in 2:N){
    TRnum1[i]= (yt[i]-yt[i-1]-TRlo1-b-(c*i/N))
    TRden1[i]= (((TRlamda1)^2)+((yt[i]-yt[i-1]-TRlo1-b-(c*i/N))^2))
    TRht[i]= (1/2)*log(TRden1[i])+((1/2)*((TRv1+(TRdel1*log((TRnum1[i]+sqrt(TRden1[i]))/TRlamda1)))^2))
  }  
  mle= sum(TRht[2:N]) 
  
  return(mle)
  
}

ans1= optim(par=c(1,1),fn = fr,lower = c(-Inf,-Inf),upper = c(Inf,Inf),method = "L-BFGS-B")
bet = ans1$par
beta0= bet[1]
beta1= bet[2]
TRerror[1]=0
for(i in 2:N){
  TRerror[i]= yt[i]-yt[i-1]-beta0- (beta1*i/N)
}
TRmeanerr= sum(TRerror[1:N])/N
TRlink[1]=0
TRnum[1]=0
TRden[1]=0
for(i in 2:N){ 
  TRnum[i]= (TRerror[i]-TRlo1)
  TRden[i]= (((TRlamda1)^2)+((TRerror[i]-TRlo1)^2))
  
  TRlink[i]= (TRnum[i]/(TRden[i]))+(((TRv1+TRdel1*log((TRnum[i]+(sqrt(TRden[i])))/TRlamda1)))*((TRdel1)/(sqrt(TRden[i]))))
}
TRmeanlink= sum(TRlink[1:N])/N 
TRlink1= TRlink[-1]
TRmeanlink= sum(TRlink[1:N])/N
reg1= lm(TRf1t~TRv1t)
res1= reg1$residuals
for(i in 1:N-1){
  TRcal[i]= res1[i]*(TRlink1[i])
}
TRexp= sum(TRcal[1:N])

for(i in 1:N-1){
  TRsumsqt[i]= (res1[i]^2)
}
TRexp1= sum(TRsumsqt[1:N])
TRtest= (TRexp)/(sqrt(TRexp1))


TRvarerr[1]=0
for(i in 2:N){
  TRvarerr[i]= (TRerror[i]-TRmeanerr)^2
}
TRvar= sum(TRvarerr[1:N])/N
TRsderr= sqrt(TRvar)
TRvarl[1]=0
for(i in 2:N){
  TRvarl[i]= (TRlink[i]-TRmeanlink)^2
}
TRvarlink= sum(TRvarl[1:N])/N
TRsdlink= sqrt(TRvarlink)
TRcor[1]=0
for(i in 2:N){
  TRcor[i]= (TRerror[i]-TRmeanerr)*(TRlink[i]-TRmeanlink)
}
TRcov= sum(TRcor[1:N])/N
TRco= TRcov/(TRsderr*TRsdlink)

# critical value for t-trend test
TRcri= (-3.389)+(1- TRco)*(1.213)+ ((1- TRco)^2)*(0.658)+ ((1- TRco)^3)*(-0.126) 

# test statistic t-trend test
TRtest1= (TRtest)/TRsdlink 

# Trtest1 > critical value, so t-trend is not able to reject


# H method by Harvey

library(urca)
L0adf= ur.ers(yt,type = "DF-GLS",model ="constant",lag.max=0)

# test statistic drift version ERS test
erstest= L0adf@teststat

# critical value
erscri= L0adf@cval[,2]

# scaling constant

scaers= 1.095*erscri

# erstest> scaers, so drift version is not able to reject

trL0adf= ur.ers(yt,type = "DF-GLS",model ="trend",lag.max=0)

# test statistic ERS test trend version
trerstest= trL0adf@teststat

#critical value
trerscri= trL0adf@cval[,2]

# scaling constant

trscaers= 1.095*trerscri

# trerstest > trerscri, so trend version not able to reject, Harvey strategy not able to reject



