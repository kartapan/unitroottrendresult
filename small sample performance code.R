rm(list=ls())
N= 100 
DIM= 1000
yt= c(rep(0,N))
zt= c(rep(0,N))
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
test1=c(rep(0,DIM))
cri= c(rep(0,DIM))
count= c(rep(0,DIM))

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
TRtest1=c(rep(0,DIM))
TRcri= c(rep(0,DIM))

TRnum= c(rep(0,N)) 
TRden=  c(rep(0,N)) 
TRnum1= c(rep(0,N)) 
TRnum2= c(rep(0,N-1))
TRden2= c(rep(0,N-1)) 
TRden1 = c(rep(0,N))
rd= c(rep(0,DIM))
s= c(rep(0,DIM))
co= c(rep(0,DIM))
TRco= c(rep(0,DIM))
co1= c(rep(0,DIM))
lamda= c(rep(0,DIM))
TRlamda= c(rep(0,DIM))
rd= c(rep(0,DIM))
cou= c(rep(0,DIM))
trecount= c(rep(0,DIM))
trerstest= c(rep(0,DIM))
trerscri= c(rep(0,DIM))

trscaers= c(rep(0,DIM))
erstest= c(rep(0,DIM))
erscri= c(rep(0,DIM))

scaers= c(rep(0,DIM))
sd= c(rep(0,DIM))
Hcount= c(rep(0,DIM))
lamda_drift = c(rep(0,DIM))
lamda_trend = c(rep(0,DIM))
C= vector()

U= vector()
drift= vector()
tre= vector()
H= vector()
set= seq(from=0, to=21, by=1)
library(urca)
library(fGarch)

for(m in 1:22) {
  C1 = m-1

for(l in 1:DIM){    
  
  # data generation
  zt= rt(N,df = 3)
  zt= (zt-mean(zt))/sd(zt)
  ut[1]=zt[1]
  
  for(p in 2:N){
    ut[p]  = (1+(-C1/N))*ut[p-1]+zt[p]
  }  
  yt[1]= 1+0.025+ut[1]
  for( j in 2:N){
    
    # here we select a=0.025 and similarly other combinations are selected
    yt[j]= 1+0.025*j+ut[j]
  }
  
  # t-drift test
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
  co[l]= cov/(sderr*sdlink)
  ratio= sdlink/sderr 
  
  # t-drift test critical value
  cri[l]= (-2.853)+(1- co[l])*0.780+ ((1- co[l])^2)*0.548+ ((1- co[l])^3)*(-0.118) 
  
  # t-drift test statistic
  test1[l]= (test)/sdlink 
  
  if(test1[l]< cri[l]){
    cou[l]=1
  }
  
  # t-trend test
  
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
  TRco[l]= TRcov/(TRsderr*TRsdlink)
  
  # t-trend test critical value
  TRcri[l]= (-3.389)+(1- TRco[l])*(1.213)+ ((1- TRco[l])^2)*(0.658)+ ((1- TRco[l])^3)*(-0.126) 
  
  # t-trend test statistic
  TRtest1[l]= (TRtest)/TRsdlink 
  
  if(TRtest1[l]< TRcri[l]){
    trecount[l]=1
  }
  
  # scaling constant for U method
  lamda_drift[l]= 1.066	+0.040*(1-co[l])-0.032*(1-co[l])^2+	0.047*(1-co[l])^3
  
  # employing U strategy
  
  if(test1[l]< lamda_drift[l]*cri[l])
    rd[l]=1
  if(test1[l]>= lamda_drift[l]*cri[l])
    rd[l]=2
  if(rd[l]==1){
    count[l]=1
  }
  
  if(rd[l]==2){
    lamda_trend[l]= 1.066	+0.040*(1-TRco[l])-0.032*(1-TRco[l])^2+	0.047*(1-TRco[l])^3
    if(TRtest1[l]< lamda_trend[l]*TRcri[l]){
      count[l]=1
    }
  }
  
  # H-test by Harvey
  
    L0adf= ur.ers(yt,type = "DF-GLS",model ="constant",lag.max=0)
    erstest[l]= L0adf@teststat
    erscri= L0adf@cval[,2]
    
    # scaling constant
    
    scaers[l]= 1.095*erscri
    
    if(erstest[l]< scaers[l])
      sd[l]=1
    if(erstest[l]>= scaers[l])
      sd[l]=2
    if(sd[l]==1){
      Hcount[l]=1
    }
    
    if(sd[l]==2){
      
      trL0adf= ur.ers(yt,type = "DF-GLS",model ="trend",lag.max=0)
      trerstest[l]= trL0adf@teststat
      trerscri= trL0adf@cval[,2]
      
      trscaers[l]= 1.095*trerscri
      if(trerstest[l]< trscaers[l])
        Hcount[l]=1
      
    }
  }
  
  U[m]= sum(count)/DIM
  drift[m]= sum(cou)/DIM
  tre[m]= sum(trecount)/DIM
  H[m]= sum(Hcount)/DIM
  C[m]= -C1
}


# drift and tre denote t-drift and t-trend test

df  = data.frame(C, U, drift, tre,H)

#write.csv(df,"a1_tdist_0.025.csv")

# similarly save other results with different a values by write.csv command

df_t_0.1= read.csv("a1_tdist_0.1.csv")


df_t_0.2= read.csv("a1_tdist_0.2.csv")


df_t_0.4= read.csv("a1_tdist_0.4.csv")


df_t_0.05= read.csv("a1_tdist_0.05.csv")



df_t_0.0= read.csv("a1_tdist_0.00.csv")


df_t_0.025= read.csv("a1_tdist_0.025.csv")



df_norm_0.0= read.csv("a1_normdist_0.00.csv")


df_norm_0.1= read.csv("a1_normdist_0.1.csv")


df_norm_0.2= read.csv("a1_normdist_0.2.csv")


df_norm_0.4= read.csv("a1_normdist_0.4.csv")


df_norm_0.05= read.csv("a1_normdist_0.05.csv")

#write.csv(df,"a1_normdist_0.025.csv")
df_norm_0.025= read.csv("a1_normdist_0.025.csv")



df_sstd_0.1= read.csv("a1_sstddist_0.1.csv")


df_sstd_0.2= read.csv("a1_sstddist_0.2.csv")


df_sstd_0.4= read.csv("a1_sstddist_0.4.csv")


df_sstd_0.05= read.csv("a1_sstddist_0.05.csv")


df_sstd_0.00= read.csv("a1_sstddist_0.00.csv")


df_sstd_0.025= read.csv("a1_sstddist_0.025.csv")

