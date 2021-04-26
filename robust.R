library(prevtoinc)
cost_binom=function(theta,x,N){
  G=ecdf(x)
  s=0
  for(i in 0:N){
    a=G(i)-G(i-0.5)-dbinom(i,N,theta)
    s=s+a*atan(a)
  }
  return(s)
}
y=rbinom(100,50,0.7)
plot(table(y))
cost_binom(0.45,y,50)
s=seq(0.1,0.9,length.out = 50)
val=numeric(50)
for(i in 1:50){
  val[i]=cost_binom(s[i],y,50)
}
plot(s,val,ty='l')
est_atan_binom=function(x,N){
  f=function(theta){
    G=ecdf(x)
    s=0
    for(i in 0:N){
      a=G(i)-G(i-0.2)-dbinom(i,N,theta)
      s=s+a*atan(a)
    }
    return(s)
  }
  o=optimize(f,lower=0.01,upper=0.99)
  return(o$minimum)
}

est_hellinger=function(x,N){
  f=function(theta){
    G=ecdf(x)
    s=0
    for(i in 0:N){
      a=G(i)-G(i-0.2)
      b=dbinom(i,N,theta)
      s=s+(sqrt(a)-sqrt(b))^2
    }
    return(s)
  }
  o=optimize(f,lower=0.01,upper=0.99)
  return(o$minimum)
}

est_TV=function(x,N){
  f=function(theta){
    G=ecdf(x)
    s=0
    for(i in 0:N){
      a=G(i)-G(i-0.2)
      b=dbinom(i,N,theta)
      s=s+abs(a-b)
    }
    return(s)
  }
  o=optimize(f,lower=0.01,upper=0.99)
  return(o$minimum)
}

est_atan_binom(y,50)
simu_atan=numeric(100)
simu_l1=numeric(100)
simu_kl=numeric(100)
simu_hellinger=numeric(100)
simu_TV=numeric(100)
for(i in 1:100){
  x=rbinom(60,50,0.5)
  x=c(x,sample(40:50,40,replace = TRUE)-1)
  simu_atan[i]=est_atan_binom(x,50)
  simu_kl[i]=mean(x)/50
  simu_l1[i]=median(x)/50
  simu_hellinger[i]=est_hellinger(x,50)
}
hist(simu_atan,freq = FALSE,col=4,lwd=2,las=1)
hist(simu_kl,freq = FALSE,col=4)
hist(simu_hellinger)
hist(simu_TV)
sd(simu_atan)
mean(simu_atan)
mean(simu_hellinger)
sd(simu_hellinger)
sqrt(mean((simu_atan-0.5)^2))
sqrt(mean((simu_l1-0.5)^2))
sqrt(mean((simu_hellinger-0.5)^2))
outlier_per=seq(0,49,by=5)
TV_se=hel_se=l1_se=kl_se=atan_se=numeric(10)
TV_mean=hel_mean=l1_mean=kl_mean=atan_mean=numeric(10)
for(j in 1:10){
  for(i in 1:100){
    x=rbinom(100-outlier_per[j],50,0.5)
    x=c(x,sample(40:50,outlier_per[j],replace = TRUE)-1)
    simu_atan[i]=est_atan_binom(x,50)
    simu_TV[i]=est_TV(x,50)
    simu_kl[i]=mean(x)/50
    simu_l1[i]=median(x)/50
    simu_hellinger[i]=est_hellinger(x,50)
  }
  atan_se[j]=sqrt(mean((simu_atan-0.5)^2))
  l1_se[j]=sqrt(mean((simu_l1-0.5)^2))
  hel_se[j]=sqrt(mean((simu_hellinger-0.5)^2))
  kl_se[j]=sqrt(mean((simu_kl-0.5)^2))
  TV_se[j]=sqrt(mean((simu_TV-0.5)^2))

  ME[j,2]=mean(simu_kl)
  ME[j,3]=mean(simu_l1)
  ME[j,4]=mean(simu_hellinger)
  ME[j,5]=mean(simu_TV)
  ME[j,6]=mean(simu_atan)

  cat(j)
  cat('\n')
}
plot(outlier_per,atan_se,ty='l',lwd=2,ylim=c(0,0.1))
points(outlier_per,atan_se,ty='l',col=2,lwd=2)
points(outlier_per,TV_se,ty='l',col=3,lwd=2)
points(outlier_per,l1_se,ty='l',lwd=2,col=4)
points(outlier_per,hel_se,ty='l',col=6,lwd=2)
SE=matrix(0,10,6)
ME=matrix(0,10,6)
SE[,1]=outlier_per
ME[,1]=outlier_per
SE[,2]=kl_se
SE[,3]=l1_se
SE[,4]=hel_se
SE[,5]=TV_se
SE[,6]=atan_se
write.matrix(SE,'SE.csv',sep=',')
write.matrix(ME,'ME.csv',sep=',')
