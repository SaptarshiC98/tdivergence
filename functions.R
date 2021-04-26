wt_t_dist=function(x,y,w){
  return(sum(w*(x-y)*atan(x-y)))
}
vec_t_wt=function(x,y,w){
  return((x-y)*atan(x-y))
}

Update_mu_1=function(x,tmax=20){
  mu=median(x)
  for(t in 1:tmax){
    f1=sum(atan(mu-x)+(mu-x)/(1+(mu-x)^2))
    f2=sum(2/(1+(mu-x)^2))
    mu=mu-f1/f2
  }
  return(mu)
}

Update_mu=function(X,tmax=20){
  p=dim(X)[2]
  mu=numeric(p)
  for(l in 1:p){
    mu[l]=Update_mu_1(X[,l],tmax)
  }
  return(mu)
}

t_wkmeans=function(X,M,beta=4,tmax=30){
  
  if(is.vector(M)==TRUE){
    M=as.matrix(M)
    M=t(M)
  }
  
  n=dim(X)[1]
  d=dim(X)[2]
  c=dim(M)[1]
  weight=rep(1/d,d)
  label=numeric(n)
  dist=numeric(c)
  t=0
  D=numeric(d)
  #update membership
  repeat{
    t=t+1
    
    for(i in 1 : n){
      for(j in 1 : c){
        dist[j]=wt_t_dist(X[i,],M[j,],weight^beta)
      }
      label[i]=which.min(dist)
    }
    
    #update centres
    for(i in 1:c){
      I=which(label==i)
      M[i,]=Update_mu(X[I,])#colMeans(X[I,])
    }
    
    #update weights
    for(j in 1:d){
      D[j]=0
    }
    for(i in 1:c){
      I=which(label==i)
      for(k in I){
        D=D+vec_t_wt(X[k,],M[i,])
      }
    }
    
    for(i in 1:d){
      if(D[i]!=0){
        D[i]=1/D[i]
        D[i]=D[i]^(1/(beta-1))
      }
    }
    sum=sum(D)
    weight=D/sum
    if(t>tmax){
      break
    }
    
  }
  return(list(label,M,weight))
  
}
