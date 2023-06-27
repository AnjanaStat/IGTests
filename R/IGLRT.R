#' Find test statistic value from a given data and corresponding critical value
#'
#' More detailed description
#'
#' @param data a real data frame
#' @param k a positive integer
#' @param alpha real number between 0 and 1 called significance level
#'
#' @return numeric vector
#'
#' @examples
#' k=4;N=c(20,10,10,30);mu=c(1,1,1,1);l=c(4,5,5,6)
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=mgcv::rig(N[i],mu[i],1/l[i])
#' }
#' data=g
#' IGLRT(data,k,0.05)
#' @export
IGLRT<-function(data,k,alpha)
{
  fun1<-function(data,k)
  {
    g=data
    N=unlist(rbind(lapply(data,length)))
    M=unlist(rbind(lapply(data,mean)))
    l=NULL
    for(i in 1:k)
    {
      tm1=sum(1/g[[i]]-1/M[i])/N[i]
      l[i]=1/tm1
    }
    mu0=sum(N*M)/sum(N)
    l0=l;x1=mu0;x2=M;l1=l
    repeat
    {
      u=NULL
      for(i in 1:k)
      {
        tmp1=N[i]*l0[i]
        u[i]=tmp1
      }
      neu=NULL
      for(i in 1:k)
      {
        tmp2=u[i]*M[i]
        neu[i]=tmp2
      }
      mun=sum(neu)/sum(u)
      L0=NULL
      for(i in 1:k)
      {
        tmp4=sum((g[[i]]-mun)^2/(mun^2*g[[i]]))/N[i]
        L0[i]=1/tmp4
      }
      diff=abs(mun-x1)
      if(diff<=0.00001)
      {
        break
      }
      x1=mun
      l0=L0
    }
    repeat
    {
      w=NULL
      for(i in 1:k)
      {
        tmp5=N[i]*l1[i]
        w[i]=tmp5
      }
      tmp6=Iso::pava(M,w)
      mu=tmp6
      L1=NULL
      for(i in 1:k)
      {
        tmp7=sum((g[[i]]-mu[i])^2/(g[[i]]*mu[i]^2))/N[i]
        L1[i]=1/tmp7
      }
      diff2=max(abs(mu-x2))
      if(diff2<=0.00001)
      {
        break
      }
      x2=mu;l1=L1
    }
    rto=L0/L1;ratio=NULL
    for(i in 1:k)
    {
      rslt=(rto[i])^(N[i]/2)
      ratio[i]=rslt
    }
    value=prod(ratio)
    return(value)
  }
  fun2<-function(N,mu,l,k)
  {
    g=NULL
    for(i in 1:k)
    {
      g[[i]]=mgcv::rig(N[i],mu,1/l[i])
    }
    value=fun1(g,k)
    return(value)
  }
  fun3<-function(N,mu,l,k,alpha)
  {
    x<-replicate(5000,fun2(N,mu,l,k))
    y<-sort(x,decreasing=FALSE)
    m=(alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  M=unlist(rbind(lapply(data1,mean)))
  g=data1
  l=NULL
  for(i in 1:k)
  {
    tm1=sum(1/g[[i]]-1/M[i])/N[i]
    l[i]=1/tm1
  }
  mu0=sum(N*M)/sum(N)
  set.seed(924)
  test_statistic<-fun1(data1,k)
  crit_value<-fun3(N,mu0,l,k,alpha)
  #result<-c(test_statistic,crit_value)
  return(list(Test_Statistic=test_statistic, Crit_Value=crit_value))
  #return(result)
}
