#' Find test statistic value from a given data and corresponding critical value
#'
#' More detailed description
#'
#' @param data a real matrix
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
#' IGAMax(data,k,0.05)
#' @export
IGAMax<-function(data,k,alpha)
{
  fun1<-function(data)
  {
    N=unlist(rbind(lapply(data,length)))
    M=unlist(rbind(lapply(data,mean)))
    S=NULL
    for(i in 1:k)
    {
      tm1=sum((data[[i]]-M[i])^2)/(N[i]-1)
      S[i]=tm1
    }
    V=NULL
    T=NULL
    for(i in 1:k-1)
    {
      V[i]=sqrt((S[i]/N[i])+(S[i+1]/N[i+1]))
      T[i]=(M[i+1]-M[i])/V[i]
    }
    value=max(T)
    return(value)
  }
  fun2<-function(N,mu0,nu,d,k,alpha)
  {
    x=NULL
    for ( i in 1:k-2)
    {
      temp=-mu0^2/nu[i+1]
      x[i]=temp
    }
    m <- diag(d^2)
    m[row(m)-col(m)== 1]<-x
    m[row(m)-col(m)==-1]<-x
    D=diag(d)
    R=solve(D)%*%m%*%solve(D)
    q=mvtnorm::qmvnorm(1-alpha,tail="lower.tail",mean=0,sigma=R)$quantile
    return(q)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  M=unlist(rbind(lapply(data1,mean)))
  N=unlist(rbind(lapply(data1,length)))
  l=NULL
  for(i in 1:k)
  {
    tm1=sum(1/data1[[i]]-1/M[i])/N[i]
    l[i]=1/tm1
  }
  Ntot=sum(N)
  mu0=sum(N*M)/Ntot
  nu=(N*l)/(mu0*Ntot)
  d=NULL
  for(i in 1:k-1)
  {
    tmp1=mu0^2*(1/nu[i]+1/nu[i+1])
    temp=sqrt(tmp1)
    d[i]=temp
  }
  set.seed(492)
  statistic_value<-fun1(data1)
  crit_value<-fun2(N,mu0,nu,d,k,alpha)
  return(list(Test_Statistic=statistic_value, Crit_Value=crit_value))
}

