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
#' IGMax(data,k,0.05)
#' @export
IGMax<-function(data,k,alpha)
{
  fun1<-function(data,k)
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
    m=(1-alpha)*5000
    c<-y[m]
    return(c)
  }
  data1<-lapply(data, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(data1,length)))
  M=unlist(rbind(lapply(data1,mean)))
  l=NULL
  for(i in 1:k)
  {
    tm1=sum(1/data1[[i]]-1/M[i])/N[i]
    l[i]=1/tm1
  }
  mu0=sum(N*M)/sum(N)
  set.seed(935)
  test_statistic<-fun1(data1,k)
  crit_value<-fun3(N,mu0,l,k,alpha)
  return(list(Test_Statistic=test_statistic, Crit_Value=crit_value))
}
