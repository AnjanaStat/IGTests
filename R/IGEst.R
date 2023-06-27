#' Find estimates of means and scale-like parameters when there is a prior information that means are under order restrictions
#'
#' More detailed description
#'
#' @param data a real data frame
#' @param k a positive integer
#'
#' @return numeric vector
#'
#' @examples
#' k=4;N=c(20,10,10,30);mu=c(1,1.2,1.3,1.5);l=c(4,5,5,6)
#' g=NULL
#' for(i in 1:k)
#' {
#'  g[[i]]=mgcv::rig(N[i],mu[i],1/l[i])
#' }
#' data=g
#' IGEst(data,k)
#' @export
IGEst<-function(data,k)
{
  data1<-lapply(data, function(col)col[!is.na(col)])
  g=data1
  N=unlist(rbind(lapply(data1,length)))
  M=unlist(rbind(lapply(data1,mean)))
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
  #res=c(x2,l1)
  #return(res)
  return(list(Means = x2,  scales = l1))
}
