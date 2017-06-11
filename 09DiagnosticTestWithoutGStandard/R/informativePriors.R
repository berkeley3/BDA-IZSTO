pparameter=function(pi, w90, nprior=10, 
                    prior.method=c('median-informative','mean-informative','mode-informative')){
  
  median.inf=function(x, pi, w90){
    y=numeric(2)
    y[1]=pbeta(pi, x[1], x[2])-0.5
    y[2] = qbeta(0.95, x[1], x[2])-qbeta(0.05, x[1], x[2])-w90
    y} 
  
  mean.inf=function(x, pi, w90){
    y=numeric(2)
    y[1]=x[1]/(x[1]+ x[2])-pi
    y[2] = qbeta(0.95, x[1], x[2])-qbeta(0.05, x[1], x[2])-w90
    y} 
  
  if(prior.method=='median-informative'){
    require(nleqslv)
    x.start=c(1,1)
    param = nleqslv(x.start, median.inf, pi=pi, w90=w90)
    res=list(alpha=param$x[1], beta=param$x[2])
    #res = param$x
  }
  
  if(prior.method=='mean-informative'){
    require(nleqslv)
    x.start=c(1,1)
    param = nleqslv(x.start, mean.inf, pi=pi, w90=w90)
    res=list(alpha=param$x[1], beta=param$x[2])
    #res = param$x
  }
  if(prior.method=='mode-informative'){
    alpha=pi+1+nprior*pi
    beta=1-pi+1+nprior*(1-pi)
    res=list(alpha=alpha, beta=beta)
  }
  res
  
}

# Usage
# pparameter(pi=0.6, w90=0.3, prior.method='median-informative')
# pparameter(pi=0.6, w90=0.3, prior.method='mean-informative')




