#' Main estimation function
#'
#' @param xob observed survival time
#' @param s.ob surrogate information at time t.0
#' @param deltaob event indicator
#' @param aob treatment indicator
#' @param t time at which the primary outcome is measured
#' @param t.0 time at which the surrogate is measured
#' @param varind whether to estimate variance (yes=0, no=1)
#' @param re number of replications for resampling
#'
#' @return estimated PTE, g1, and g2(s) at equally spaced s point
#' @importFrom stats pnorm rexp
#'
pte.survival = function(xob, s.ob, deltaob, aob, t, t.0, varind=0, re=100){

  n=length(xob)
  yob=as.numeric(xob>t)

  s.ob=(s.ob-mean(s.ob))/sd(s.ob)
  s.ob=pnorm(s.ob, mean = 0, sd = 1)
  nn=199
  from = 0.01; to = .99; step=((to - from)/nn)
  s=seq(from, to, by = step)

  data=cbind(xob,yob,deltaob,s.ob,aob)
  set.seed(2021)
  indexindex=sample(n, n/2, replace = FALSE)
  data1=data[indexindex,]
  data2=data[-indexindex,]
  n1=nrow(data1)
  n2=nrow(data2)

  #### pte2 given data1
  xob=data1[,1];deltaob=data1[,3];s.ob=data1[,4];aob=data1[,5];n=n1

  temp=WEIGHT(xob,deltaob,aob,n=n,t,t.0)
  wei.t0=temp[,1];wei.t=temp[,2]
  bw = 1.06*sd(s.ob)*n^(-1/5)/(n^0.11)
  kern = Kern.FUN(zz=s,zi=s.ob,bw)
  p0.t0.s.hat=apply((aob==0)*(xob>t.0)*kern*wei.t0,2,sum)/sum((aob==0)*wei.t0)
  p.t0.s.hat=apply((xob>t.0)*kern*wei.t0,2,mean)/mean(wei.t0)
  m.s.hat=(apply((xob>t)*kern*wei.t,2,sum)/sum(wei.t))/(apply((xob>t.0)*kern*wei.t0,2,sum)/sum(wei.t0))
  kern2 = Kern.FUN(zz=s.ob,zi=s.ob,bw)
  m.sob.hat=(apply((xob>t)*kern2*wei.t,2,sum)/sum(wei.t))/(apply((xob>t.0)*kern2*wei.t0,2,sum)/sum(wei.t0))
  c.hat=sum((aob==0)*(xob>t.0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)-sum((aob==0)*m.sob.hat*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
  integrand<-p0.t0.s.hat^2/p.t0.s.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  p0t=sum((aob==0)*(xob<=t.0)*wei.t0)/sum((aob==0)*wei.t0)
  pt=mean((xob<=t.0)*wei.t0)/mean(wei.t0)
  g.s.hat=m.s.hat+p0.t0.s.hat/p.t0.s.hat *c.hat/(temp+p0t^2/pt)
  g1.hat=p0t/pt *c.hat/(temp+p0t^2/pt)
  g.s.hat.1=g.s.hat
  g1.hat.1=g1.hat

  ##
  xob=data2[,1];yob=data2[,2];deltaob=data2[,3];s.ob=data2[,4];aob=data2[,5];n=n2

  temp=WEIGHT(xob,deltaob,aob,n=n,t,t.0)
  wei.t0=temp[,1];wei.t=temp[,2]
  causal=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(s.ob[kk]-s))}))
  causals=sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*aob*wei.t0)/sum(aob*wei.t0)-
    sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
  pte2=causals/causal

  ######## pte1 given data2
  xob=data2[,1];yob=data2[,2];deltaob=data2[,3];s.ob=data2[,4];aob=data2[,5];n=n2

  temp=WEIGHT(xob,deltaob,aob,n=n,t,t.0)
  wei.t0=temp[,1];wei.t=temp[,2]
  bw = 1.06*sd(s.ob)*n^(-1/5)/(n^0.11)
  kern = Kern.FUN(zz=s,zi=s.ob,bw)
  p0.t0.s.hat=apply((aob==0)*(xob>t.0)*kern*wei.t0,2,sum)/sum((aob==0)*wei.t0)
  p.t0.s.hat=apply((xob>t.0)*kern*wei.t0,2,mean)/mean(wei.t0)
  m.s.hat=(apply((xob>t)*kern*wei.t,2,sum)/sum(wei.t))/(apply((xob>t.0)*kern*wei.t0,2,sum)/sum(wei.t0))
  kern2 = Kern.FUN(zz=s.ob,zi=s.ob,bw)
  m.sob.hat=(apply((xob>t)*kern2*wei.t,2,sum)/sum(wei.t))/(apply((xob>t.0)*kern2*wei.t0,2,sum)/sum(wei.t0))
  c.hat=sum((aob==0)*(xob>t.0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)-sum((aob==0)*m.sob.hat*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
  integrand<-p0.t0.s.hat^2/p.t0.s.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  p0t=sum((aob==0)*(xob<=t.0)*wei.t0)/sum((aob==0)*wei.t0)
  pt=mean((xob<=t.0)*wei.t0)/mean(wei.t0)
  g.s.hat=m.s.hat+p0.t0.s.hat/p.t0.s.hat *c.hat/(temp+p0t^2/pt)
  g1.hat=p0t/pt *c.hat/(temp+p0t^2/pt)
  g.s.hat.2=g.s.hat
  g1.hat.2=g1.hat

  ##
  xob=data1[,1];yob=data1[,2];deltaob=data1[,3];s.ob=data1[,4];aob=data1[,5];n=n1

  temp=WEIGHT(xob,deltaob,aob,n=n,t,t.0)
  wei.t0=temp[,1];wei.t=temp[,2]
  causal=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(s.ob[kk]-s))}))
  causals=sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*aob*wei.t0)/sum(aob*wei.t0)-
    sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
  pte1=causals/causal

  pte.es=(pte1+pte2)/2
  g1.es=(g1.hat.1+g1.hat.2)/2
  gs.es=(g.s.hat.1+g.s.hat.2)/2


  if (varind==1){
    n=nrow(data)
    vv=matrix(rexp(n*re),nrow=n)
    temp=apply(vv,2,resam,t,t.0,data,data1,data2,s,nn,step)
    pte2.re=temp[1,]
    aa=which(pte2.re>1 | pte2.re<0)
    pte1.re=temp[2,]
    aa=c(aa, which(pte1.re>1 | pte1.re<0))
    aa=unique(aa)
    if (length(aa)>= re-1){
      pte.se=NA
      g1.se=NA
      gs.se=rep(NA,nn+1)
    } else {
      pte.se=sd(temp[3,!(1:re) %in% aa])
      g1.se=sd(temp[4,!(1:re) %in% aa])
      gs.se=apply( temp[-(1:4),!(1:re) %in% aa],1,sd)
    }
    output=list("pte.est"=pte.es,"pte.ese"=pte.se,"g1.est"=g1.es,"g1.ese"=g1.se,
                "sgrid"=s,"gs.est"=gs.es,"gs.ese"=gs.se)
  }else{
    output=list("pte.est"=pte.es,"g1.est"=g1.es,"sgrid"=s,"gs.est"=gs.es)
  }

}
