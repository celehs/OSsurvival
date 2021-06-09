
# Internal functions
#' @importFrom stats dnorm
Kern.FUN <- function(zz, zi, bw)
{
  out = (VTM(zz,length(zi))- zi)/bw
  dnorm(out)/bw
}

VTM <- function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

WEIGHT <- function(xob, deltaob, aob, n, t, t.0){
  x = xob[aob==0]; delta = 1-deltaob[aob==0]
  xsort = sort(x[delta==1])
  risk = VTM(x, length(xsort))>=xsort
  risk.n = apply(risk, 1, sum)
  Lam = cumsum(1/risk.n)
  s0 = exp(-Lam)
  sur0 = data.frame("time"=xsort, "surv"=s0)
  x = xob[aob==1]; delta = 1-deltaob[aob==1]
  xsort = sort(x[delta==1])
  risk = VTM(x, length(xsort))>=xsort
  risk.n = apply(risk, 1, sum)
  Lam = cumsum(1/risk.n)
  s1 = exp(-Lam)
  sur1 = data.frame("time"=xsort, "surv"=s1)
  ind0 = c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0 = sur0$surv[ind0]
  ind1 = c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1 = sur1$surv[ind1]
  G = (1-aob)*G0+aob*G1
  G.t0 = (1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t = (1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0 = (xob<=t.0)*deltaob/G+(xob>t.0)/G.t0; #wei.t0[is.nan(wei.t0)]=max(wei.t0[!is.nan(wei.t0)])
  wei.t = (xob<=t)*deltaob/G+(xob>t)/G.t; #wei.t[is.nan(wei.t)]=max(wei.t[!is.nan(wei.t)])

  out = cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}

WEIGHT.p <- function(xob, deltaob, aob, n, v, t, t.0){
  x = xob[aob==0];delta = 1-deltaob[aob==0]
  xsort = sort(x[delta==1])
  risk = VTM(x, length(xsort))>=xsort
  risk.n = apply(risk*VTM(v[aob==0],length(xsort)), 1, sum)
  N = (VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN = rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu = apply(VTM(v[aob==0],length(xsort))*dN, 1, sum)
  Lam = cumsum(nu/risk.n)
  s0 = exp(-Lam)
  sur0 = data.frame("time"=xsort, "surv"=s0)
  x = xob[aob==1]; delta = 1-deltaob[aob==1]
  xsort = sort(x[delta==1])
  risk = VTM(x,length(xsort))>=xsort
  risk.n = apply(risk*VTM(v[aob==1],length(xsort)), 1, sum)
  N = (VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN = rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu = apply(dN*VTM(v[aob==1],length(xsort)), 1, sum)
  Lam = cumsum(nu/risk.n)
  s1 = exp(-Lam)
  sur1 = data.frame("time"=xsort,"surv"=s1)
  ind0 = c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0 = sur0$surv[ind0]
  ind1 = c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1 = sur1$surv[ind1]
  G = (1-aob)*G0+aob*G1
  G.t0 = (1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t = (1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0 = (xob<=t.0)*deltaob/G+(xob>t.0)/G.t0; #wei.t0[is.nan(wei.t0)]=max(wei.t0[!is.nan(wei.t0)])
  wei.t = (xob<=t)*deltaob/G+(xob>t)/G.t; #wei.t[is.nan(wei.t)]=max(wei.t[!is.nan(wei.t)])

  out = cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}

#' @importFrom stats sd
resam <- function(vv, t, t.0, data, data1, data2, s, nn, step){
  n1 = nrow(data1); n2 = nrow(data2)

  ################ pte2
  xob = data1[,1]; yob = data1[,2]; deltaob=data1[,3]; s.ob = data1[,4]; aob = data1[,5]; n = n1
  v = vv[1:nrow(data1)]

  temp = WEIGHT.p(xob, deltaob, aob, n=n, v, t, t.0)
  wei.t0 = temp[, 1]; wei.t = temp[, 2]
  bw = 1.06*sd(s.ob)*n^(-1/5)/(n^0.11)
  kern = Kern.FUN(zz=s, zi=s.ob, bw)
  p0.t0.s.hat = apply(as.numeric(v)*(aob==0)*(xob>t.0)*kern*wei.t0,2,sum)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p.t0.s.hat = apply(as.numeric(v)*(xob>t.0)*kern*wei.t0,2,mean)/mean(as.numeric(v)*wei.t0)
  m.s.hat = (apply(as.numeric(v)*(xob>t)*kern*wei.t,2,sum)/sum(as.numeric(v)*wei.t))/
    (apply(as.numeric(v)*(xob>t.0)*kern*wei.t0,2,sum)/sum(as.numeric(v)*wei.t0))
  kern2 = Kern.FUN(zz=s.ob, zi=s.ob, bw)
  m.sob.hat = (apply(as.numeric(v)*(xob>t)*kern2*wei.t,2,sum)/sum(as.numeric(v)*wei.t))/
    (apply(as.numeric(v)*(xob>t.0)*kern2*wei.t0,2,sum)/sum(as.numeric(v)*wei.t0))
  c.hat = sum(as.numeric(v)*(aob==0)*(xob>t.0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)-
    sum(as.numeric(v)*(aob==0)*m.sob.hat*(xob>t.0)*wei.t0)/sum(as.numeric(v)*(aob==0)*wei.t0)
  integrand <- p0.t0.s.hat^2/p.t0.s.hat
  temp = (integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  p0t = sum(as.numeric(v)*(aob==0)*(xob<=t.0)*wei.t0)/sum(as.numeric(v)*(aob==0)*wei.t0)
  pt = mean(as.numeric(v)*(xob<=t.0)*wei.t0)/mean(as.numeric(v)*wei.t0)
  g.s.hat = m.s.hat+p0.t0.s.hat/p.t0.s.hat *c.hat/(temp+p0t^2/pt)
  g1.hat = p0t/pt *c.hat/(temp+p0t^2/pt)
  g.s.hat.1 = g.s.hat
  g1.hat.1 = g1.hat

  ####
  xob = data2[, 1]; yob = data2[, 2]; deltaob = data2[, 3]; s.ob = data2[, 4]; aob = data2[, 5]
  n = n2; v = vv[(nrow(data1)+1):nrow(data)]

  temp = WEIGHT.p(xob, deltaob, aob, n=n, v, t, t.0)
  wei.t0 = temp[, 1]; wei.t = temp[, 2]
  causal = sum(v*(xob>t)*aob*wei.t)/sum(v*aob*wei.t)-sum(v*(xob>t)*(1-aob)*wei.t)/sum(v*(1-aob)*wei.t)
  tempind = c(sapply(1:n, function(kk){which.min(abs(s.ob[kk]-s))}))
  causals = sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*v*aob*wei.t0)/sum(v*aob*wei.t0)-
    sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*(1-aob)*v*wei.t0)/sum(v*(1-aob)*wei.t0)
  pte2 = causals/causal

  ################ pte1
  xob = data2[, 1]; yob = data2[, 2]; deltaob = data2[, 3]; s.ob = data2[, 4]; aob = data2[, 5]
  n = n2; v = vv[(nrow(data1)+1):nrow(data)]

  temp = WEIGHT.p(xob, deltaob, aob, n=n, v, t, t.0)
  wei.t0 = temp[, 1]; wei.t = temp[, 2]
  bw = 1.06*sd(s.ob)*n^(-1/5)/(n^0.11)
  kern = Kern.FUN(zz=s, zi=s.ob, bw)
  p0.t0.s.hat = apply(as.numeric(v)*(aob==0)*(xob>t.0)*kern*wei.t0,2,sum)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p.t0.s.hat = apply(as.numeric(v)*(xob>t.0)*kern*wei.t0,2,mean)/mean(as.numeric(v)*wei.t0)
  m.s.hat = (apply(as.numeric(v)*(xob>t)*kern*wei.t,2,sum)/sum(as.numeric(v)*wei.t))/
    (apply(as.numeric(v)*(xob>t.0)*kern*wei.t0,2,sum)/sum(as.numeric(v)*wei.t0))
  kern2 = Kern.FUN(zz=s.ob,zi=s.ob,bw)
  m.sob.hat = (apply(as.numeric(v)*(xob>t)*kern2*wei.t,2,sum)/sum(as.numeric(v)*wei.t))/
    (apply(as.numeric(v)*(xob>t.0)*kern2*wei.t0,2,sum)/sum(as.numeric(v)*wei.t0))
  c.hat = sum(as.numeric(v)*(aob==0)*(xob>t.0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)-
    sum(as.numeric(v)*(aob==0)*m.sob.hat*(xob>t.0)*wei.t0)/sum(as.numeric(v)*(aob==0)*wei.t0)
  integrand<-p0.t0.s.hat^2/p.t0.s.hat
  temp = (integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  p0t = sum(as.numeric(v)*(aob==0)*(xob<=t.0)*wei.t0)/sum(as.numeric(v)*(aob==0)*wei.t0)
  pt = mean(as.numeric(v)*(xob<=t.0)*wei.t0)/mean(as.numeric(v)*wei.t0)
  g.s.hat = m.s.hat+p0.t0.s.hat/p.t0.s.hat *c.hat/(temp+p0t^2/pt)
  g1.hat = p0t/pt *c.hat/(temp+p0t^2/pt)
  g.s.hat.2 = g.s.hat
  g1.hat.2 = g1.hat

  ####
  xob = data1[, 1]; yob = data1[, 2]; deltaob = data1[, 3]; s.ob = data1[, 4]; aob = data1[, 5]
  n = n1; v = vv[1:nrow(data1)]

  temp = WEIGHT.p(xob, deltaob, aob, n=n, v, t, t.0)
  wei.t0 = temp[, 1]; wei.t = temp[, 2]
  causal = sum(v*(xob>t)*aob*wei.t)/sum(v*aob*wei.t)-sum(v*(xob>t)*(1-aob)*wei.t)/sum(v*(1-aob)*wei.t)
  tempind = c(sapply(1:n, function(kk){which.min(abs(s.ob[kk]-s))}))
  causals = sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*v*aob*wei.t0)/sum(v*aob*wei.t0)-
    sum(((xob<=t.0)*g1.hat+(xob>t.0)*g.s.hat[tempind])*(1-aob)*v*wei.t0)/sum(v*(1-aob)*wei.t0)
  pte1 = causals/causal

  pte = (pte1+pte2)/2
  g1.es = (g1.hat.1+g1.hat.2)/2
  gs.es = (g.s.hat.1+g.s.hat.2)/2

  out = c(pte2, pte1, pte, g1.es, gs.es)
}
