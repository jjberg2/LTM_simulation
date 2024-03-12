bgamma = function(gamma) {
  ifelse(gamma>20,1,(exp(gamma)-1)/(exp(gamma)+1))
}
normDel= function(a,thr){
  pnorm(thr) - pnorm(thr - a)
}
sigmaa = function(a,thr,nc){
  del = normDel(a,thr)
  2*a^2 * bgamma(4*nc*del) / (4*nc*del)
}
sigmaa2 = function(as,ay,thr,nc){
  del = normDel(as,thr)
  2*ay^2 * bgamma(4*nc*del) / (4*nc*del)
}
sigmaRisk = function(a,thr,nc){
  del = normDel(a,thr)
  2*del^2 * bgamma(4*nc*del) / (4*nc*del)
}
sigmaaSmall = function(a,thr,nc){
  approxDel = a*dnorm(thr)
  2*a^2 * bgamma(4*nc*approxDel) / (4*nc*approxDel)
}
sigmaaSmall2 = function(ay){
  2*ay * bgamma(ay) 
}