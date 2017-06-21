RChoiceHB <- function(data,xcoding) UseMethod("RChoiceHB")
RChoiceHB.default<-function (data, xcoding) 
{
  natts = ncol(data) - 4
  npw = sum(xcoding == 0)
  ID = unique(data[, 1])
  nunits = length(ID)
  none =F
  save=T
  keep=5
  restart=F
  wgt=5
  drawdelta=F
  maxsets = max(data[, 2])
  maxalts = max(data[, 3]) + ifelse(none, 1, 0)
  nsets = rep(0, nunits)
  nalts = matrix(0, nrow = nunits, ncol = maxsets)
  for (i in 1:nunits) {
    temp = data[data[, 1] == ID[i], ]
    nsets[i] = max(temp[, 2])
    for (j in 1:nsets[i]) {
      nalts[i, j] = max(temp[temp[, 2] == j, 3]) + ifelse(none, 
                                                          1, 0)
    }
  }
  setsmap = (nalts > 0)
  altsmap = matrix(FALSE, nrow = sum(nsets), ncol = maxalts)
  index = 1
  nalts = t(nalts)
  for (i in 1:length(nalts)) {
    if (nalts[i] != 0) {
      altsmap[index, 1:nalts[i]] = TRUE
      index = index + 1
    }
  }
  avgalts = mean(nalts[nalts > 0])
  avgsets = mean(nsets)
  nalts = colSums(nalts)
  info = list(nunits = nunits, maxsets = maxsets, maxalts = maxalts, 
              nsets = nsets, nalts = nalts, setsmap = setsmap, altsmap = altsmap)
  share = FALSE
  Xt = data[, c(-1, -2, -3, -ncol(data)), drop = FALSE]
  y = data[data[, 3] == 1, ncol(data)]
  ytab = table(y)
  ytab = rbind(ytab, round(ytab/sum(ytab) * 100, 2))
  nlevels = c(rep(0, natts))
  for (i in 1:natts) {
    if (xcoding[i] == 0) {
      nlevels[i] = length(unique(Xt[, i]))
      if (any(unique(Xt[, i]) == 0)) {
        nlevels[i] = nlevels[i] - 1
      }
    }
    else {
      nlevels[i] = 1
    }
  }
  npar = sum(nlevels) - npw + ifelse(none, 1, 0)
  X = matrix(0, nrow = nrow(Xt), ncol = npar)
  Xind = 1
  nn = rep(TRUE, nrow(Xt))
  for (i in 1:natts) {
        for (j in 1:(nlevels[i] - 1)) {
        X[nn, Xind] = (Xt[nn, i] == j) * 1 - (Xt[nn, 
                                                 i] == nlevels[i]) * 1
        Xind = Xind + 1
        }
  }
  effectsmap = matrix(FALSE, nrow = natts, ncol = npar)
  index = 1
  for (i in 1:natts) {
    count = 1
    repeat {
      effectsmap[i, index] = TRUE
      index = index + 1
      count = count + 1
      if (count > max(1, nlevels[i] - 1)) {
        break
      }
    }
  }
  prior = list()
  mubar = matrix(rep(0, npar), nrow = 1)
  Amu = matrix(0.01, ncol = 1)
  df = 5
  nu = npar + df
  v = 2
  V = matrix(0, nrow = npar, ncol = npar)
  v.ind = 1
  for (i in 1:natts) {
    if (nlevels[i] == 1) {
      V[v.ind, v.ind] = 1
      v.ind = v.ind + 1
    }
    else {
      pcov = -1/nlevels[i]
      pvar = (nlevels[i] - 1)/nlevels[i]
      temp = pvar * diag(nlevels[i] - 1)
      temp[upper.tri(temp) | lower.tri(temp)] = pcov
      V[v.ind:(v.ind + nlevels[i] - 2), v.ind:(v.ind + 
                                                 nlevels[i] - 2)] = temp
      v.ind = v.ind + nlevels[i] - 1
    }
  }
  V = V * nu * v
  s = 0.1
  adjust.s = TRUE
  R = 123
  use = 5000
  constrain = FALSE
  fR = 0
  betadraw = array(0, dim = c(nunits, npar, floor(use/keep)))
  compdraw = NULL
  loglike = rep(0, floor(use/keep))
  betaout = matrix(0, nrow = nunits, ncol = npar)
  oldbetas = matrix(0, nrow = nunits, ncol = npar)
  oldll = rep(0, nunits)
  oldcomp = NULL
  muplot = matrix(0, nrow = R, ncol = npar + npw)
  xplot = (1 + fR):(R + fR)
  getLLMnl = function(beta, y, X, info) {
    nunits = info$nunits
    nalts = info$maxalts
    nsets = info$maxsets
    map.sets = t(info$setsmap)
    map.alts = t(info$altsmap)
    tmpXbeta = X * beta[rep(1:nunits, info$nalts), ]
    tmpXbeta = rowSums(tmpXbeta)
    Xbeta = matrix(0, nrow = nalts, ncol = sum(info$nsets))
    Xbeta[map.alts] = tmpXbeta
    if (is.matrix(y)) {
      Xbeta[map.alts] = exp(Xbeta[map.alts])
      denom = colSums(Xbeta)
      probs = (t(Xbeta)/denom)^y
      tmpProbs = apply(probs, 1, prod)
      probs = matrix(1, nrow = nsets, ncol = nunits)
      probs[map.sets] = tmpProbs
      ll = log(apply(probs, 2, prod))
    }
    else {
      ind = cbind(y, 1:sum(info$nsets))
      xby = matrix(0, nrow = nsets, ncol = nunits)
      xby[map.sets] = Xbeta[ind]
      Xbeta[map.alts] = exp(Xbeta[map.alts])
      denom = matrix(0, nrow = nsets, ncol = nunits)
      denom[map.sets] = log(colSums(Xbeta))
      ll = colSums(xby - denom)
    }
    return(ll)
  }
  getLndMvn = function(x, mu, rooti) {
    npar = ncol(x)
    if (is.matrix(mu)) {
      z = (x - mu) %*% rooti
    }
    else {
      z = crossprod(t(x) - mu, rooti)
    }
    logs = -(npar/2) * log(2 * pi) - 0.5 * rowSums(z * z) + 
      sum(log(diag(rooti)))
    return(logs)
  }
  drawDelta = function(x, y, comps, deltabar, Ad) {
    yy = t(t(y) - comps$mu)
    sig = tcrossprod(comps$rooti)
    xtx = crossprod(x) %x% sig
    xty = matrix(sig %*% crossprod(yy, x), ncol = 1)
    cov = chol2inv(chol(xtx + Ad))
    return(cov %*% (xty + Ad %*% deltabar) + t(chol(cov)) %*% 
             rnorm(length(deltabar)))
  }
  mnlRwMetropOnce = function(y, X, oldbeta, oldll, s, inc.root, 
                             betabar, rootpi, info, constraints, oldbeta.c) {
    stay = 0
    nunits = info$nunits
    npar = ncol(oldbeta)
    increment = s * crossprod(matrix(rnorm(nunits * npar), 
                                     ncol = nunits), inc.root)
    newbeta = oldbeta + increment
    newll = getLLMnl(newbeta, y, X, info)
    newlpost = newll + getLndMvn(newbeta, betabar, rootpi)
    ldiff = newlpost - oldll - getLndMvn(oldbeta, betabar, 
                                         rootpi)
    alpha = exp(ldiff)
    alpha[alpha > 1] = 1
    unif = runif(nunits)
    unif[alpha == 1] = 0
    good = (unif <= alpha)
    betadraw = oldbeta
    betadraw[good, ] = newbeta[good, ]
    if (!missing(constraints)) {
      betadraw.c = oldbeta.c
      betadraw.c[good, ] = newbeta.c[good, ]
    }
    oldll[good] = newll[good]
    stay = sum(!good)
    if (!missing(constraints)) {
      return(list(betadraw = betadraw, betadraw.c = betadraw.c, 
                  stay = stay, oldll = oldll))
    }
    else {
      return(list(betadraw = betadraw, stay = stay, oldll = oldll))
    }
  }
  rGibbs = function(y, betabar, A, nu, V) {
    temp = rmultireg(y, matrix(rep(1, nrow(y)), ncol = 1), 
                     betabar, A, nu, V)
    comps = list(mu = as.vector(temp$B), rooti = backsolve(chol(temp$Sigma), 
                                                           diag(ncol(temp$Sigma))))
    return(comps)
  }
  rmultireg = function(Y, X, Bbar, A, nu, V) {
    n = nrow(Y)
    m = ncol(Y)
    k = ncol(X)
    RA = chol(A)
    W = rbind(X, RA)
    Z = rbind(Y, RA %*% Bbar)
    IR = backsolve(chol(crossprod(W)), diag(k))
    Btilde = crossprod(t(IR)) %*% crossprod(W, Z)
    S = crossprod(Z - W %*% Btilde)
    rwout = rwishart(nu + n, chol2inv(chol(V + S)))
    B = Btilde + IR %*% matrix(rnorm(m * k), ncol = m) %*% 
      t(rwout$CI)
    return(list(B = B, Sigma = rwout$IW))
  }
  rwishart = function(nu, V) {
    m = nrow(V)
    df = (nu + nu - m + 1) - (nu - m + 1):nu
    if (m > 1) {
      T = diag(sqrt(rchisq(c(rep(1, m)), df)))
      T[lower.tri(T)] = rnorm((m * (m + 1)/2 - m))
    }
    else {
      T = sqrt(rchisq(1, df))
    }
    U = chol(V)
    C = t(T) %*% U
    CI = backsolve(C, diag(m))
    return(list(W = crossprod(C), IW = crossprod(t(CI)), 
                C = C, CI = CI))
  }
  getEffectsCodedParameters = function(parms, map, xcoding) {
    out = NULL
    if (nrow(map) > length(xcoding)) {
      xcoding = c(xcoding, 1)
    }
    for (i in 1:nrow(map)) {
      if (xcoding[i] == 0) {
        out = cbind(out, parms[, map[i, ], drop = FALSE], 
                    -1 * rowSums(parms[, map[i, ], drop = FALSE]))
      }
      else {
        out = cbind(out, parms[, map[i, ], drop = FALSE])
      }
    }
    return(out)
  }
  fsh = function() {
    if (Sys.info()[1] == "Windows") {
      flush.console()
    }
    return()
  }
  
  RMSas<-c()
  converged<-F
  convergeStartedAt<-0
  itime = proc.time()[3]
  acceptr.t = 0
  rep<-1
  removeLastEl<-F
  while(rep<500000){
    mgout = rGibbs(oldbetas, mubar, Amu, nu, V)
    oldcomp = mgout
    if (rep == 1) {
      oldll = getLLMnl(oldbetas, y, X, info)
    }
    rootpi = oldcomp$rooti
    inc.root = chol(chol2inv(chol(tcrossprod(rootpi))))
    
    betabar = oldcomp$mu
    metropout = mnlRwMetropOnce(y, X, oldbetas, oldll, 
                                s, inc.root, betabar, rootpi, info)
    oldbetas = metropout$betadraw
    oldll = metropout$oldll
    RLH = exp(mean(oldll))^(1/avgsets)
    PctCert = (mean(oldll) - log(1/avgalts) * avgsets)/(-log(1/avgalts) * 
                                                          avgsets)
    AvgVar = mean(diag(crossprod(inc.root)))
    RMS = sqrt(mean(oldbetas^2))
    if (rep == 1) {
      RLH.a = RLH
      PctCert.a = PctCert
      AvgVar.a = AvgVar
      RMS.a = RMS
    } else {
      RLH.a = 0.99 * RLH.a + 0.01 * RLH
      PctCert.a = 0.99 * PctCert.a + 0.01 * PctCert
      AvgVar.a = 0.99 * AvgVar.a + 0.01 * AvgVar
      RMS.a = 0.99 * RMS.a + 0.01 * RMS
    }
    acceptr = nunits - metropout$stay
    acceptr.t = acceptr.t + acceptr
    if (acceptr/nunits < 0.3) {
      s = s * 0.9
    } else if (acceptr/nunits > 0.3) {
      s = s * 1.1
    }
    acceptr = 0
    mutemp = matrix(oldcomp$mu, nrow = 1)
    if (rep%%100 == 0) {
      if (rep == 100) {
         cat("Iteration ", "Acceptance  ", "RLH    ", 
             "Pct. Cert.  ", "Avg. Var.  ", "RMS    ", "Time to End", 
             fill = TRUE)
      }
       cat(sprintf("%9.0i  %-5.3f        %-5.3f   %-5.3f        %-5.2f       %-5.2f   %-6s", 
                   rep + fR, acceptr.t/(100 * nunits), RLH.a, PctCert.a, 
                   AvgVar.a, RMS.a, 1), 
          fill = TRUE)
      fsh()
      acceptr.t = 0
      if(!converged){
        RMSas<-c(RMSas,RMS.a)
        if(removeLastEl){
          RMSas<-RMSas[-1]
          removeLastEl<-F
        }else{
          removeLastEl<-T
        }
        if(length(RMSas)>6){
          if(!converged)
          {
            if(abs(mean(RMSas[1:floor(length(RMSas)/2)])-mean(RMSas[(floor(length(RMSas)/2)+1):length(RMSas)]))<0.05){
              converged<-T
            }
          }
        }
        
      }
    }
    #if (rep > R - use) {
    if (converged) { 
      if(convergeStartedAt==0){
        convergeStartedAt<-rep
      }
      if(rep-convergeStartedAt>5000){
        rep<-9999999	
      }else{
        if(rep>convergeStartedAt){
          mkeep = (rep - convergeStartedAt)/keep
          if (rep%%keep == 0) {
            betadraw[, , mkeep] = oldbetas
            loglike[mkeep] = sum(oldll)
            compdraw[[mkeep]] = oldcomp
          }
          betaout = betaout + oldbetas
        }
      }
    }
    rep<-rep+1
  }
  betaout = betaout/use
  betawrite = cbind(matrix(ID, ncol = 1), getEffectsCodedParameters(betaout, 
                                                                    effectsmap, xcoding))
  betanames = "ID"
  for (i in 1:natts) {
    for (j in 1:nlevels[i]) {
      betanames = c(betanames, paste("A", i, "B", j, sep = ""))
    }
  }
  if (none) {
    betanames = c(betanames, "NONE")
  }
  cat("", fill = TRUE)
  if (save) {
    switch(1 + 1 * constrain + 2 * drawdelta, return(list(betadraw = betadraw, 
                                                          compdraw = compdraw, loglike = loglike)), return(list(betadraw = betadraw, 
                                                                                                                betadraw.c = betadraw.c, compdraw = compdraw, loglike = loglike)), 
           return(list(betadraw = betadraw, deltadraw = deltadraw, 
                       compdraw = compdraw, loglike = loglike)), return(list(betadraw = betadraw, 
                                                                             betadraw.c = betadraw.c, deltadraw = deltadraw, 
                                                                             compdraw = compdraw, loglike = loglike)))
  }
  else {
    return(NULL)
  }
}
