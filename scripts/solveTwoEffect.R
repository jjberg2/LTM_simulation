source('scripts/solveSingleEffect.R')

fourPoissonDiffs <-
  function(X,
           as = as,
           al = al,
           bt = bt,
           L = L,
           h2 = h2,
           C,
           var.ratio = var.ratio,
           equalize.observed.vars = equalize.observed.vars) {
    ## inputs
    my.deltal <- 1 / (1 + exp(X[1]))
    my.tstar <- X[2]
    my.gs <- 1 / (1 + exp(X[3]))
    my.Ll <- L * (1 - my.gs)
    my.Ls <- L * my.gs
    ##Ls <- L-Ll
    my.bs <- 1 / (1 + exp(X[4]))
    
    ##my.L  <- my.Ll + my.Ls
    
    ## unscaled small effect stuff
    my.ys <- log((1 + my.bs) / (1 - my.bs))
    my.raw.Vas <-
      8 * Ne * u * my.Ls * as ^ 2 * my.bs / my.ys # originall 4*Ne*u*Ls*as^2*bs/ys
    
    
    ## raw large effect stuff
    my.s <- my.deltal * C
    my.yl <- 4 * Ne * my.s #originally 2*Ne*my.s
    my.mean.nl <- 2 * my.Ll * u / my.s
    my.raw.Val <- al ^ 2 * my.mean.nl
    
    ## get total variance in original units
    my.raw.Vt <- (my.raw.Vas + my.raw.Val) / h2
    
    ## standardized effects
    my.std.as <- as / sqrt (my.raw.Vt)
    my.std.al <- al / sqrt (my.raw.Vt)
    
    ## standardize variances
    my.std.Vas <- my.raw.Vas / my.raw.Vt
    my.std.Val <- my.raw.Val / my.raw.Vt
    my.std.norm.sd <- sqrt(1 - h2 + my.std.Vas)
    
    ## standardized thr dens
    my.ft <- my.ys / (4 * Ne * C * my.std.as)
    
    my.range <-
      seq(qpois(1e-8, my.mean.nl), qpois(1 - 1e-8, my.mean.nl))
    my.prot <-
      pPoisConv(my.tstar,
                my.mean.nl,
                my.std.norm.sd,
                alphal = my.std.al,
                risk.allele = FALSE)
    my.risk <-
      pPoisConv(my.tstar,
                my.mean.nl,
                my.std.norm.sd,
                alphal = my.std.al,
                risk.allele = TRUE)
    
    ## global stuff
    my.meana <- (my.std.as * my.Ls + my.std.al * my.Ll) / L
    my.maxg <- 2 * my.meana * L
    my.prev <-
      as.numeric(pPoisConv(my.tstar, my.mean.nl, my.std.norm.sd, alphal = my.std.al))
    my.risk.var <- (my.prev * (1 - my.prev))
    
    ## differences
    my.deltal.tild <- my.risk - my.prot
    diff.one <- my.deltal - my.deltal.tild
    
    my.ft.tild <-
      dPoisConv(my.tstar,
                my.mean.nl,
                my.std.norm.sd,
                alphal = my.std.al,
                risk.allele = FALSE)
    diff.two <- my.ft - my.ft.tild
    
    if (equalize.observed.vars) {
      my.Vol <- my.deltal ^ 2 * my.mean.nl
      my.Vos <- my.std.Vas * my.ft ^ 2
      diff.three <- (my.Vol - var.ratio * my.Vos) / my.risk.var
    } else {
      diff.three <- my.std.Val - var.ratio * my.std.Vas
    }
    
    bt.tild <-
      (2 * my.Ls * my.bs * my.std.as + 2 * my.Ll * 1 * my.std.al) / (my.maxg)
    diff.four <- bt - bt.tild
    
    ##diff.five <- my.L*my.meana - Lmeana
    
    return(c(diff.one, diff.two, diff.three, diff.four))
  }

fourNormalDiffs <-
  function(X,
           as = as,
           al = al,
           bt = bt,
           L = L,
           h2 = h2,
           C,
           var.ratio = var.ratio,
           equalize.observed.vars = equalize.observed.vars) {
    ## inputs
    my.deltal <- 1 / (1 + exp(X[1]))
    my.tstar <- X[2]^2
    my.gs <- 1 / (1 + exp(X[3]))
    my.Ll <- L * (1 - my.gs)
    my.Ls <- L * my.gs
    
    ##Ls <- L-Ll
    my.bs <- 1 / (1 + exp(X[4]))
    my.Ls <- L - my.Ll
    ##my.L  <- my.Ll + my.Ls
    
    ## unscaled small effect stuff
    my.ys <- log((1 + my.bs) / (1 - my.bs))
    my.raw.Vas <- 8 * Ne * u * my.Ls * as ^ 2 * my.bs / my.ys
    
    
    ## raw large effect stuff
    my.s <- my.deltal * C
    my.yl <- 4 * Ne * my.s #originally 2*Ne*my.s
    bal <- ifelse(my.yl > 20, 1, (exp(my.yl) - 1) / (exp(my.yl) + 1))
    my.raw.Val <- 8 * Ne * u * my.Ll * al ^ 2 * bal / my.yl
    
    ## get total variance in original units
    my.raw.Vt <- (my.raw.Vas + my.raw.Val) / h2
    
    ## standardized effects
    my.std.as <- as / sqrt (my.raw.Vt)
    my.std.al <- al / sqrt (my.raw.Vt)
    
    ## standardize variances
    my.std.Vas <- my.raw.Vas / my.raw.Vt
    my.std.Val <- my.raw.Val / my.raw.Vt
    ##  my.std.norm.sd <- sqrt(1 - h2 + my.std.Vas)
    
    ## standardized thr dens
    my.ft <- my.ys / (4 * Ne * C * my.std.as)
    
    ##my.range <- seq(qpois(1e-8,my.mean.nl),qpois(1-1e-8,my.mean.nl))
    my.prot <- 1 - pnorm(my.tstar, 0, 1)
    my.risk <- 1 - pnorm(my.tstar - my.std.al, 0, 1)
    
    ## global stuff
    my.meana <- (my.std.as * my.Ls + my.std.al * my.Ll) / L
    my.maxg <- 2 * my.meana * L
    my.prev <- 1 - pnorm(my.tstar, 0, 1)
    my.risk.var <- (my.prev * (1 - my.prev))
    
    ## differences
    my.deltal.tild <- my.risk - my.prot
    diff.one <- my.deltal - my.deltal.tild
    
    my.ft.tild <- dnorm(my.tstar, 0, 1)
    diff.two <- my.ft - my.ft.tild
    
    if (equalize.observed.vars) {
      my.Vol <- 2 * my.deltal ^ 2 * my.Ll * u * bal / my.s
      my.Vos <- my.std.Vas * my.ft ^ 2
      diff.three <- (my.Vol - var.ratio * my.Vos) / my.risk.var
    } else {
      diff.three <- my.std.Val - var.ratio * my.std.Vas
    }
    
    bt.tild <-
      (2 * my.Ls * my.bs * my.std.as + 2 * my.Ll * 1 * my.std.al) / (my.maxg)
    diff.four <- bt - bt.tild
    
    ##diff.five <- my.L*my.meana - Lmeana
    
    return(c(diff.one, diff.two, diff.three, diff.four))
  }

solveTwoEffect <-
  function(bt,
           bs = bt,
           Ne,
           as,
           al,
           L,
           gs,
           last.tstar = NULL,
           h2 = NULL,
           Ve = NULL,
           u,
           C,
           LL.soln = FALSE,
           var.ratio = NULL,
           equalize.observed.vars = FALSE) {
    ## bs:  asymmetry for small effects
    ## Ne:  pop size
    ## Ls:  number of small effect loci
    ## Ll:  number of large effect loci
    ## as:  small effect size
    ## al:  large effect size
    ## r2n: ratio of Ve to Vas (small effect variance)
    ## u:   mutation rate
    ## T:   threshold position
    ## C:   cost
    ## h2l: ratio of large effect variance to total genetic variance (not currently used)
    
    ## transformed params
    Ls <- L * gs
    Ll <- L * (1 - gs)
    
    ## single effect solution
    single.norm.y <- log((1 + bt / (1 - bt)))
    single.norm.ft <- 1 / (4 * Ne * C) * single.norm.y
    single.norm.Vg <- 8 * Ne * L * u * bt / log((1 + bt / (1 - bt)))
    if (is.null(Ve)) {
      norm.Vt <- single.norm.Vg / h2
    } else {
      norm.Vt <- single.norm.Vg + Ve
    }
    std.dev.tot <- sqrt (norm.Vt)
    single.norm.as.std <- 1 / std.dev.tot
    single.norm.dens <-
      2 * L * u * bt * single.norm.as.std / (h2 * C)
    single.norm.std.thr <- dnorminv(single.norm.dens)
    single.norm.prev <- 1 - pnorm(single.norm.std.thr)
    
    
    init.al.std <- single.norm.as.std * al
    init.deltas <- single.norm.as.std * single.norm.dens
    init.deltal <-
      pnorm(single.norm.std.thr) - pnorm(single.norm.std.thr - init.al.std)
    init.sl <- init.deltal * C
    init.yl <- 4 * Ne * init.sl
    init.mean.nl <- 2 * Ll * u / init.sl
    init.Vas <-
      8 * Ne * u * Ls * single.norm.as.std^2 * bs / single.norm.y
    init.Val <- init.al.std ^ 2 * init.mean.nl
    init.Vt <- (init.Vas + init.Val) / h2
    
    sec.as.std <- single.norm.as.std / sqrt(init.Vt)
    sec.al.std <- init.al.std / sqrt(init.Vt)
    sec.Vas <- 8 * Ne * u * Ls * sec.as.std ^ 2 * bs / single.norm.y
    sec.Val <- sec.al.std ^ 2 * init.mean.nl
    sec.norm.sd <- sqrt (1 - sec.Val)
    sec.deltal <-
      pnorm(single.norm.std.thr) - pnorm(single.norm.std.thr - sec.al.std)
    sec.sl <- sec.deltal * C
    sec.yl <- 4 * Ne * sec.sl
    sec.mean.nl <- 2 * Ll * u / sec.sl
    sec.norm.dens <-
      log((1 + bs) / (1 - bs)) / (4 * Ne * C) * sqrt(init.Vt)
    sec.tstar <- dnorminv(sec.norm.dens)
    sec.prev <- 1 - pnorm(sec.tstar)
    
    my.seq <- seq(0, 12*sec.tstar, length.out = 1000)
    pois.dens.seq <- sapply(my.seq, function(X)
      dPoisConv(X, sec.mean.nl, sec.norm.sd, sec.al.std))
    single.norm.dens - dPoisConv(12*sec.tstar, sec.mean.nl, sec.norm.sd, sec.al.std)
    this.diff <- single.norm.dens - pois.dens.seq
    this.tf <- this.diff < 0
    this.min <- my.seq[min(which(this.tf))]
    
    init.soln <- uniroot(
      f = function(THR) {
        ftild <- dPoisConv(
          t = THR,
          lambda = sec.mean.nl,
          norm.sd = sec.norm.sd,
          alphal = sec.al.std
        )
        single.norm.dens - ftild
      },
      interval = c(this.min, 100 * sec.norm.sd)
    )
    init.tstar <- init.soln$root
    
    
    ## get 2d solution
    soln <- nleqslv(
      x = c(log((1 - sec.deltal) / sec.deltal), init.tstar),
      #x=c(0,init.tstar),
      fn = function(X) {
        my.deltal <- 1 / (1 + exp(X[1]))
        my.tstar <- X[2]
        
        my.s <- my.deltal * C
        ## my.y <- 4*Ne*my.s ## originally 2*Ne*my.s
        my.mean.nl <- 2 * Ll * u / my.s
        my.Val <- sec.al.std ^ 2 * my.mean.nl
        my.norm.sd <- sqrt (1 - my.Val)
        my.range <-
          seq(qpois(1e-8, my.mean.nl), qpois(1 - 1e-8, my.mean.nl))
        my.prot <-
          pPoisConv(my.tstar,
                    my.mean.nl,
                    my.norm.sd,
                    alphal = sec.al.std,
                    risk.allele = FALSE)
        my.risk <-
          pPoisConv(my.tstar,
                    my.mean.nl,
                    my.norm.sd,
                    alphal = sec.al.std,
                    risk.allele = TRUE)
        my.deltal.tild <- my.risk - my.prot
        ft.tild <-
          dPoisConv(my.tstar,
                    my.mean.nl,
                    my.norm.sd,
                    alphal = sec.al.std,
                    risk.allele = FALSE)
        diff.one <- my.deltal - my.deltal.tild
        diff.two <- single.norm.dens - ft.tild
        c(diff.one, diff.two)
      },
      control = list(scalex = c(1, 1), maxit = 400)
    )
    trans.deltal <- soln$x[1]
    deltal <- 1 / (1 + exp(trans.deltal))
    new.init.tstar <- soln$x[2]
    
    
    
    if (LL.soln) {
      init.trans.deltal <- trans.deltal
##      init.Ll <- Ll
##      init.Ls <- Ls
      init.trans.gs <- log((1 - gs) / gs)
      init.trans.bs <- log((1 - bs) / bs)
      
      
      ## get 4d Poisson convolution solution
      soln <- nleqslv(
        x = c(
          init.trans.deltal,
          new.init.tstar,
          init.trans.gs,
          init.trans.bs
        ),
        fn = function(X)
          fourPoissonDiffs(X, as, al, bt, L, h2, C, var.ratio, equalize.observed.vars),
        control = list(
          scalex = c(1, 1, 1, 1),
          allowSingular = FALSE
        )
      )
      if (soln$termcd == 4)
        return(NA)
      trans.deltal <- soln$x[1]
      deltal <- 1 / (1 + exp(trans.deltal))
      tstar <- soln$x[2]
      trans.gs <- soln$x[3]
      gs <- 1 / (1 + exp(trans.gs))
      Ls <- L * gs
      Ll <- L * (1 - gs)
      trans.bs <- soln$x[4]
      bs <- 1 / (1 + exp(trans.bs))
      ys <- log((1 + bs) / (1 - bs))
      
      if (TRUE) {
        ## get 4d Normal solution
        soln <- nleqslv(
          x = c(trans.deltal, sqrt(tstar), trans.gs, trans.bs),
          fn = function(X)
            fourNormalDiffs(X, as, al, bt, L, h2, C, var.ratio, equalize.observed.vars),
          control = list(
            scalex = c(1, 1 , 1, 1),
            allowSingular = FALSE
          )
        )
        multi.norm.trans.deltal <- soln$x[1]
        multi.norm.deltal <- 1 / (1 + exp(trans.deltal))
        multi.norm.tstar <- soln$x[2]^2
        multi.norm.trans.gs <- soln$x[3]
        multi.norm.gs <- 1 / (1 + exp(multi.norm.trans.gs))
        multi.norm.Ls <- L * multi.norm.gs
        multi.norm.Ll <- L * (1 - multi.norm.gs)
        multi.norm.bs <- 1 / (1 + exp(soln$x[4]))
        multi.norm.prev <- 1 - pnorm(multi.norm.tstar, 0, 1)
        
        
        multi.norm.ys <-
          log((1 + multi.norm.bs) / (1 - multi.norm.bs))
        multi.norm.Vas <-
          8 * Ne * u * multi.norm.Ls * as ^ 2 * multi.norm.bs / multi.norm.ys
        multi.norm.mean.nl <-
          2 * multi.norm.Ll * u / (multi.norm.deltal * C)
        multi.norm.Val <- al ^ 2 * multi.norm.mean.nl
        multi.norm.Vt <- (multi.norm.Vas + multi.norm.Val) / h2
        
        
      }
    }
    
    ##    L <- Ls+Ll
    ##    meana <- (as*Ls + al*Ll)/L
    raw.Vas <-
      8 * Ne * u * Ls * as ^ 2 * bs / ys # originall 4*Ne*u*Ls*as^2*bs/ys
    mean.nl <- 2 * Ll * u / (deltal * C)
    raw.Val <- al ^ 2 * mean.nl
    raw.Vt <- (raw.Vas + raw.Val) / h2
    std.as <- as / sqrt (raw.Vt)
    std.al <- al / sqrt (raw.Vt)
    std.Vas <- raw.Vas / raw.Vt
    std.Val <- raw.Val / raw.Vt
    
    ## Val <- 2*al^2*mean.nl
    ft <- ys / (4 * Ne * C)
    std.ft <- ys / (4 * Ne * C * std.as)
    deltas <- std.as * std.ft
    std.Vg <- std.Vas + std.Val
    norm.Va <- 1 - h2 + std.Vas
    norm.sd <- sqrt(norm.Va)
    
    ## risk scale variances
    Vos <- std.Vas * std.ft ^ 2
    Vol <- deltal ^ 2 * mean.nl
    Vrg <- Vos + Vol
    
    ## bulk stuff
    prev <-
      as.numeric(pPoisConv(tstar, mean.nl, norm.sd, alphal = std.al))
    my.range <- seq(qpois(1e-8, mean.nl), qpois(1 - 1e-8, mean.nl))
    ldens <-
      dpois(my.range, mean.nl) ## density on large effect liability
    pdil <-
      pnorm(tstar, std.al * my.range, norm.sd, lower.tail = FALSE) ## prev among inds with i large effect alleles
    pidl <-
      ldens * pdil / prev ## density on number of large effect alleles conditional on having disease
    
    
    ## does mut-sel balance actually hold??
    mutl <- 2 * Ll * std.al * u
    muts <- 2 * Ls * std.as * u * bs
    mutall <- muts + mutl
    
    sell <- mean.nl * std.al * deltal * C
    sels <- std.Vas * std.ft * C
    selall <- sell + sels
    stopifnot(abs((mutall - selall) / (mutall + selall)) < 1e-8)
    
    
    tol <- 1e-4
    min.gl <- uniroot(function(X)
      tol - (1 - pPoisConv(X, mean.nl, norm.sd, alphal = std.al)),
      interval = c(-10 * tstar, 10 * tstar))$root
    max.gl <- uniroot(
      function(X)
        tol - pPoisConv(X, mean.nl, norm.sd, alphal = std.al),
      interval = c(0, 12 * tstar)
    )$root
    seq.li <- seq(min.gl, max.gl, length.out = 1000)
    li.dense <-
      sapply(seq.li, function(G)
        dPoisConv(G, mean.nl, norm.sd, alphal = std.al))
    
    
    ## heritabilies on liability scale
    h2s <- std.Vas / (std.Val + std.Vas + (1 - h2))
    h2l <- std.Val / (std.Val + std.Vas + (1 - h2))
    
    
    ## liability scale h2 estimation via normal approx
    phit <- dnorm(qnorm(1 - prev))
    risk.var <- (prev * (1 - prev))
    
    ## observed scale
    h2os <- Vos / risk.var
    h2ol <- Vol / risk.var
    h2o <- Vrg / risk.var
    prs <- Vos / (Vos + Vol)
    
    ## naive estimates on liability scale
    h2s.est <- h2os * risk.var / phit ^ 2
    h2l.est <- h2ol * risk.var / phit ^ 2
    h2all.est <- h2o * risk.var / phit ^ 2
    
    
    
    implied.al.norm <- deltal / phit
    implied.al.true <- deltal / std.ft
    
    
    ## other stuff
    #L <- Ll+Ls
    
    meana <- (std.as * Ls + std.al * Ll) / L
    maxg <- 2 * meana * L
    bt <- (2 * Ls * bs * std.as + 2 * Ll * 1 * std.al) / maxg
    
    if (multi.norm.prev>0.5) recover()
    
    return(
      list(
        soln = soln,
        ft = ft,
        std.ft = std.ft,
        phit = phit,
        bs = bs,
        bt = bt,
        as = as,
        al = al,
        std.as = std.as,
        std.al = std.al,
        meana = meana,
        maxg = maxg,
        bt = bt,
        my.range = my.range,
        deltas = deltas,
        deltal = deltal,
        ss = deltas * C,
        sl = deltal * C,
        ys = 4 * Ne * deltas * C,
        #originally 2*Ne*deltas*C
        yl = 4 * Ne * deltal * C,
        #originally 2*Ne*deltal*C
        Ne = Ne,
        u = u,
        C = C,
        ors = (prev + deltas) / prev * (1 - prev) / (1 - prev - deltas),
        orl = (prev + deltal) / prev * (1 - prev) / (1 - prev - deltal),
        tstar = tstar,
        multi.norm.tstar = multi.norm.tstar,
        mean.nl = mean.nl,
        prev = prev,
        single.norm.prev = single.norm.prev,
        multi.norm.prev = multi.norm.prev,
        raw.Ve = (1 - h2) * (as / std.as) ^ 2,
        std.Vas = std.Vas,
        std.Val = std.Val,
        std.Vg = std.Vg,
        norm.sd = norm.sd,
        Vos = Vos,
        Vol = Vol,
        h2s = h2s,
        h2l = h2l,
        h2 = h2,
        h2os = h2os,
        h2ol = h2ol,
        h2o = h2o,
        h2s.est = h2s.est,
        ## small effect h2 you would infer by naive application of ltm methods
        h2l.est = h2l.est,
        ## large effect h2 you would infer by naive application of ltm methods
        h2all.est = h2all.est,
        ## total h2 you would infer by naive application of ltm methods
        pgs.est = h2s.est / h2all.est,
        implied.al.norm = implied.al.norm,
        implied.al.true = implied.al.true,
        pgs = std.Vas / (std.Val + std.Vas),
        pgl = std.Val / (std.Val + std.Vas),
        Vos = Vos,
        Vol = Vol,
        pos = Vos / (Vos + Vol),
        pol = Vol / (Vos + Vol),
        ldens = ldens,
        pdil = pdil,
        pidl = pidl,
        gs = gs,
        Ls = Ls,
        Ll = Ll,
        li.dense = li.dense,
        var.ratio = var.ratio
      )
    )
  }




makeOutput = function(soln) {
  output = numeric()
  output['as'] <- soln$as  ## small effect size
  output['al'] <- soln$al  ## small effect size
  output['std.as'] <-
    soln$std.as  ## standardized small effect size
  output['std.al'] <-
    soln$std.al  ## standardized large effect size
  output['phit'] <-
    soln$phit  ## normal equivalent standardized density
  output['ft'] <-
    soln$ft  ## raw density
  output['std.ft'] <-
    soln$std.ft  ## std density
  output['ys'] <-
    soln$ys  ## small scaled selection coefficient
  output['yl'] <-
    soln$yl  ## large scaled selection coefficient
  output['Ls'] <- soln$Ls  ## numer of small effect loci
  output['gs'] <- soln$gs  ## fraction of small effect loci
  output['bs'] <- soln$bs  ## b for small effect loci
  output['rhos'] <-
    1 / 2 - output['bs'] / 2  ## rho for small effect loci
  output['al'] <- soln$al  ## large effect size
  output['Ll'] <-
    round(soln$Ll, 0)  ## numer of large effect loci
  output['bl'] <- 1              ## b for large effect loci
  output['rhol'] <- 0            ## rho for large effect loci
  output['bt'] <- soln$bt  ## total b
  output['rhot'] <-  1 / 2 - soln$bt / 2  ## total rho
  
  output['Ne'] <- soln$Ne  ## population size
  output['u'] <- soln$u    ## mutation rate
  output['C'] <- soln$C    ## cost of disease
  output['Ve'] <- soln$raw.Ve  ## environmental variance
  output['norm.sd'] <- soln$norm.sd  ## environmental variance
  output['h2s'] <-
    soln$h2s  ## small effect h2 on liability scale
  output['h2l'] <-
    soln$h2l  ## large effect h2 on liability scale
  output['h2'] <- output['h2s'] + output['h2l']
  output['h2os'] <-
    soln$h2os  ## small effect h2 on observed scale
  output['h2ol'] <-
    soln$h2ol  ## large effect h2 on observed scale
  output['pgl'] <-
    soln$pgl  ## large effect h2 on observed scale
  output['h2s.est'] <-
    soln$h2s.est  
  output['h2l.est'] <-
    soln$h2l.est  
  output['h2all.est'] <-
    soln$h2all.est  ## large effect h2 on observed scale
  output['prev'] <- soln$prev ## prevalence
  output['multi.norm.prev'] <- soln$multi.norm.prev ## prevalence
  output['deltas'] <- soln$deltas ## risk small effect
  output['deltal'] <- soln$deltal ## risk small effect
  output['maxG'] <-
    2 * output['Ls'] * output['as'] + 2 * output['Ll'] * output['al']
  output['bigU_small']  <-
    2 * output['Ls'] * u * output['as'] * output['bs']
  output['bigU_large']  <-
    2 * output['Ll'] * u * output['al']
  output['bigU'] <- output['maxG'] * u * output['bt']
  output['thr'] <-
    2 * output['rhos'] * output['Ls'] * output['as'] + 1
  output['mean.nl'] <-
    soln$mean.nl  ## average number of large effect alleles per individual
  output['tstar'] <-
    soln$tstar  ## threshold distance in standardized units
  output['multi.norm.tstar'] <-
    soln$multi.norm.tstar
  output['var.ratio'] <-
    soln$var.ratio
  return(output)
}

# 
# solveTwoEffectNew <-
#   function(bt,
#            bs = bt,
#            Ne,
#            as,
#            al,
#            L,
#            gs,
#            normalSoln = FALSE,
#            last.tstar = NULL,
#            h2 = NULL,
#            Ve = NULL,
#            u,
#            C,
#            LL.soln = FALSE,
#            var.ratio = NULL,
#            equalize.observed.vars = FALSE) {
#     ## bs:  asymmetry for small effects
#     ## Ne:  pop size
#     ## Ls:  number of small effect loci
#     ## Ll:  number of large effect loci
#     ## as:  small effect size
#     ## al:  large effect size
#     ## r2n: ratio of Ve to Vas (small effect variance)
#     ## u:   mutation rate
#     ## T:   threshold position
#     ## C:   cost
#     ## h2l: ratio of large effect variance to total genetic variance (not currently used)
#     
#     if (recover.flag)
#       recover()
#     
#     ## transformed params
#     Ls <- L * gs
#     Ll <- L * (1 - gs)
#     
#     ## single effect solution
#     single.norm.y <- log((1 + bt / (1 - bt)))
#     single.norm.ft <- 1 / (4 * Ne * C) * single.norm.y
#     single.norm.Vg <- 8 * Ne * L * u * bt / log((1 + bt / (1 - bt)))
#     if (is.null(Ve)) {
#       norm.Vt <- single.norm.Vg / h2
#     } else {
#       norm.Vt <- single.norm.Vg + Ve
#     }
#     std.dev.tot <- sqrt (norm.Vt)
#     single.norm.as.std <- 1 / std.dev.tot
#     single.norm.dens <-
#       2 * L * u * bt * single.norm.as.std / (h2 * C)
#     
#     
#     single.norm.std.thr <- dnorminv(single.norm.dens)
#     single.norm.prev <- 1 - pnorm(single.norm.std.thr)
#     init.as.std <- single.norm.as.std 
#     init.al.std <- single.norm.as.std * al
#     init.deltas <- single.norm.as.std * single.norm.dens
#     init.deltal <- pnorm(single.norm.std.thr) - pnorm(single.norm.std.thr - init.al.std)
#     init.sl <- init.deltal * C
#     init.yl <- 4 * Ne * init.sl 
#     init.mean.nl <- 2 * Ll * u / init.sl
#     init.Vas <-
#       8 * Ne * u * Ls * single.norm.as.std ^ 2 * bs / single.norm.y
#     init.Val <- init.al.std ^ 2 * init.mean.nl
#     init.Vt <- (init.Vas + init.Val) / h2
#     
#     ## rescale to achieve desired heritability
#     sec.as.std <- init.as.std
#     sec.al.std <- init.al.std
#     sec.ys <- single.norm.y
#     sec.dens <- single.norm.dens
#     sec.mean.nl <- init.mean.nl
#     sec.std.thr <- single.norm.std.thr
#     sec.Vt <- init.Vt
#     
#     for(i in 1:5) {
#       sec.as.std <- sec.as.std / sqrt(sec.Vt)
#       sec.al.std <- sec.al.std / sqrt(sec.Vt)
#       sec.ys <- 4 * Ne * sec.as.std * sec.dens * C
#       sec.Vas <- 8 * Ne * u * Ls * sec.as.std ^ 2 * bs / sec.ys
#       sec.Val <- sec.al.std ^ 2 * sec.mean.nl
#       sec.Vt <- (sec.Vas + sec.Val) / h2
#       sec.norm.sd <- sqrt (1 - sec.Val)
#       sec.deltal <-
#         pnorm(sec.std.thr) - pnorm(sec.std.thr - sec.al.std)
#       sec.sl <- sec.deltal * C
#       sec.yl <- 4 * Ne * sec.sl
#       sec.mean.nl <- 2 * Ll * u / sec.sl
#       sec.norm.dens <-
#         log((1 + bs) / (1 - bs)) / (4 * Ne * C) * sqrt(sec.Vt)
#       sec.tstar <- dnorminv(sec.norm.dens)
#       sec.prev <- 1 - pnorm(sec.tstar)
#     }
#     
#     my.seq <- seq(0, 12, length.out = 100)
#     this.diff <-
#       sapply(my.seq, function(X)
#         dPoisConv(X, sec.mean.nl, sec.norm.sd, sec.al.std)) > single.norm.dens
#     this.min <- my.seq[min(which(this.diff))]
#     
#     init.soln <- uniroot(
#       f = function(THR) {
#         ftild <- dPoisConv(
#           t = THR,
#           lambda = sec.mean.nl,
#           norm.sd = sec.norm.sd,
#           alphal = sec.al.std
#         )
#         single.norm.dens - ftild
#       },
#       interval = c(this.min, 12 * sec.norm.sd)
#     )
#     init.tstar <- init.soln$root
#     
#     
#     ## get 2d solution
#     soln <- nleqslv(
#       x = c(log((1 - sec.deltal) / sec.deltal), init.tstar),
#       #x=c(0,init.tstar),
#       fn = function(X) {
#         my.deltal <- 1 / (1 + exp(X[1]))
#         my.tstar <- X[2]
#         
#         my.s <- my.deltal * C
#         ## my.y <- 4*Ne*my.s ## originally 2*Ne*my.s
#         my.mean.nl <- 2 * Ll * u / my.s
#         my.Val <- sec.al.std ^ 2 * my.mean.nl
#         my.norm.sd <- sqrt (1 - my.Val)
#         my.range <-
#           seq(qpois(1e-8, my.mean.nl), qpois(1 - 1e-8, my.mean.nl))
#         my.prot <-
#           pPoisConv(my.tstar,
#                     my.mean.nl,
#                     my.norm.sd,
#                     alphal = sec.al.std,
#                     risk.allele = FALSE)
#         my.risk <-
#           pPoisConv(my.tstar,
#                     my.mean.nl,
#                     my.norm.sd,
#                     alphal = sec.al.std,
#                     risk.allele = TRUE)
#         my.deltal.tild <- my.risk - my.prot
#         ft.tild <-
#           dPoisConv(my.tstar,
#                     my.mean.nl,
#                     my.norm.sd,
#                     alphal = sec.al.std,
#                     risk.allele = FALSE)
#         diff.one <- my.deltal - my.deltal.tild
#         diff.two <- single.norm.dens - ft.tild
#         c(diff.one, diff.two)
#       },
#       control = list(scalex = c(1, 1), maxit = 400)
#     )
#     trans.deltal <- soln$x[1]
#     deltal <- 1 / (1 + exp(trans.deltal))
#     new.init.tstar <- soln$x[2]
#     
#     
#     
#     if (LL.soln) {
#       ## init.trans.deltal <- trans.deltal
#       ## init.Ll <- Ll
#       ## init.Ls <- Ls
#       
#       
#       
#       
#       if (normalSoln) {
#         
#         ## get 4d Normal solution
#         sec.trans.deltal <- log((1 - sec.deltal) / sec.deltal)
#         init.trans.gs <- log((1 - gs) / gs)
#         init.trans.bs <- log((1 - bs) / bs)
#         soln <- nleqslv(
#           x = c(
#             sec.trans.deltal,
#             sec.tstar,
#             init.trans.gs,
#             init.trans.bs
#           ),
#           fn = function(X)
#             fourNormalDiffs(X, as, al, bt, L, h2, C, var.ratio, equalize.observed.vars),
#           control = list(
#             scalex = c(1, 1, 1, 1),
#             allowSingular = FALSE
#           )
#         )
#         trans.deltal <- soln$x[1]
#         deltal <- 1 / (1 + exp(trans.deltal))
#         tstar <- soln$x[2]
#         trans.gs <- soln$x[3]
#         gs <- 1 / (1 + exp(trans.gs))
#         Ls <- L * gs
#         Ll <- L * (1 - gs)
#         bs <- 1 / (1 + exp(soln$x[4]))
#         prev <- 1 - pnorm(tstar, 0, 1)
#         ys <- log((1 + bs) / (1 - bs))
#       } else {
#         ## get 4d Poisson convolution solution
#         soln <- nleqslv(
#           x = c(
#             init.trans.deltal,
#             new.init.tstar,
#             init.trans.gs,
#             init.trans.bs
#           ),
#           fn = function(X)
#             fourPoissonDiffs(X, as, al, bt, L, h2, C, var.ratio, equalize.observed.vars),
#           control = list(
#             scalex = c(1, 1 / new.init.tstar, 1, 1),
#             allowSingular = FALSE
#           )
#         )
#         if (soln$termcd == 4)
#           return(NA)
#         trans.deltal <- soln$x[1]
#         deltal <- 1 / (1 + exp(trans.deltal))
#         tstar <- soln$x[2]
#         trans.gs <- soln$x[3]
#         gs <- 1 / (1 + exp(trans.gs))
#         Ls <- L * gs
#         Ll <- L * (1 - gs)
#         trans.bs <- soln$x[4]
#         bs <- 1 / (1 + exp(trans.bs))
#         ys <- log((1 + bs) / (1 - bs))
#         
#         
#         tol <- 1e-5
#         min.gl <- uniroot(function(X)
#           tol - (1 - pPoisConv(X, mean.nl, norm.sd, alphal = std.al)),
#           interval = c(-10 * tstar, 10 * tstar))$root
#         max.gl <- uniroot(
#           function(X)
#             tol - pPoisConv(X, mean.nl, norm.sd, alphal = std.al),
#           interval = c(0, 12 * tstar)
#         )$root
#         seq.li <- seq(min.gl, max.gl, length.out = 1000)
#         li.dense <-
#           sapply(seq.li, function(G)
#             dPoisConv(G, mean.nl, norm.sd, alphal = std.al))
#         ## bulk stuff
#         prev <-
#           as.numeric(pPoisConv(tstar, mean.nl, norm.sd, alphal = std.al))
#         my.range <-
#           seq(qpois(1e-8, mean.nl), qpois(1 - 1e-8, mean.nl))
#         ldens <-
#           dpois(my.range, mean.nl) ## density on large effect liability
#         pdil <-
#           pnorm(tstar, std.al * my.range, norm.sd, lower.tail = FALSE) ## prev among inds with i large effect alleles
#         pidl <-
#           ldens * pdil / prev ## density on number of large effect alleles conditional on having disease
#         
#         
#       }
#       
#     }
#     
#     ##    L <- Ls+Ll
#     ##    meana <- (as*Ls + al*Ll)/L
#     raw.Vas <-
#       8 * Ne * u * Ls * as ^ 2 * bs / ys # originall 4*Ne*u*Ls*as^2*bs/ys
#     mean.nl <- 2 * Ll * u / (deltal * C)
#     raw.Val <- al ^ 2 * mean.nl
#     raw.Vt <- (raw.Vas + raw.Val) / h2
#     std.as <- as / sqrt (raw.Vt)
#     std.al <- al / sqrt (raw.Vt)
#     std.Vas <- raw.Vas / raw.Vt
#     std.Val <- raw.Val / raw.Vt
#     
#     ## Val <- 2*al^2*mean.nl
#     ft <- ys / (4 * Ne * C)
#     std.ft <- ys / (4 * Ne * C * std.as)
#     deltas <- std.as * std.ft
#     std.Vg <- std.Vas + std.Val
#     norm.Va <- 1 - h2 + std.Vas
#     norm.sd <- sqrt(norm.Va)
#     ## as.std  <- as / sqrt( Vt )
#     ## al.std  <- al / sqrt( Vt )
#     ## ft.std <- ft * sqrt ( Vt )
#     
#     ## risk scale variances
#     Vos <- std.Vas * std.ft ^ 2
#     Vol <- deltal ^ 2 * mean.nl
#     Vrg <- Vos + Vol
#     
#     
#     ## does mut-sel balance actually hold??
#     mutl <- 2 * Ll * std.al * u
#     muts <- 2 * Ls * std.as * u * bs
#     mutall <- muts + mutl
#     
#     sell <- mean.nl * std.al * deltal * C
#     sels <- std.Vas * std.ft * C
#     selall <- sell + sels
#     stopifnot(abs((mutall - selall) / (mutall + selall)) < 1e-8)
#     
#     
#     ##recover()
#     ## heritabilies on liability scale
#     h2s <- std.Vas / (std.Val + std.Vas + (1 - h2))
#     h2l <- std.Val / (std.Val + std.Vas + (1 - h2))
#     
#     
#     ## liability scale h2 estimation via normal approx
#     phit <- dnorm(qnorm(1 - prev))
#     risk.var <- (prev * (1 - prev))
#     
#     ## observed scale
#     h2os <- Vos / risk.var
#     h2ol <- Vol / risk.var
#     h2o <- Vrg / risk.var
#     prs <- Vos / (Vos + Vol)
#     
#     ## naive estimates on liability scale
#     h2s.est <- h2os * risk.var / phit ^ 2
#     h2l.est <- h2ol * risk.var / phit ^ 2
#     h2all.est <- h2o * risk.var / phit ^ 2
#     
#     
#     
#     implied.al.norm <- deltal / phit
#     implied.al.true <- deltal / std.ft
#     
#     
#     ## other stuff
#     #L <- Ll+Ls
#     
#     meana <- (std.as * Ls + std.al * Ll) / L
#     maxg <- 2 * meana * L
#     bt <- (2 * Ls * bs * std.as + 2 * Ll * 1 * std.al) / maxg
#     
#     return(
#       list(
#         soln = soln,
#         ft = ft,
#         std.ft = std.ft,
#         phit = phit,
#         bs = bs,
#         bt = bt,
#         as = as,
#         al = al,
#         std.as = std.as,
#         std.al = std.al,
#         meana = meana,
#         maxg = maxg,
#         bt = bt,
#         ##my.range = my.range,
#         deltas = deltas,
#         deltal = deltal,
#         ss = deltas * C,
#         sl = deltal * C,
#         ys = 4 * Ne * deltas * C,
#         yl = 4 * Ne * deltal * C,
#         Ne = Ne,
#         u = u,
#         C = C,
#         ors = (prev + deltas) / prev * (1 - prev) / (1 - prev - deltas),
#         orl = (prev + deltal) / prev * (1 - prev) / (1 - prev - deltal),
#         tstar = tstar,
#         mean.nl = mean.nl,
#         prev = prev,
#         single.norm.prev = single.norm.prev,
#         raw.Ve = (1 - h2) * (as / std.as) ^ 2,
#         std.Vas = std.Vas,
#         std.Val = std.Val,
#         std.Vg = std.Vg,
#         norm.sd = norm.sd,
#         Vos = Vos,
#         Vol = Vol,
#         h2s = h2s,
#         h2l = h2l,
#         h2 = h2,
#         h2os = h2os,
#         h2ol = h2ol,
#         h2o = h2o,
#         h2s.est = h2s.est,
#         ## small effect h2 you would infer by naive application of ltm methods
#         h2l.est = h2l.est,
#         ## large effect h2 you would infer by naive application of ltm methods
#         h2all.est = h2all.est,
#         ## total h2 you would infer by naive application of ltm methods
#         pgs.est = h2s.est / h2all.est,
#         implied.al.norm = implied.al.norm,
#         implied.al.true = implied.al.true,
#         pgs = std.Vas / (std.Val + std.Vas),
#         pgl = std.Val / (std.Val + std.Vas),
#         Vos = Vos,
#         Vol = Vol,
#         pos = Vos / (Vos + Vol),
#         pol = Vol / (Vos + Vol),
#         # ldens = ldens,
#         # pdil = pdil,
#         # pidl = pidl,
#         gs = gs,
#         Ls = Ls,
#         Ll = Ll
#         ##li.dense = li.dense
#       )
#     )
#   }
# 
# 
# 




## solveTwoEffectFixed_bs <- function(bs,bt,Ne,Ls,Ll,as,al,Lmeana,L.init,r2n=NULL,Ve=NULL,u,C,LL.soln=FALSE,var.ratio=NULL,equalize.observed.vars=FALSE){

##     ## bs:  asymmetry for small effects
##     ## Ne:  pop size
##     ## Ls:  number of small effect loci
##     ## Ll:  number of large effect loci
##     ## as:  small effect size
##     ## al:  large effect size
##     ## r2n: ratio of Ve to Vas (small effect variance)
##     ## u:   mutation rate
##     ## T:   threshold position
##     ## C:   cost
##     ## h2l: ratio of large effect variance to total genetic variance (not currently used)

##     ## bt <- bs

##     if(recover.flag) recover()


##     if(LL.soln)
##         stopifnot(!is.null(var.ratio))

##     if((!is.null(r2n))&(!is.null(Ve))){
##         stop('either r2n or Ve must be NULL')
##     }

##     init.Ls <- Ls
##     ## solve for ft
##     ys <- log((1+bs)/(1-bs))
##     ft <- ys*(4*Ne*C)^-1  #originally ys*(2*Ne*C)^-1

##     ## get small effect additive genetic variance
##     Vas <- 8*Ne*u*L.init*as^2*bs/ys  ## originally 4*Ne*u*Ls*as^2*bs/ys
##     if(is.null(Ve))
##         Ve <- Vas*(1-r2n)/r2n
##     norm.sd <- sqrt(Vas+Ve)

##     ## initial guess for deltal
##     init.deltal <- 1##al*ft



##     ## compute selection coefficient
##     init.sl <- init.deltal*C
##     init.yl <- 4*Ne*init.sl  #originally 2*Ne*init.sl

##     if(init.yl<3)
##         warning('initial scaled selection coefficient is too small')


##     ## mean number of large effect alleles per individual
##     init.mean.nl <- 2*Ll*u/init.sl

##     ## solve for tstart conditional on initial guess of deltal

##     init.soln <- uniroot(
##         f=function(THR){
##             ftild <- dPoisConv(
##                 t=THR,
##                 lambda=init.mean.nl,
##                 norm.sd=norm.sd,
##                 alphal=al
##             )
##             ft-ftild
##         },
##         interval=c(0,12*norm.sd)
##     )
##     init.tstar <- init.soln$root
##     init.Ll <- Ll
##     init.L <- init.Ls + init.Ll
##     init.meana <- (as*init.Ls + al*init.Ll)/init.L
##     init.maxg <- 2*init.L*init.meana
##     init.bt <- (2*init.Ls*bs*as + 2*init.Ll*1*al)/(init.maxg)


##     ## get 2d solution
##     soln <- nleqslv(
##         x=c(0,init.tstar),
##         fn=function(X){
##             deltal <- 1/(1+exp(X[1]))
##             tstar <- X[2]
##             my.s <- deltal*C
##             my.y <- 4*Ne*my.s ## originally 2*Ne*my.s
##             mean.nl <- 2*Ll*u/my.s
##             my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##             prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##             risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)
##             deltal.tild <- risk - prot
##             ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##             Val <- 2*al*mean.nl
##             diff.one <- deltal-deltal.tild
##             diff.two <- ft-ft.tild
##             ##diff.three <- Val-Vas
##             ##diff.three <- 4*al^2*mean.nl - Vas*h2l/(1-h2l)
##             c(diff.one,diff.two)
##         },
##         control=list(scalex=c(1,1/init.tstar),maxit=400)
##     )
##     trans.deltal <- soln$x[1]
##     deltal <- 1/(1+exp(trans.deltal))
##     tstar <- soln$x[2]


##     if(LL.soln){

##         ## added the normal jitters because giving the solution as the initial condition led
##         ## it to somehow find other pathological solutions for some reason. Confused..
##         init.trans.deltal <- trans.deltal + rnorm(1,0,abs(trans.deltal)/1000)
##         new.init.tstar <- soln$x[2] + rnorm(1,0,tstar/1000)
##         init.Ll <- Ll + rnorm(1,0,Ll/1000)
##         init.Ls <- Ls
##         init.trans.bt <- log((1-init.bt)/init.bt)

##         j <- 1
##         ## get 4d solution
##         soln <- nleqslv(
##             x=c(init.trans.deltal,new.init.tstar,init.Ll,init.Ls),
##             fn=function(X){
##                 ## recover()
##                 j <- j+1
##                 ## inputs
##                 deltal <- 1/(1+exp(X[1]))
##                 tstar <- X[2]
##                 Ll <- X[3]
##                 ##Ls <- L-Ll
##                 ## bt <- 1/(1+exp(X[4]))
##                 Ls <- X[4]
##                 L  <- Ll + Ls

##                 ## small effect stuff
##                 ys <- log((1+bs)/(1-bs))
##                 Vas <- 8*Ne*u*Ls*as^2*bs/ys # originall 4*Ne*u*Ls*as^2*bs/ys
##                 if(is.null(Ve))
##                     Ve <- Vas*(1-r2n)/r2n
##                 norm.sd <- sqrt(Vas + Ve)

##                 ## large effect stuff
##                 my.s <- deltal*C
##                 my.y <- 4*Ne*my.s #originally 2*Ne*my.s
##                 mean.nl <- 2*Ll*u/my.s
##                 my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##                 prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##                 risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)

##                 ## global stuff
##                 meana <- (as*Ls + al*Ll)/L
##                 maxg <- 2*meana*L
##                 prev <- as.numeric(pPoisConv(tstar,mean.nl,sqrt(Vas+Ve),alphal=al))
##                 risk.var <- (prev*(1-prev))

##                 ## differences
##                 deltal.tild <- risk - prot
##                 diff.one <- deltal-deltal.tild

##                 ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##                 diff.two <- ft-ft.tild

##                 if(equalize.observed.vars){
##                     Vol <- 2*deltal^2*mean.nl
##                     Vos <- Vas*ft^2
##                     diff.three <- (Vol-var.ratio*Vos)/risk.var
##                 } else {
##                     Val <- 2*al^2*mean.nl
##                     diff.three <- Val-var.ratio*Vas
##                 }

##                 ##bt.tild <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)
##                 ## diff.four <- bt - bt.tild

##                 diff.four <- L*meana - Lmeana

##                 return(c(diff.one,diff.two,diff.three,diff.four))
##             },
##             control=list(scalex=c(1,1/new.init.tstar,1/init.Ll,1/init.maxg),maxit=800,allowSingular=FALSE)
##         )
##         if(recover.flag) recover()
##         trans.deltal <- soln$x[1]
##         deltal <- 1/(1+exp(trans.deltal))
##         tstar <- soln$x[2]
##         Ll <- soln$x[3]
##         ##bs <- 1/(1+exp(soln$x[4]))
##         Ls <- soln$x[4]

##         ys <- log((1+bs)/(1-bs))
##         ft <- ys*(4*Ne*C)^-1 #originally ys*(2*Ne*C)^-1

##         ## get small effect additive genetic variance
##         Vas <- 8*Ne*u*Ls*as^2*bs/ys  #originally 4*Ne*u*Ls*as^2*bs/ys
##         if(is.null(Ve))
##             Ve <- Vas*(1-r2n)/r2n
##         norm.sd <- sqrt(Vas+Ve)

##     }

##     L <- Ls+Ll
##     deltas <- as*ft
##     mean.nl <- 2*Ll*u/(deltal*C)
##     Val <- 2*al^2*mean.nl
##     Vt <- Vas + Ve + Val


##     ## risk scale variances
##     Vos <- Vas*ft^2
##     Vol <- 2*deltal^2*mean.nl
##     Vrg <- Vos + Vol

##     ## bulk stuff
##     prev <- as.numeric(pPoisConv(tstar,mean.nl,sqrt(Vas+Ve),alphal=al))
##     my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##     ldens <- dpois(my.range,mean.nl) ## density on large effect liability
##     pdil <- pnorm(tstar, al*my.range,norm.sd,lower.tail=FALSE) ## prev among inds with i large effect alleles
##     pidl <- ldens*pdil/prev ## density on number of large effect alleles conditional on having disease


##     ## does mut-sel balance actually hold??
##     mutl <- 2*Ll*al*u
##     muts <- 2*Ls*as*u*bs
##     mutall <- muts+mutl

##     sell <- mean.nl*al*deltal*C
##     sels <- Vas*ft*C
##     selall <- sell+sels
##     stopifnot(abs((mutall-selall)/(mutall+selall))<1e-8)


##     tol <- 1e-5
##     min.gl <- uniroot(
##         function(X)
##             tol-(1-pPoisConv(X,mean.nl,norm.sd,alphal=al)),
##         interval=c(-10*tstar,10*tstar)
##     )$root
##     max.gl <- uniroot(
##         function(X)
##             tol-pPoisConv(X,mean.nl,norm.sd,alphal=al),
##         interval=c(-10*tstar,10*tstar)
##     )$root
##     seq.li <- seq(min.gl,max.gl,length.out=1000)
##     li.dense <- sapply(seq.li,function(G) dPoisConv(G,mean.nl,norm.sd,alphal=al))

##     ## heritabilies on liability scale
##     h2s <- Vas/(Val+Vas/r2n)
##     h2l <- Val/(Val+Vas/r2n)
##     h2 <- (Vas+Val)/(Val+Vas/r2n)



##     ## liability scale h2 estimation via normal approx
##     phit <- dnorm(qnorm(1-prev))
##     risk.var <- (prev*(1-prev))

##     ## observed scale
##     h2os <- Vos/risk.var
##     h2ol <- Vol/risk.var
##     h2o <- Vrg/risk.var
##     prs <- Vos/(Vos+Vol)

##     ## naive estimates on liability scale
##     h2s.est <- h2os*risk.var/Vt*phit^-2
##     h2l.est <- h2ol*risk.var/Vt*phit^-2
##     h2all.est <- h2o*risk.var/Vt*phit^-2

##     implied.al.norm <- deltal/phit
##     implied.al.true <- deltal/ft


##     ## other stuff
##     #L <- Ll+Ls

##     meana <- (as*Ls + al*Ll)/L
##     maxg <- 2*meana*L
##     bt <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)

##     return(
##         list(
##             soln=soln,
##             ft=ft,
##             phit=phit,
##             bs=bs,
##             bt=bt,
##             as=as,
##             al=al,
##             meana=meana,
##             maxg=maxg,
##             bt=bt,
##             my.range=my.range,
##             deltas=deltas,
##             deltal=deltal,
##             ss=deltas*C,
##             sl=deltal*C,
##             ys=4*Ne*deltas*C,  #originally 2*Ne*deltas*C
##             yl=4*Ne*deltal*C,  #originally 2*Ne*deltal*C
##             Ne=Ne,
##             u=u,
##             C=C,
##             ors=(prev+deltas)/prev*(1-prev)/(1-prev-deltas),
##             orl=(prev+deltal)/prev*(1-prev)/(1-prev-deltal),
##             tstar=tstar,
##             mean.nl=mean.nl,
##             prev=prev,
##             Vas=Vas,
##             Ve=Ve,
##             Val=Val,
##             Vos=Vos,
##             Vol=Vol,
##             h2s=h2s,
##             h2l=h2l,
##             h2=h2,
##             h2os=h2os,
##             h2ol=h2ol,
##             h2o=h2o,
##             h2s.est=h2s.est, ## small effect h2 you would infer by naive application of ltm methods
##             h2l.est=h2l.est, ## large effect h2 you would infer by naive application of ltm methods
##             h2all.est=h2all.est, ## total h2 you would infer by naive application of ltm methods
##             pgs.est=h2s.est/h2all.est,
##             implied.al.norm=implied.al.norm,
##             implied.al.true=implied.al.true,
##             pgs=Vas/(Val+Vas),
##             pgl=Val/(Val+Vas),
##             Vos=Vos,
##             Vol=Vol,
##             pos=Vos/(Vos+Vol),
##             pol=Vol/(Vos+Vol),
##             ldens=ldens,
##             pdil=pdil,
##             pidl=pidl,
##             Ls=Ls,
##             Ll=Ll,
##             li.dense=li.dense
##         )
##     )
## }



## solveTwoEffectOld <- function(bs,bt,Ne,Ls,Ll,as,al,Lmeana,L.init,last.tstar=NULL,r2n=NULL,Ve=NULL,u,C,LL.soln=FALSE,var.ratio=NULL,equalize.observed.vars=FALSE){

##     ## bs:  asymmetry for small effects
##     ## Ne:  pop size
##     ## Ls:  number of small effect loci
##     ## Ll:  number of large effect loci
##     ## as:  small effect size
##     ## al:  large effect size
##     ## r2n: ratio of Ve to Vas (small effect variance)
##     ## u:   mutation rate
##     ## T:   threshold position
##     ## C:   cost
##     ## h2l: ratio of large effect variance to total genetic variance (not currently used)

##     ## bt <- bs

##     if(recover.flag) recover()


##     if(LL.soln)
##         stopifnot(!is.null(var.ratio))

##     if((!is.null(r2n))&(!is.null(Ve))){
##         stop('either r2n or Ve must be NULL')
##     }


##     ## solve for ft
##     ys <- log((1+bs)/(1-bs))
##     ft <- ys*(4*Ne*C)^-1  #originally ys*(2*Ne*C)^-1

##     ## get small effect additive genetic variance
##     Vas <- 8*Ne*u*L.init*as^2*bs/ys  ## originally 4*Ne*u*Ls*as^2*bs/ys
##     if(is.null(Ve))
##         Ve <- Vas*(1-r2n)/r2n
##     norm.sd <- sqrt(Vas+Ve)

##     ## initial guess for deltal
##     init.deltal <- 1##al*ft



##     ## compute selection coefficient
##     init.sl <- init.deltal*C
##     init.yl <- 4*Ne*init.sl  #originally 2*Ne*init.sl

##     if(init.yl<3)
##         warning('initial scaled selection coefficient is too small')


##     ## mean number of large effect alleles per individual
##     init.mean.nl <- 2*Ll*u/init.sl

##     ## solve for tstart conditional on initial guess of deltal

##     init.soln <- uniroot(
##         f=function(THR){
##             ftild <- dPoisConv(
##                 t=THR,
##                 lambda=init.mean.nl,
##                 norm.sd=norm.sd,
##                 alphal=al
##             )
##             ft-ftild
##         },
##         interval=c(0,12*norm.sd)
##     )
##     init.tstar <- init.soln$root


##     init.Ll <- Ll


##     ## get 2d solution
##     soln <- nleqslv(
##         x=c(0,init.tstar),
##         fn=function(X){
##             deltal <- 1/(1+exp(X[1]))
##             tstar <- X[2]

##             my.s <- deltal*C
##             my.y <- 4*Ne*my.s ## originally 2*Ne*my.s
##             mean.nl <- 2*Ll*u/my.s
##             my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##             prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##             risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)
##             deltal.tild <- risk - prot
##             ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##             Val <- 2*al*mean.nl
##             diff.one <- deltal-deltal.tild
##             diff.two <- ft-ft.tild
##             ##diff.three <- Val-Vas
##             ##diff.three <- 4*al^2*mean.nl - Vas*h2l/(1-h2l)
##             c(diff.one,diff.two)
##         },
##         control=list(scalex=c(1,1/init.tstar),maxit=400)
##     )
##     trans.deltal <- soln$x[1]
##     deltal <- 1/(1+exp(trans.deltal))
##     tstar <- soln$x[2]



##     if(LL.soln){

##         ## added the normal jitters because giving the solution as the initial condition led
##         ## it to somehow find other pathological solutions for some reason. Confused..
##         init.trans.deltal <- trans.deltal + rnorm(1,0,abs(trans.deltal)/1000)
##         if(!is.null(last.tstar)){
##           new.init.tstar <- last.tstar - abs(rnorm(1,0,tstar/1000))
##         } else{
##           new.init.tstar <- soln$x[2] + rnorm(1,0,tstar/1000)
##         }
##         init.Ll <- Ll + rnorm(1,0,Ll/1000)
##         init.Ls <- Ls
##         init.trans.bs <- log((1-bs)/bs)


##         j <- 1
##         ## get 4d solution
##         soln <- nleqslv(
##             x=c(init.trans.deltal,new.init.tstar,init.Ll,init.trans.bs,init.Ls),
##             fn=function(X){
##                 #recover()
##                 j <- j+1
##                 ## inputs
##                 deltal <- 1/(1+exp(X[1]))
##                 tstar <- X[2]
##                 Ll <- X[3]
##                 ##Ls <- L-Ll
##                 bs <- 1/(1+exp(X[4]))
##                 Ls <- X[5]
##                 L  <- Ll + Ls

##                 ## small effect stuff
##                 ys <- log((1+bs)/(1-bs))
##                 Vas <- 8*Ne*u*Ls*as^2*bs/ys # originall 4*Ne*u*Ls*as^2*bs/ys
##                 if(is.null(Ve))
##                     Ve <- Vas*(1-r2n)/r2n
##                 norm.sd <- sqrt(Vas + Ve)

##                 ## large effect stuff
##                 my.s <- deltal*C
##                 my.y <- 4*Ne*my.s #originally 2*Ne*my.s
##                 mean.nl <- 2*Ll*u/my.s
##                 my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##                 prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##                 risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)

##                 ## global stuff
##                 meana <- (as*Ls + al*Ll)/L
##                 maxg <- 2*meana*L
##                 prev <- as.numeric(pPoisConv(tstar,mean.nl,sqrt(Vas+Ve),alphal=al))
##                 risk.var <- (prev*(1-prev))

##                 ## differences
##                 deltal.tild <- risk - prot
##                 diff.one <- deltal-deltal.tild

##                 ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
##                 diff.two <- ft-ft.tild

##                 if(equalize.observed.vars){
##                     Vol <- 2*deltal^2*mean.nl
##                     Vos <- Vas*ft^2
##                     diff.three <- (Vol-var.ratio*Vos)/risk.var
##                 } else {
##                     Val <- 2*al^2*mean.nl
##                     diff.three <- Val-var.ratio*Vas
##                 }

##                 bt.tild <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)
##                 diff.four <- bt - bt.tild

##                 diff.five <- L*meana - Lmeana

##                 return(c(diff.one,diff.two,diff.three,diff.four,diff.five))
##             },
##             control=list(scalex=c(1,1/new.init.tstar,1/init.Ll,1,1/init.Ls),maxit=800,allowSingular=FALSE)
##         )
##         if(recover.flag) recover()
##         trans.deltal <- soln$x[1]
##         deltal <- 1/(1+exp(trans.deltal))
##         tstar <- soln$x[2]
##         Ll <- soln$x[3]
##         bs <- 1/(1+exp(soln$x[4]))
##         Ls <- soln$x[5]

##         ys <- log((1+bs)/(1-bs))
##         ft <- ys*(4*Ne*C)^-1 #originally ys*(2*Ne*C)^-1

##         ## get small effect additive genetic variance
##         Vas <- 8*Ne*u*Ls*as^2*bs/ys  #originally 4*Ne*u*Ls*as^2*bs/ys
##         if(is.null(Ve))
##             Ve <- Vas*(1-r2n)/r2n
##         norm.sd <- sqrt(Vas+Ve)

##     }

##     L <- Ls+Ll
##     meana <- (as*Ls + al*Ll)/L
##     deltas <- as*ft
##     mean.nl <- 2*Ll*u/(deltal*C)
##     Val <- 2*al^2*mean.nl
##     Vg <- Vas + Val
##     Vt <- Vg + Ve
##     as.std  <- as / sqrt( Vt )
##     al.std  <- al / sqrt( Vt )
##     ft.std <- ft * sqrt ( Vt )

##     ## risk scale variances
##     Vos <- Vas*ft^2
##     Vol <- 2*deltal^2*mean.nl
##     Vrg <- Vos + Vol

##     ## bulk stuff
##     prev <- as.numeric(pPoisConv(tstar,mean.nl,sqrt(Vas+Ve),alphal=al))
##     my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
##     ldens <- dpois(my.range,mean.nl) ## density on large effect liability
##     pdil <- pnorm(tstar, al*my.range,norm.sd,lower.tail=FALSE) ## prev among inds with i large effect alleles
##     pidl <- ldens*pdil/prev ## density on number of large effect alleles conditional on having disease


##     ## does mut-sel balance actually hold??
##     mutl <- 2*Ll*al*u
##     muts <- 2*Ls*as*u*bs
##     mutall <- muts+mutl

##     sell <- mean.nl*al*deltal*C
##     sels <- Vas*ft*C
##     selall <- sell+sels
##     stopifnot(abs((mutall-selall)/(mutall+selall))<1e-8)


##     tol <- 1e-5
##     min.gl <- uniroot(
##         function(X)
##             tol-(1-pPoisConv(X,mean.nl,norm.sd,alphal=al)),
##         interval=c(-10*tstar,10*tstar)
##     )$root
##     max.gl <- uniroot(
##         function(X)
##             tol-pPoisConv(X,mean.nl,norm.sd,alphal=al),
##         interval=c(-10*tstar,10*tstar)
##     )$root
##     seq.li <- seq(min.gl,max.gl,length.out=1000)
##     li.dense <- sapply(seq.li,function(G) dPoisConv(G,mean.nl,norm.sd,alphal=al))






##     ## heritabilies on liability scale
##     h2s <- Vas/(Val+Vas+Ve)
##     h2l <- Val/(Val+Vas+Ve)
##     h2 <- (Vas+Val)/(Val+Vas+Ve)



##     ## liability scale h2 estimation via normal approx
##     phit <- dnorm(qnorm(1-prev))
##     risk.var <- (prev*(1-prev))

##     ## observed scale
##     h2os <- Vos/risk.var
##     h2ol <- Vol/risk.var
##     h2o <- Vrg/risk.var
##     prs <- Vos/(Vos+Vol)

##     ## naive estimates on liability scale
##     h2s.est <- h2os*risk.var/Vt*phit^-2
##     h2l.est <- h2ol*risk.var/Vt*phit^-2
##     h2all.est <- h2o*risk.var/Vt*phit^-2

##     implied.al.norm <- deltal/phit
##     implied.al.true <- deltal/ft


##     ## other stuff
##     #L <- Ll+Ls

##     meana <- (as*Ls + al*Ll)/L
##     maxg <- 2*meana*L
##     bt <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)

##     return(
##         list(
##             soln=soln,
##             ft=ft,
##             ft.std=ft.std,
##             phit=phit,
##             bs=bs,
##             bt=bt,
##             as=as,
##             al=al,
##             as.std = as.std,
##             al.std = al.std,
##             meana=meana,
##             maxg=maxg,
##             bt=bt,
##             my.range=my.range,
##             deltas=deltas,
##             deltal=deltal,
##             ss=deltas*C,
##             sl=deltal*C,
##             ys=4*Ne*deltas*C,  #originally 2*Ne*deltas*C
##             yl=4*Ne*deltal*C,  #originally 2*Ne*deltal*C
##             Ne=Ne,
##             u=u,
##             C=C,
##             ors=(prev+deltas)/prev*(1-prev)/(1-prev-deltas),
##             orl=(prev+deltal)/prev*(1-prev)/(1-prev-deltal),
##             tstar=tstar,
##             mean.nl=mean.nl,
##             prev=prev,
##             Vas=Vas,
##             Val=Val,
##             Vg=Vg,
##             Ve=Ve,
##             norm.sd=norm.sd,
##             Vos=Vos,
##             Vol=Vol,
##             h2s=h2s,
##             h2l=h2l,
##             h2=h2,
##             h2os=h2os,
##             h2ol=h2ol,
##             h2o=h2o,
##             h2s.est=h2s.est, ## small effect h2 you would infer by naive application of ltm methods
##             h2l.est=h2l.est, ## large effect h2 you would infer by naive application of ltm methods
##             h2all.est=h2all.est, ## total h2 you would infer by naive application of ltm methods
##             pgs.est=h2s.est/h2all.est,
##             implied.al.norm=implied.al.norm,
##             implied.al.true=implied.al.true,
##             pgs=Vas/(Val+Vas),
##             pgl=Val/(Val+Vas),
##             Vos=Vos,
##             Vol=Vol,
##             pos=Vos/(Vos+Vol),
##             pol=Vol/(Vos+Vol),
##             ldens=ldens,
##             pdil=pdil,
##             pidl=pidl,
##             Ls=Ls,
##             Ll=Ll,
##             li.dense=li.dense
##         )
##     )
## }
