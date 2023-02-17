.hello <- function(){
  #' @export
  message("Hello")
}


# classes
measures <- function(...){
  #' @title Data structure for repeated measures
  #' @export
  call = match.call(expand.dots = TRUE)
  charargs = as.character(call)[-1]
  lens = sapply(1:length(charargs), function(i) length(eval(as.name(charargs[i]), envir = NULL)))
  if(length(unique(lens))==1){
    ret = cbind(...)
    class(ret) <- c("multimeasure", "matrix")
  }
  if(length(unique(lens))>1){
    ret = list(...)
    class(ret) <- c("multimeasure", "list")
  }
ret
}

as.measures <- function(x,...){
  #' @title Data structure for repeated measures
  #' @export
  if(inherits(x, "list")) class(x) <- c("multimeasure", "list")
  if(inherits(x, "matrix")) class(x) <- c("multimeasure", "matrix")
  x
}

dgm2 <- function(n=50, mu=50, sd=5, errsd=1, f2=m2 ~ m1 + m1^2 + err){
  #' @title Generating measures data
  #' @export
  e = new.env()
  #fun = f2[[3]]
  #print(fun)
  e$m1 <- rnorm(n, mu, sd)
  e$err <- rnorm(n, 0, errsd)
  #ff = eval(f2, env=e)
  fun = f2[[3]]
  e$m2 <- eval(fun, e)
  m1= e$m1
  m2 = e$m2
  measures(m1, m2)
}


# generics
bland_altman = function(object, ...){
  #' @title Bland altman measure of reliability
  #' @export
  UseMethod("bland_altman")
}

correlation = function(object, ...){
  #' @title Pearson correlation coefficient
  #' @export
  UseMethod("correlation")
}

regression = function(object, ...){
  #' @title Linear model coefficient
  #' @export
  UseMethod("regression")
}

rankcorrelation = function(object, ...){
  #' @title Spearman rank correlation coefficient
  #' @export
  UseMethod("rankcorrelation")
}

icc = function(object, ...){
  #' @title Intraclass correlation coefficient
  #' @export
  UseMethod("icc")
}

.makelongdata <- function(object, ...){
  #' @title Making measures data into a dataframe
  #' @export
  UseMethod(".makelongdata")
}

.bland_altman2 <- function(object, refcol=NA, ...){
  #' @export
  if(is.null(dim(object)) || dim(object)[2]!=2 ) stop("only implemented for 2 column matrixes")
  if(is.na(refcol)){
    x = apply(object, 1, mean)
  }
  if(!is.na(refcol)){
    x = object[,refcol]
  }
  n = length(x)
  y = (object[,1]-object[,2])
  u = mean(y)
  s = sd(y)
  tcrit = qt(.975, n-1)
  ret = list(x=x,y=y, u=u, s=s, sm=s/sqrt(n), tcrit=tcrit, measures = object)
  class(ret) <- "blandaltman"
  attr(ret, "refcol") <- refcol
  ret
}

.bland_altman3 <- function(object, refcol=NA, ...){
  #' @export
  if(is.na(refcol)){
    x = apply(object, 1, mean)
  }
  if(!is.na(refcol)){
    x = object[,refcol]
  }
  n = length(x)
  y = apply(object, 1, sd)
  u = mean(y)
  s = sd(y)
  tcrit = qt(.975, n-1)
  ret = list(x=x,y=y, u=u, s=s, sm=s/sqrt(n), tcrit=tcrit, measures = object)
  class(ret) <- "blandaltman"
  attr(ret, "refcol") <- refcol
  ret
}

bland_altman.multimeasure <- function(object, refcol=NA, method=NULL, ...){
  #' @export
  if(is.null(dim(object))) stop("only implemented for 2+ column matrixes")
  if(!is.null(method) && method<2) method=dim(object)[2]
  if(!is.null(method) && method>=2) {
    method==max(method, 3)
    }
  if(is.null(method)) method = min(dim(object)[2], 3)
  if(method==3) ret = .bland_altman3(object, refcol=refcol, ...)
  if(method==2) ret = .bland_altman2(object, refcol=refcol, ...)
  ret
}


correlation.multimeasure <- function(object, ...){
  #' @export
  if(is.null(dim(object)) || dim(object)[2]!=2 ) stop("only implemented for 2 column matrixes")
  x = object[,1]
  #n = length(x)
  y = object[,2]
  num = sum((x - mean(x))*(y-mean(y)))
  den = sqrt(sum((x - mean(x))^2))*sqrt(sum((y-mean(y))^2))
  ret = num/den
  class(ret) <- "correlation"
  attr(ret, 'method') = "Pearson"
  ret
}

.myrank <- function(X, tol = 1e-07){
  #qr(zapsmall(X, digits = -log10(tol) + 5), tol = tol, LAPACK = FALSE)$rank
  rank(X, na.last="keep")
}

rankcorrelation.multimeasure <- function(object, ...){
  #' @export
  if(is.null(dim(object)) || dim(object)[2]!=2 ) stop("only implemented for 2 column matrixes")
  #object[,] = apply(object,2, order)
  object[,] = apply(object,2, .myrank)
  ret = correlation(object, ...)
  attr(ret, 'method') = "Spearman"
  ret
}

regression.multimeasure <- function(object, ...){
  #' @export
  if(is.null(dim(object)) || dim(object)[2]!=2 ) stop("only implemented for 2 column matrixes")
  x = object[,1]
  #n = length(x)
  y = object[,2]
  xy = x-y
  #xy= c(x,y)
  #method = c(rep(1, length(x)), rep(0, length(y)))
  m = lm(y~x)
  m2 = lm(x~y)
  m3 = lm(xy~1)
  ret = as.numeric(c(m$coefficients, sd(m$residuals),
                   m2$coefficients, sd(m2$residuals),
                   m3$coefficients, sd(m3$residuals)))
  class(ret) <- "regression"
  attr(ret, 'method') = "Linear model"
  ret
}


.fishericc <- function(object, ...){
  if(is.null(dim(object)) || dim(object)[2]==1 ) stop("only implemented for 2+ column matrixes")
  nvar = ncol(object)
  n = nrow(object)
  #xy = as.numeric(object)
  #xybar = mean(xy) # also 0.5*(mean of x+y)
  xybar = mean(object) # also 0.5*(mean of x+y)
  nump = object - xybar
  combs = combn(nvar,2)
  num = 0
  for(comb in 1:ncol(combs)){
    num = num + sum(nump[,combs[1,comb]]*nump[,combs[2,comb]])
  }
  dencols = (object - xybar)^2
  s2 = sum(apply(dencols, 2, sum))/(nvar*n)
  #s2 = (sum((x - xybar)^2) + sum((y-xybar)^2))/(2*n)  # original
  ret = num/(s2*n*ncol(combs))                        # original
  class(ret) <- "icc"
  attr(ret, 'method') = "Fisher"
  ret
}

.makelongdata.multimeasure <- function(object, ...){
  #' @export
  if(inherits(object, "matrix")){
    ids = 1:dim(object)[1]
    id=rep(ids, ncol(object))
    val=as.numeric(object)
    valtype = rep(1:ncol(object), each=nrow(object))
  }
  if(inherits(object, "list")){
    id = c()
    val = c()
    valtype = c()
    vt = 1
    for(item in object){
      id = c(id, 1:length(item))
      val = c(val, item)
      valtype = c(valtype, vt)
      vt = vt+1
    }

  }
  data.frame(id = id, val=val, valtype=valtype)
}

.randefficc <- function(object, ...){
  #' @import lme4
  df = .makelongdata(object)
  # one way random effect model, method 1 of mcgraw, wong
  fit = lmer(val ~ 1 + (1 | id), data=df)
  vc = VarCorr(fit)
  saj = attr(vc$id, "stddev")
  sr = attr(vc, "sc")
  ret = as.numeric(saj^2/(saj^2 + sr^2))
  class(ret) <- "icc"
  attr(ret, 'method') = "Random effects (method 1)"
  ret
}

.randeffrowicc <- function(object, ...){
  #' @import lme4
  df = .makelongdata(object)
  # two way random effect model, method 2 of mcgraw, wong
  fit = lmer(val ~ 1 + (1 | id) + (1 | valtype), data=df)
  vc = VarCorr(fit)
  sai = attr(vc$valtype, "stddev")
  saj = attr(vc$id, "stddev")
  sr = attr(vc, "sc")
  ret = as.numeric(saj^2/(saj^2 + sr^2))
  class(ret) <- "icc"
  attr(ret, 'method') = "Random effects (method 2)"
  ret
}

icc.multimeasure <- function(object,
                             method="fisher",
                             ...){
  #' @title ICC calculation for a
  #' @param object multimeasure object
  #' @param method "fisher" "randeff" or "randeff2"
  #' @export
  if(inherits(object, "matrix") & method=="fisher") return(.fishericc(object))
  if(inherits(object, "list") || method=="randeff" || method=="randeff1") return(.randefficc(object))
  if(inherits(object, "list") || method=="randeff2") return(.randeffrowicc(object))
}



as.data.frame.multimeasure <- function(object,...){
  #' @export
  as.data.frame.matrix(object,...)
}

print.blandaltman <- function(x, ...){
  #' @export
  cat(paste0("mu = ",x$u, "; s^2 = ",x$s^2, "\n"))
}

print.correlation <- function(x, ...){
  #' @export
  cat(paste0(attr(x, "method"), " correlation = ",x, "\n"))
}

print.regression <- function(x, ...){
  #' @export
  cat(paste0(attr(x, "method"),
             "\nM2 ~ intercept + slope*M1 + error:\nintercept = ",x[1],
             "\nslope = ", x[2],
             "\nerror sd  = ", x[3],
             "\n\nM1 ~ intercept + slope*M2 + error:\nintercept = ",x[4],
             "\nslope = ", x[5],
             "\nerror sd  = ", x[6],
             "\n\nM1-M2 ~ intercept + error:\nintercept = ",x[7],
             "\nerror sd  = ", x[8],
             "\n"))
}


print.icc <- function(x, ...){
  #' @export
  cat(paste0(attr(x, "method"), " icc = ",x, "\n"))
}

plot.blandaltman <- function(x,
                             data=TRUE,
                             bias=TRUE,
                             error=TRUE,
                             meanerror=TRUE,
                             regline=TRUE,
                             ideal=TRUE,
                             alpha=1.0,
                             m1nm = NULL,
                             m2nm = "Difference between measures",
                             ...){
 #' @import ggplot2
 #' @export
 # base
 lim= max(abs(min(x$y)), abs(max(x$y)), abs(x$u+c(1.96, -1.96)*x$s))
 p <- ggplot() +
   theme_classic() +
   scale_color_grey() +
   scale_y_continuous(name=m2nm, limits = c(-lim, lim)) +
   guides(shape = "none", size = "none")
 if(ideal) p <- p + geom_hline(aes(yintercept=0, linetype="Ideal difference"))

 if(is.na(attr(x, "refcol"))){
   p <- p + scale_x_continuous(name="Mean across measures")
 }
 if(!is.na(attr(x, "refcol"))){
   p <- p + scale_x_continuous(name=ifelse(attr(x, "refcol")==1,
                                           "Measure 1",
                                           "Measure 2"))
 }
 if(!is.null(m1nm)){
   p <- p + scale_x_continuous(name=m1nm)
 }
 # the data
 if(data) p <- p +  geom_point(aes(x=x,y=y), data=data.frame(x=x$x, y=x$y), alpha=alpha)

 if(bias | error | regline | ideal) p <- p + scale_linetype_discrete(name="")

 if(bias){
   # overall error
   p <- p +
     geom_hline(aes(yintercept=u, linetype="Mean difference (bias)"),
                data=data.frame(u=x$u))
 }
 if(error){
   p <- p +
     geom_hline(aes(yintercept=u, linetype="95% limits of agreement"),
                data=data.frame(u=x$u+c(1.96, -1.96)*x$s))
 }
 if(regline){
   # regline
   p <- p +
     geom_smooth(aes(x=x,y=y, linetype="Regression fit"), se=FALSE, method="lm",
                 data=data.frame(x=x$x, y=x$y), color="black", size=0.5)

 }
 if(meanerror){
   # mean error
   p <- p +
     geom_ribbon(aes(ymin=umin, ymax=umax, x=x, fill="95% bias confidence intervals"),
                 alpha=0.2, data=data.frame(x=c(min(x$x), max(x$x)),
                                 umin = rep(x$u-x$tcrit*x$sm, 2),
                                 umax = rep(x$u+x$tcrit*x$sm, 2))) +
     scale_fill_grey(name="")
 }
 print(p)
}

plot.multimeasure <- function(x,
                              data=TRUE,
                              best=TRUE,
                              smooth=TRUE,
                              m1nm = "Measure 1",
                              m2nm = "Measure 2",
                              alpha=1.0,
                              ...){
  #' @import ggplot2
  #' @export
  # base
  upl = max(x)
  lowl = min(x)
  pdat = as.data.frame(x)
  p = ncol(pdat)
  names(pdat)[1:2] <- paste0("m", 1:p)
  p <- ggplot(data=pdat, aes(x=m1, y=m2)) +
    theme_classic() +
    labs(y=m2nm, x=m1nm) +
    coord_fixed(ratio=1) +
    lims(x=c(lowl, upl),y =c(lowl, upl)) +
    scale_color_grey(name="")
  if(data) p <- p + geom_point(aes(color="Data"), alpha=alpha)
  if(smooth) p <- p + geom_smooth(aes(linetype="Smooth fit"), se=FALSE, color="gray50")
  if(best) p <- p + geom_abline(aes(linetype="Line of agreement", intercept = 0,slope=1), color="gray50")
  if(smooth | best) p <- p + scale_linetype_discrete(name="") +
    guides(fill = "none",  shape = "none", size = "none")
  print(p)
}