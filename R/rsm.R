### Functions to facilitate response-surface analysis

# First-order model
FO = function(...) {
  fo = sapply(list(...), I)
  if (is.null(nrow(fo))) fo = matrix(fo, nrow=1)
  fo
}

# Pure quadratic
PQ = function(...) {
  X = FO(...)^2
  nm = dimnames(X)[[2]]
  if (is.null(nm)) nm = 1:ncol(X)
  dimnames(X) = list(NULL, paste(nm,"2",sep="^"))
  X
}

# Two-way interactions
TWI = function(...) {
  fo = FO(...)
  k = ncol(fo)
  fon = dimnames(fo)[[2]]
  if (is.null(fon)) fon=1:k
  X = matrix(0, nrow=nrow(fo), ncol=k*(k-1)/2)
  nm = rep("", k*(k-1)/2)
  col = 1
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      X[, col] = fo[ ,i] * fo[ ,j]
      nm[col] = paste(fon[i],fon[j],sep="")
      col = col+1
    }
  }
  dimnames(X) = list(NULL,nm)
  X
}

# Second-order model.  But in rsm(), this will get replaced by FO()+TWI()+PQ()
SO = function(...)
  X = cbind(FO(...), TWI(...), PQ(...))


# Pure-error model
PE = function(...)
  factor(paste(...))



# Fit a response-surface model
rsm = function (..., data) {
    CALL = match.call(lm)
    CALL[[1]] = as.name("lm")
    oc = as.character(deparse(CALL$formula))
    nc = sub("SO\\(([a-zA-Z0-9, ._]+)\\)", "FO\\(\\1\\) + TWI\\(\\1\\) + PQ\\(\\1\\)", oc)
  # no comma -> only 1 var -> no TWI ...
    nc = sub("TWI\\([a-zA-Z0-9 ._]+\\)", "", nc)
    CALL$formula = formula(nc)
    LM = eval(CALL)
    LM$call[[1]] = as.name("rsm")
    LM$call$formula = formula(oc)
    
    nm = names(LM$coef)
    i.fo = grep("FO\\(", nm)
    if (length(i.fo) == 0) {
        warning("No FO() terms in model; cannot use RSM methods\nAn 'lm' object has been returned.")
        return(LM)
    }
    k = length(i.fo)
    LM$b = LM$coef[i.fo]
    LM$order = 1
    foterm = as.list(LM$terms[LM$assign[min(i.fo)]][[3]])
    fonm = names(LM$b) = sapply(foterm, as.character)[-1]
    LM$labels = list(FO=list(idx=i.fo, lab=fonm))
    names(LM$coef)[i.fo] = LM$labels

    i.twi = grep("TWI\\(", nm)
    if ((k > 1) & (length(i.twi) == k*(k-1)/2)) {
        btwi = LM$coef[i.twi]
        LM$order = 1.5
        LM$B = diag(rep(0,k))
        col = 1
        twi.lab = rep("", length(i.twi))
        for (i in 1:(k-1)) {
          rng = 1:(k-i)
          LM$B[i,rng+i] = LM$B[rng+i,i] = 0.5 * btwi[col+rng-1]
          twi.lab[col+rng-1] = paste(fonm[i], fonm[rng+i], sep=":")
          col = col + k-i
        }
        dimnames(LM$B) = list(fonm, fonm)
        LM$labels$TWI = list(idx=i.twi, lab=twi.lab)
    }
    else if (length(i.twi) > 0)
        warning(paste("TWI() term not usable because it has", length(i.twi), "d.f. instead of", k*(k-1)/2))
        
    i.pq = grep("PQ\\(", nm)
    if (length(i.pq) == k) {
        LM$order = 2
        if (is.null(LM$B)) {
            if (k > 1)  LM$B = diag(LM$coef[i.pq])
            else        LM$B = matrix(LM$coef[i.pq], nrow=1)
            dimnames(LM$B) = list(fonm, fonm)
        }
        else
            diag(LM$B) = LM$coef[i.pq]
        LM$labels$PQ = list(idx=i.pq, lab=paste(fonm,2,sep="^"))
    }
    else if (length(i.pq) > 0)
        warning(paste("PQ() term not usable because it has", length(i.pq), "d.f. instead of", k))
        
    if (!missing(data)) 
        if (inherits(data, "coded.data")) 
            LM$coding = attr(data, "coding")
    class(LM) = c("rsm", "lm")
    LM
}

# do a lack-of-fit test
loftest = function (object) {
    cl = match.call(lm, call = object$call)
    cl[[1]] = as.name("lm")
    pieces = as.character(object$call$formula)
    pieces[3] = sub("(FO)|(SO)", "PE", pieces[3])
    cl$formula = formula(paste(pieces[2], "~", pieces[3]))
    lof = anova(object, eval(cl))
    df = c(lof[1,1], lof[2,3], lof[2,1])
    ss = c(lof[1,2], lof[2,4], lof[2,2])
    ans = data.frame (
        df, ss, ss/df, c(NA, lof[2,5], NA), c(NA, lof[2,6], NA),
        row.names = c("Model residual", "Lack of fit", "Pure error"))
    names(ans) = c("Df","Sum Sq","Mean Sq", "F value","Pr(>F)")
    class(ans) = class(lof)
    ans
}

# Summary method
summary.rsm = function (object, ...) {
    SUM = summary.lm(object, ...)
    if (object$order > 0) {
        for (lst in object$labels)
            row.names(SUM$coefficients)[lst$idx] = lst$lab
        if (object$order > 1) {
            SUM$canonical = list(xs = -0.5 * solve(object$B, object$b), 
                eigen = eigen(object$B))
        }
        else SUM$sa = object$b/sqrt(sum(object$b^2))
        SUM$lof = rbind(anova(object), loftest(object)[-1,])
        SUM$coding = object$coding
        class(SUM) = c("summary.rsm", "summary.lm")
    }
    SUM
}

# Print method for summary
print.summary.rsm = function(x,...) {
  getS3method("print", "summary.lm") (x, ...)
  print(x$lof, signif.stars=FALSE, ...)
  cat("\n")
  can = x$canonical
  if (!is.null(can)) {
    cat("Stationary point of response surface:\n")
    print(can$xs)
    if(!is.null(x$coding)) {
      cat("\nStationary point in original units:\n")
      print (code2val (can$xs, x$coding))
    }
    cat("\nEigenanalysis:\n")
    print(can$eigen)
  }
  else {
    cat("Direction of steepest ascent (at radius 1):\n")
    print(x$sa)
    cat("\nCorresponding increment in original units:\n")
    temp = code2val (rbind(x$sa, 2*x$sa), x$coding)
    print (temp[2,] - temp[1,])
  }
  cat("\n")
}


# Steepest ascent (and ridge analysis)
steepest = function (object, dist=seq(0,5,by=.5), descent=FALSE) {
    goal = ifelse(descent, "descent", "ascent")
    dist = abs (dist)
    if (is.null(object$B)) {
        d = object$b / sqrt (sum (object$b^2))
        if (descent) d = -d
        path = t(sapply(dist, function(x) d*x))
        cat(paste("Linear path of steepest", goal, "\n"))
    }
    else {
        iden = diag (rep (1, length(object$b)))
        rng = range (eigen (object$B) $values)
          
        soln = function (gam) {
            -0.5 * solve (object$B - gam*iden, object$b)
        }
        deldist = function (gam, d) {
            xx = soln (gam)
            sqrt (sum (xx^2)) - d
        }
        find.pt = function(d, bd) {
            if (abs(d) < .01) return (0 * object$b)
            gamma = uniroot (deldist, bd, d)$root
            soln (gamma)
        }
        incr.out = function(bd, inc, mind) {
            while (deldist(bd, mind) > 0) {
              bd = bd + inc
              inc = 2*inc
            }
            bd
        }
        
        mind = min(dist[dist>.009])
        if (descent)
          bds = c(incr.out(rng[1]-5, -2, mind), rng[1]-.001)
        else 
          bds = c(rng[2]+.001, incr.out(rng[2]+5, 2, mind))

        path = t(sapply(dist, find.pt, bds))
        cat(paste("Path of steepest", goal, "from ridge analysis:\n"))
    }
    
    path = newdata = as.data.frame (round (path, 3))
    md = model.data(object)
    for (vn in names(md))
        if (is.null(newdata[[vn]])) {
            v = md[[vn]]
            if(is.factor(v)) 
                newdata[[vn]] = factor(levels(v)[1], levels=levels(v))
                else newdata[[vn]] = mean(v)
        }
    yhat = predict(object, newdata=newdata)
    
    path[["|"]] = factor("|")
    if (!is.null(object$coding)) {
        orig = code2val(path, object$coding)
        path = cbind(path, orig)
    }
    ans = cbind(data.frame(dist=dist), path, yhat=round(yhat,3))
    ans
}

canonical.path = function(object, 
  which = ifelse(descent, length(object$b), 1),
  dist = seq(-5, 5, by=.5),
  descent = FALSE)
{
  if (!inherits(object, "rsm"))
    stop(paste(as.character(substitute(object)),"is not an 'rsm' object"))
  if (object$order == 1)
    stop("Requires a seconnd-order response surface")
  can = summary(object)$canonical
  dir = can$eigen$vectors[ , which]
  path = t(sapply(dist, function(d) can$xs + d*dir))

  path = newdata = as.data.frame(round(path, 3))
  md = model.data(object)
  for (vn in names(md)) if (is.null(newdata[[vn]])) {
    v = md[[vn]]
    if (is.factor(v)) 
      newdata[[vn]] = factor(levels(v)[1], levels = levels(v))
    else newdata[[vn]] = mean(v)
  }
  yhat = predict(object, newdata = newdata)
  path[["|"]] = factor("|")
   if (!is.null(object$coding)) {
    orig = code2val(path, object$coding)
    path = cbind(path, orig)
  }
  ans = cbind(data.frame(dist = dist), path, yhat = round(yhat, 3))
  ans
}
