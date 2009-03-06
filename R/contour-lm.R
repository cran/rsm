### reconstructs the data set used in a linear model, and
### returns it as a data.frame

model.data = function (lmobj, lhs=FALSE) {
    form = lmobj$call$formula
    if (is.name(form)) {
        lmobj$call$data = form
        form = formula(lmobj)
    }
    if (lhs)
        nm = all.vars(form)
    else
        nm = all.vars(form[[3]])
    form = as.formula(paste("~", paste(nm, collapse = "+")))
    model.frame(form, eval(lmobj$call$data), subset = eval(lmobj$call$subset))
}



### contour plot(s) for a lm
contour.lm = function(x, form, at, bounds, zlim,
    image=TRUE, img.col=terrain.colors(50), ...) 
{
  lmobj = x   # generic wants it named 'x', I don't!
  
  
  data = model.data(lmobj)
  
  # make list of formulas if not already
  if(inherits(form,"formula")) {
    if (length(form)==2) { # rhs only
      vars = all.vars(form[[2]])
      n = length(vars)
      if (n < 2) stop("Need at least two variables")
      form = list()
      elt = 1
      for (i in 1:(n-1))
        for (j in (i+1):n) {
          form[[elt]] = formula(paste(vars[j],vars[i],sep="~"))
          elt = elt + 1
        }
    }
    else {
      yvars = all.vars(form[[2]])
      xvars = all.vars(form[[3]])
      form = list()
      elt = 1
      for (i in 1:length(xvars))
        for (j in 1:length(yvars)) {
          form[[elt]] = formula(paste(yvars[j],xvars[i],sep="~"))
          elt = elt + 1
        }
    }
  }

  # gather 'at' info
  tmp = lapply(data, function(var) {
    if (is.factor(var)) factor(levels(var)[1], levels=levels(var))
    else mean(var)
  })
  if (!missing(at))
    for (nm in names(at))
      tmp[[nm]] = at[[nm]]
  at = tmp
  
  # gather 'bounds' info -- elts can be vectors of length 2, 3, or n
  tmp = lapply(data, function(x) if (is.numeric(x)) range(x))
  if (!missing(bounds))
    for (nm in names(bounds))
      if (length(bounds[[nm]]) > 1)
        tmp[[nm]] = bounds[[nm]]
  bounds = lapply(tmp, function(x) {
    if (length(x) == 2) seq(x[1], x[2], length=26)
    else if (length(x) == 3) seq(x[1], x[2], length=x[3])
    else x
  })
  
  
  ### Accumulate the z values
  plot.data = list()
  z.rng = NULL
  for (i in 1:length(form)) {
    AT = at
    v = all.vars(form[[i]])
    if (length(unique(v)) == 1) next
    y = AT[[v[1]]] = bounds[[v[1]]]
    x = AT[[v[2]]] = bounds[[v[2]]]
    newdata = do.call(expand.grid, AT)
    ord = order(newdata[[v[1]]], newdata[[v[2]]])
    newdata = newdata[ord, ]
    ###z = matrix (predict (lmobj, newdata = newdata), nrow = length(x))
    z = plot.data[[i]] =  predict (lmobj, newdata = newdata)
    z.rng = range (c (z.rng, z))
  }
  
  if (missing (zlim)) zlim = z.rng
  
  ### Do then plots with a common image scale
  for (i in 1:length(form)) {
    v = all.vars(form[[i]])
    x = bounds[[v[2]]]
    y = bounds[[v[1]]]
    z = matrix (plot.data[[i]], nrow = length (x))
    if (image) {
      image(x, y, z, col=img.col, xlab = v[2], ylab = v[1], zlim = zlim, ...)
      contour(x, y, z, add=TRUE, ...)
    }
    else {
      contour(x, y, z, ,xlab = v[1], ylab = v[2], ...)
    }
  }
}


