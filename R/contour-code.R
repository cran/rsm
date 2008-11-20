### reconstructs the data set used in a linear model, and
### returns it as a data.frame

model.data = function(lmobj) {
  find.names = function(obj) {
    if (is.call(obj)) {
      obj[[1]] = ""
      sapply(obj, find.names)
    }
    else {
      if (is.name(obj)) as.character(obj)  
      else NULL
    }
  }
  vars = attr(lmobj$terms, "variables")
  vars[[1]] = ""
  varnames = unlist(lapply(vars, find.names), recursive=TRUE)
  form = as.formula(paste(c("~-1",varnames), collapse="+"))
  model.frame (form, eval(lmobj$call$data), subset=eval(lmobj$call$subset))
}



### Uses 1st 2 vars in list for x and y, others allowed at fixed values
contour.lm = function(x, varlist, image=TRUE, img.col=terrain.colors(50), ...) {
  lmobj = x   # generic wants it named 'x', I don't!
  data = model.data(lmobj)
  fcn = function(var) {
    if (is.factor(var)) factor(levels(var)[1], levels=levels(var))
    else mean(var)
  }
  typical = lapply(data, fcn)
  nm = names(varlist)
  axv = list()
  for (n in nm[1:2]) {
    x = varlist[[n]]
    if (is.null(x)) x = 26
    if (length(x) == 1) {
      r = range(data[[n]])
      x = seq(r[1], r[2], length=x)
    }
    else if (length(x) == 2)
      x = seq(x[1], x[2], length=26)
    else if (length(x) == 3)
      x = seq(x[1], x[2], length=x[3])
    axv[[n]] = typical[[n]] = x
  }
  newdata = do.call(expand.grid, typical)
  ord = order(newdata[[nm[2]]], newdata[[nm[1]]])
  newdata = newdata[ord, ]
  if (length(nm) > 2) 
    for (n in nm[3:length(nm)]) {
      if (is.factor(data[[n]]))
        newdata[[n]] = factor(varlist[[n]], levels=levels(data[[n]]))
      else
        newdata[[n]] = varlist[[n]]
    }
  z = matrix (predict (lmobj, newdata = newdata), nrow = length(axv[[1]]))

  if (image) {
    image(axv[[1]], axv[[2]], z, col=img.col, ,xlab = nm[1], ylab = nm[2], ...)
    contour(axv[[1]], axv[[2]], z, add=TRUE, ...)
  }
  else {
    contour(axv[[1]], axv[[2]], z, ,xlab = nm[1], ylab = nm[2], ...)
  }
}
