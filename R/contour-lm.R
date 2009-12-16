### Reconstructs the data set used in a linear model, and
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
    xlabs, hook, plot.it=TRUE, # args added Dec 09
    image=FALSE, img.col=terrain.colors(50), ...) 
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
  vars = sort(unique(as.character(sapply(form, all.vars))))
  if (missing(xlabs) && inherits(lmobj, "rsm")) {
    forms = codings(lmobj)
    if (!is.null(forms))
      xlabs = sapply(vars, function(v) paste(as.character(forms[[v]][2:3]), collapse=" = "))
  }   
  else if (missing(xlabs)) xlabs = vars

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
  lbls = rep("", length(form))
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
    z = predict (lmobj, newdata = newdata)
    z.rng = range (c (z.rng, z))
    if (!missing(zlim)) {
      z[z > zlim[2]] = NA
      z[z < zlim[1]] = NA
    }
    vnames = c(x=v[2], y=v[1])
    labs = c(xlabs[sapply(vnames, charmatch, vars)], vnames)
    plot.data[[i]] = list(x=x, y=y, 
      z=matrix(z, nrow=length(x)), labs=labs)
    lbls[i] = paste(labs[3], labs[4], sep=" ~ ")
  }
  names(plot.data) = lbls
  
  if (missing (zlim)) zlim = z.rng
  for (i in 1:length(lbls))
    plot.data[[i]]$zlim = zlim
  
  ### If plots requested, do plots with a common image scale
  if (plot.it) for (i in 1:length(form)) {
    dat = plot.data[[i]]
    if (!missing(hook))
      if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
    if (image) {
      image(dat$x, dat$y, dat$z, col=img.col, 
        xlab = dat$labs[1], ylab = dat$labs[2], zlim = dat$zlim, ...)
      contour(dat$x, dat$y, dat$z, add=TRUE, ...)
    }
    else {
      contour(dat$x, dat$y, dat$z, 
        xlab = dat$labs[1], ylab = dat$labs[2], ...)
    }
    if (!missing(hook))
      if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
  }
  
  invisible(plot.data)
}


# Image plot for a lm
image.lm = function(x, form, at, bounds, zlim, xlabs, hook,  ...)  {
  plot.data = contour.lm(x, form, at, bounds, zlim, xlabs, plot.it=FALSE)
  for (i in 1:length(plot.data)) {
    dat = plot.data[[i]]
    if (!missing(hook))
      if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
    
    image(dat$x, dat$y, dat$z, 
        xlab = dat$labs[1], ylab = dat$labs[2], zlim = dat$zlim, ...)
    
    if (!missing(hook))
      if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
  }
  
  invisible(plot.data)
}


# Perspective plot(s) for 'lm' objects
# arg notes:
# col: facet colors; if null, default, else color pallette based on z value
# contours: if TRUE, black contours.  Can also be a list with elements
#   z="bottom" (or "top" or value), col="black", lwd=1
persp.lm = function(x, form, at, bounds, zlim, zlab, 
    xlabs, col = "white", contours=NULL, hook,
    theta = -25, phi = 20, r = 4, border=NULL, box=TRUE, ticktype = "detailed", 
    ...) 
{
  draw.cont.line = function(line) {
    if (cont.varycol) {
      cont.col = col
      if (length(col) > 1) cont.col = col[cut(c(line$level, dat$zlim), length(col))][1]
    }
    lines(trans3d(line$x, line$y, cont.z, transf),
      col=cont.col, lwd=cont.lwd)
  }
  plot.data = contour.lm(x, form, at, bounds, zlim, xlabs, plot.it=FALSE)
  transf = list()
  if (missing(zlab)) zlab = ""
  
  facet.col = col

  cont = !is.null(contours)
  if (mode(contours) == "logical") cont = contours
  cont.first = cont
  cont.z = cz = plot.data[[1]]$zlim[1]
  cont.col = 1
  cont.varycol = FALSE
  cont.lwd = 1
  if (is.character(contours)) {
    idx = charmatch(contours, c("top","bottom", "colors"), 0)
    if (idx == 1) {
      cont.first = FALSE
      cont.z = plot.data[[1]]$zlim[2]
    }
    else if (idx == 2) {}
    else if (idx == 3) {
      cont.varycol = TRUE
      if (length(col) < 2) col = rainbow(40)
    }
    else
      cont.col = contours
  }
  else if (is.list(contours)) {
    if(!is.null(contours$z)) cz = contours$z
    if (is.numeric(cz)) cont.z = cz
    else if (cz=="top") {
      cont.first = FALSE
      cont.z = plot.data[[1]]$zlim[2]
    }
    if(!is.null(contours$col)) cont.col = contours$col
    if(!is.null(contours$lwd)) cont.lwd = contours$lwd
    if(charmatch(cont.col, "colors", 0) == 1) {
      cont.varycol = TRUE
      if (length(col) < 2) col = rainbow(40)
    }
  }
  
  # Loop through the plots
  for (i in 1:length(plot.data)) {
    dat = plot.data[[i]]
    cont.lines = NULL
    if (!missing(hook))
      if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
    if (cont) cont.lines = contourLines(dat$x, dat$y, dat$z)
    if (cont && cont.first) {
      transf = persp(dat$x, dat$y, dat$z, zlim=dat$zlim, theta=theta, phi=phi, r=r, col = NA, border=NA, box=FALSE, ...)
      lapply(cont.lines, draw.cont.line)
      par(new=TRUE)
    }
    if (length(col) > 1) {
      nrz = nrow(dat$z)
      ncz = ncol(dat$z)
      zfacet = dat$z[-1,-1] + dat$z[-1,-ncz] + dat$z[-nrz,-1] + dat$z[-nrz,-ncz]
      zfacet = c(zfacet/4, dat$zlim)
      facet.col = cut(zfacet, length(col))
      facet.col = col[facet.col]
    }
    transf = persp(dat$x, dat$y, dat$z,
      xlab=dat$labs[1], ylab=dat$labs[2], zlab=zlab,
      zlim=dat$zlim, col=facet.col, border=border, box=box, theta=theta, phi=phi, r=r, ticktype=ticktype, ...)
    if (cont && !cont.first)
      lapply(cont.lines, draw.cont.line)
    if (!missing(hook))
      if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
    plot.data[[i]]$transf = transf
  }
  invisible(plot.data)
}
