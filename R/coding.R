### Functions for handling coded data

# parse a coding formula.
# Normally expects the form coded ~ (orig - center) / divisor,
# But any linear expression in one variable is OK.
parse.coding = function(form) {
  if (!inherits(form, "formula")) stop("Need a formula")
  if (length(form) < 3) stop("Formula lacks a left-hand-side")
  nm = all.vars(form)
  names(nm) = c("coded", "orig")
  rhs = as.character(form)[3]
  a = eval(parse(text = sub(nm[2], "0", rhs)))
  b = eval(parse(text = sub(nm[2], "1", rhs)))
  d = 1 / (b - a)
  list(names = nm, const=c(center = signif(-a * d, 4), divisor = signif(d, 4)))
}


# Code the data in a data.frame; may specify as arguments or in a list
coded.data = function(data, ..., formulas=list(...)) {
  nm = names(data)
  codings = list()
  for (f in formulas) {
    info = parse.coding(f)
    cod = info$names[["coded"]]
    org = info$names[["orig"]]
    codings[[cod]] = f
    if (!is.null(data[[org]])) {
      data[[org]] = (data[[org]] - info$const[["center"]]) / info$const[["divisor"]]
      nm[nm==org] = cod
    }
  }
  names(data) = nm
  attr(data, "codings") = codings
  class(data) = c("coded.data", class(data))
  data
}

# Add coding attributes to already-coded data
as.coded.data = function(data, ..., formulas=list(...)) {
  codings = list()
  for (f in formulas) {
    info = parse.coding(f)
    cod = info$names[["coded"]]
    codings[[cod]] = f
  }
  attr(data, "codings") = codings
  class(data) = c("coded.data", class(data))
  data
}

print.coded.data = function(x, ...) {
  print.data.frame (x, ...)
  cat ("\nVariable codings ...\n")
  sapply (attr(x, "codings"), print)
  invisible (x)
}

# Transform coded data back to original scale
decode.data = function(data) {
  nm = names(data)
  codings = attr(data, "codings")
  if (!is.null(codings)) for (f in codings) {
    info = parse.coding(f)
    cod = info$names[["coded"]]
    org = info$names[["orig"]]
    if (!is.null(data[[cod]])) {
      data[[cod]] = info$const[["divisor"]] * data[[cod]] + info$const[["center"]]
      nm[nm==cod] = org
    }
  }
  names(data) = nm
  attr(data, "codings") = NULL
  cls = class(data)
  class(data) = cls[cls != "coded.data"]
  data
}

# Convert values in X to original based on codings
# Returns an object of the same type (data.frame, matrix, or vector)
# names (or column names) of X must match those used in codings
code2val = function(X, codings) {
  if (is.matrix(X))
    Z = as.matrix (code2val(as.data.frame(X), codings))
  else if (is.vector(X)) {
    nm = names(X)
    X = as.data.frame(as.list(X))
    names(X) = nm
    X = code2val (X, codings)
    Z = as.numeric (as.matrix (X))
    names(Z) = names(X)
  }
  else if (is.data.frame(X)) {
    attr(X, "codings") = codings
    Z = decode.data(X)
  }
  else stop ("Can't convert this object")
  Z
}

# Convert values in X to original based on codings
# Returns an object of the same type (data.frame, matrix, or vector)
# names (or column names) of X must match those used in codings
val2code = function(X, codings) {
  if (is.matrix(X))
    Z = as.matrix (val2code (as.data.frame(X), codings))
  else if (is.vector(X)) {
    nm = names(X)
    X = as.data.frame(as.list(X))
    names(X) = nm
    X = val2code (X, codings)
    Z = as.numeric (as.matrix (X))
    names(Z) = names(X)
  }
  else if (is.data.frame(X)) {
    Z = coded.data(X, formulas=codings)
    attr(Z, "codings") = NULL
    cls = class(Z)
    class(Z) = cls[cls != "coded.data"]
  }
  else stop ("Can't convert this object")
  Z
}

# Primitive accessor methods
codings = function(object)
    UseMethod("codings")

codings.coded.data = function(object)
    attr(object, "codings")


  