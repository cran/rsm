### Functions for handling coded data

# parse a formula of the form coded ~ (orig - center) / midrange
parse.coding = function(form) {
  nm = c(coded = as.character(form[[2]]), orig = "")
  const = c(center = 0, divisor = as.numeric(form[[3]][[3]]))
  form = form[[3]][[2]][[2]]
  nm[2] = as.character(form[[2]])
  const[1] = as.numeric(form[[3]])
  list(names = nm, const = const)
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

