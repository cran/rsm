### generate a CCD
# basis - formula - lhs (opt): dep var name(s);  rhs: var names for basic grid
# 

ccd = function(basis, generators, blocks="Block", n0=4, alpha="orthogonal", 
               wbreps=1, bbreps=1, randomize=TRUE, coding)
{
    if (inherits(basis, "formula"))
        xvars = all.vars(basis[[length(basis)]])
    else if (is.numeric(basis))
        xvars = paste("x", 1:basis, sep="")
    else
        stop("'basis' must be an integer or a formula")
        
    args = lapply(xvars, function(nm) c(-1,1))
    names(args) = xvars
    cube = do.call(expand.grid, args)
    
    if (!missing(generators)) {
        if (!is.list(generators)) generators = list(generators)
        for (gen in generators) {
            gen = as.character(gen)        
            cube[[gen[[2]]]] = with(cube, eval (parse (text = as.character(gen[[3]]))))
        }
    }
    
    k = ncol(cube)
    
    ### At first, star will be face-centered (as if alpha = 1); we'll scale it later
    star = as.data.frame (matrix(0, nrow = 2*k, ncol = k))
    xvars = names(star) = names (cube)
    for (j in 1:k) star[c(2*j-1,2*j), j] = c(-1, 1)
    
    if (length(wbreps) == 1) wbreps = rep(wbreps, 2)
    if (length(bbreps) == 1) bbreps = rep(bbreps, 2)
    
    # within-block reps
    if (wbreps[1] > 1) cube = cube[rep(1:nrow(cube), wbreps[1]), ]
    if (wbreps[2] > 1) star = star[rep(1:nrow(star), wbreps[2]), ]
    
    # Fractional blocking
    if (is.character(blocks)) {
        blknm = blocks
        nblev = 1
        blk = rep(1, nrow(cube))
        chkterm = ""
    }
    else if (inherits(blocks, "formula")) {
        blknm = as.character(blocks[[2]])
        what = as.character(blocks[[3]][[1]])
        if (what == "*") 
            gens = as.character(blocks[3])
        else 
            gens = as.character(blocks[[3]])[-1]
        bgen = lapply(gens, function(g) with(cube, eval(parse(text=g))))
        blk = as.numeric(factor(do.call(paste, bgen)))
        nblev = max(blk)
        chkterm = "factor(blk) + "
    }
    else
        stop("'blocks' must be a string or a formula")
        
    # Check for aliasing
    v = paste(names(cube), collapse=",")
    fake.resp = rnorm(nrow(cube))
    fstg = paste("fake.resp ~", chkterm, "FO(", v, ") + TWI(", v, ")")
    modl = lm(formula(fstg), data=cube)
    if (any(is.na(coef(modl))))
        warning("Some 1st or 2nd-order terms are aliased in the cube portion of this design")


    # center points
    zero = as.data.frame(matrix(rep(0, k), nrow=1))
    names(zero) = names(cube)
    if (length(n0) == 1) n0 = c(n0, n0)
    if (n0[1] > 0) {
        cube = rbind(cube, zero[rep(1, nblev*n0[1]), ])
        blk = c(blk, rep(unique(blk), n0[1]))
    }
    if (n0[2] > 0) 
        star = rbind(star, zero[rep(1, n0[2]), ])
        
    # Block reps
    nc = nrow(cube)
    if (bbreps[1] > 1) {
        cube = cube[rep(1:nc, bbreps[1]), ]
        blk = nblev*rep(0:(bbreps[1]-1), rep(nc, bbreps[1])) 
              + rep(blk, bbreps[1])
        nblev = max(blk)
    }
    ns = nrow(star)
    if (bbreps[2] > 1) star = star[rep(1:ns, bbreps[2]), ]
    sblk = rep( (1+nblev):(bbreps[2]+nblev), rep(ns, bbreps[2]))
    
    # Figure out alpha if given as criterion
    if (is.character(alpha)) {
        c.ii = sum(cube[[1]]^2)
        s.ii = sum(star[[1]]^2)
        what = pmatch(alpha, c("rotatable", "orthogonal"))
        if (is.na(what)) stop ("alpha must be 'rotatable', 'orthogonal', or a value")
        if (what==1)
            alpha = (2 * c.ii / s.ii) ^ .25
        else
            alpha = sqrt(nrow(star) / s.ii  * c.ii / nrow(cube))
    }
    star = star * alpha
    
    # append blocking variable
    cube = cbind(blk, cube)
    star = cbind(sblk, star)
    names(cube)[1] = names(star)[1] = blknm
    
    # Figure out row names
    ord = order(blk, 1:nrow(cube))
    cube = cube[ord, ]
    blk = blk[ord]
    row.names(cube) = paste("C", blk, ".", rep(1:(nrow(cube)/nblev), nblev), sep="")
    row.names(star) = paste("S", sblk, ".", rep(1:(nrow(star)/bbreps[2]), bbreps[2]), sep="")

    # assemble design
    des = rbind(cube, star)
    
    # Add vars from left-hand side, if any
    if (inherits(basis, "formula") & (length(basis) > 2)) {
        yvars = all.vars(basis[[2]])
        for (v in yvars)  des[[v]] = NA
    }
    
    # Figure out sort order
    if (randomize) {
        ord = order(des[[1]] + runif(nrow(des)))
        des = des[ord, ]
    }
    
    des[[1]] = factor(des[[1]])
    
    if (!missing(coding))
        des = as.coded.data (des, formulas=coding)
    
    des
}

