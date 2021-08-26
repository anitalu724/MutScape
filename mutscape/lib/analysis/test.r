lsqnonneg <- function(C, d) {
    stopifnot(is.numeric(C), is.numeric(d))
    if (!is.matrix(C) || !is.vector(d))
        stop("Argument 'C' must be a matrix, 'd' a vector.")
    m <- nrow(C); n <- ncol(C)
    if (m != length(d))
        stop("Arguments 'C' and 'd' have nonconformable dimensions.")

    tol = 10 * eps() * norm(C, type = "2") * (max(n, m) + 1)

    x  <- rep(0, n)             # initial point
    P  <- logical(n); Z <- !P   # non-active / active columns
    
    resid <- d - C %*% x
    w <- t(C) %*% resid
    wz <- numeric(n)

    # iteration parameters
    outeriter <- 0; it <- 0
    itmax <- 3 * n; exitflag <- 1

    while (any(Z) && any(w[Z] > tol)) {
        outeriter <- outeriter + 1
        z <- numeric(n)
        wz <- rep(-Inf, n)
        wz[Z] <- w[Z]
        im <- which.max(wz)
        P[im] <- TRUE; Z[im] <- FALSE
        z[P] <- qr.solve(C[, P], d)
        
        while (any(z[P] <= 0)) {
            it <- it + 1
            if (it > itmax) stop("Iteration count exceeded.")

            Q <- (z <= 0) & P
            alpha <- min(x[Q] / (x[Q] - z[Q]))
            x <- x + alpha*(z - x)
            Z <- ((abs(x) < tol) & P) | Z
            P <- !Z
            z <- numeric(n)
            z[P] <- qr.solve(C[, P], d)
        }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
    }
    return(list(x = x, resid.norm = sum(resid*resid)))
}





if tol is None:
        tol = 10*eps*norm1(C)*(max(C.shape)+1)

    C = numpy.asarray(C)

    (m,n) = C.shape
    P = numpy.zeros(n)
    Z = numpy.arange(1, n+1)

    if x0 is None:
        x=P
    else:
        if any(x0 < 0):
            x=P
        else:
            x=x0

    ZZ=Z

    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)

    outeriter=0
    it=0
    itmax=itmax_factor*n
    exitflag=1

    # outer loop to put variables into set to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol):
        outeriter += 1

        t = w[ZZ-1].argmax()
        t = ZZ[t]

        P[t-1]=t
        Z[t-1]=0

        PP = numpy.where(P <> 0)[0]+1
        ZZ = numpy.where(Z <> 0)[0]+1

        CP = numpy.zeros(C.shape)

        CP[:, PP-1] = C[:, PP-1]
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))

        z=numpy.dot(numpy.linalg.pinv(CP), d)

        z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))

        # inner loop to remove elements from the positve set which no longer belong
        while numpy.any(z[PP-1] <= tol):
            it += 1

            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the tolerance tol. (max_error=%d)' % (it, max_error))

            QQ = numpy.where((z <= tol) & (P <> 0))[0]
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
            x = x + alpha*(z-x)

            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1
            Z[ij-1] = ij
            P[ij-1] = numpy.zeros(max(ij.shape))
            PP = numpy.where(P <> 0)[0]+1
            ZZ = numpy.where(Z <> 0)[0]+1

            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))

            z=numpy.dot(numpy.linalg.pinv(CP), d)
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))

        x = z
        resid = d - numpy.dot(C, x)
        w = numpy.dot(C.T, resid)

    return (x, sum(resid * resid), resid)