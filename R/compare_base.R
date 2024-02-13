########################## BASE FUNCTIONS FOR THE SIMULATIONS ##################
#source('sbm_functions.R')

if (!require(irlba)) {
  install.packages('irlba', dependencies = T)
  require(irlba)
}

irlbaMod <-
  function (A, nu=5, nv=5, adjust=3, aug=c("ritz","harm"), sigma=c("ls","ss"), 
            maxit=1000, m_b=20, reorth=2, tol=1e-6, V=NULL,
            matmul=NULL)
  {
    # ---------------------------------------------------------------------
    # Check input parameters
    # ---------------------------------------------------------------------
    eps <- .Machine$double.eps
    # Profiling option
    options(digits.secs=3)
    m <- nrow(A[[1]])
    n <- ncol(A[[1]])
    k <- max(nu,nv)
    interchange <- FALSE
    sigma = match.arg(sigma)
    aug   = match.arg(aug)
    # Interchange dimensions m,n so that dim(A'A) = min(m,n) when seeking
    # the smallest singular values. This avoids finding zero smallest 
    # singular values.
    if (n>m && sigma=="ss") {
      t <- m
      m <- n
      n <- t
      interchange <- TRUE
    }
    
    # Increase the number of desired signular values by 'adjust' to
    # help convergence. k is re-adjusted as vectors converge--this is
    # only an initial value;
    k_org <- k;
    k <- k + adjust;
    if (k<=0)  stop ("k must be positive")
    if (k>min(m,n)) stop ("k must be less than min(m,n)+adjust")
    if (m_b<=1) stop ("m_b must be greater than 1")
    if (tol<0) stop ("tol must be non-negative")
    if (maxit<=0) stop ("maxit must be positive")
    if (m_b>= min(n,m)) {
      m_b <- floor(min(n,m)-0.1)
      if (m_b-k-1<0) {
        adjust <- 0
        k <- m_b-1
      }
    }
    if (m_b-k-1<0) m_b <- ceiling(k+1+0.1)
    if (m_b>=min(m,n)) {
      m_b <- floor(min(m,n)-0.1)
      adjust <- 0
      k <- m_b - 1
    }
    if (tol<eps) tol <- eps
    
    # Allocate memory for W and F:
    W <- matrix(0.0,m,m_b) 
    F <- matrix(0.0,n,1)
    # If starting matrix V is not given then set V to be an
    # (n x 1) matrix of normally distributed random numbers.
    # In any case, allocate V appropriate to problem size:
    if (is.null(V)) {
      V <- matrix(0.0,n,m_b)
      V[,1] <- rnorm(n)
    }
    else {
      V <- cbind(V, matrix(0.0,n,m_b-ncol(V)))
    }
    
    
    # ---------------------------------------------------------------------
    # Initialize local variables
    # ---------------------------------------------------------------------
    
    B <- NULL                  # Bidiagonal matrix
    Bsz <- NULL                # Size of B
    eps23 <- eps^(2/3)         # Used for Smax/avoids using zero
    I <- NULL                  # Indexing
    J <- NULL                  # Indexing
    iter <- 1                  # Man loop iteration count
    mprod <- 0                 # Number of matrix-vector products
    R_F <- NULL                # 2-norm of residual vector F
    sqrteps <- sqrt(eps)       #
    Smax <- 1                  # Max value of all computed singular values of
    # B est. ||A||_2
    Smin <- NULL               # Min value of all computed singular values of
    # B est. cond(A)
    SVTol <- max(sqrteps,tol)  # Tolerance for singular vector convergence
    S_B <- NULL                # Singular values of B
    U_B <- NULL                # Left singular vectors of B
    V_B <- NULL                # Right singular vectors of B
    V_B_last <- NULL           # last row of modified V_B 
    S_B2 <- NULL               # S.V. of [B ||F||]
    U_B2 <- NULL               # 
    V_B2 <- NULL               #  
    
    # ---------------------------------------------------------------------
    # Basic functions
    # ---------------------------------------------------------------------
    
    # Euclidean norm
    norm2 <- function (x) return(as.numeric(sqrt(crossprod(x))))
    
    # Orthogonalize vectors Y against vectors X. Y and X must be R matrix
    # objects (they must have a dim attribute).
    # Note: this function unnecessarily copies the contents of Y
    orthog <- function (Y,X)
    {
      if (dim(X)[2] < dim(Y)[2]) dotY <- crossprod (X,Y)
      else dotY <- t (crossprod(Y,X))
      return (Y - X %*% dotY)
    }
    
    # Convergence tests
    # Input parameters
    # Bsz            Number of rows of the bidiagonal matrix B
    # tol
    # k_org
    # U_B            Left singular vectors of small matrix B
    # S_B            Singular values of B
    # V_B            Right singular vectors of B
    # residuals      
    # k
    # SVTol
    # Smax
    #
    # Output parameter list
    # converged      TRUE/FALSE
    # U_B            Left singular vectors of small matrix B
    # S_B            Singular values of B
    # V_B            Right singular vectors of B
    # k              Number of singular vectors returned 
    convtests <- function (Bsz, tol, k_org, U_B, S_B, V_B, 
                           residuals, k, SVTol, Smax)
    {
      Len_res <- sum(residuals[1:k_org] < tol*Smax)
      if (Len_res == k_org) {
        return (list(converged=TRUE, U_B=U_B[,1:k_org, drop=FALSE], 
                     S_B=S_B[1:k_org, drop=FALSE], V_B=V_B[,1:k_org, drop=FALSE], k=k) )
      } 
      #   Not converged yet...
      #   Adjust k to include more vectors as the number of vectors converge.
      Len_res <- sum(residuals[1:k_org] < SVTol*Smax)
      k <- max(k, k_org + Len_res)
      if (k > Bsz -3) k <- Bsz -3
      return (list(converged=FALSE, U_B=U_B, S_B=S_B, V_B=V_B, k=k) )
    }
    
    # ---------------------------------------------------------------------
    # Main iteration
    # ---------------------------------------------------------------------
    
    while (iter <= maxit) {
      
      # ---------------------------------------------------------------------
      # Lanczos bidiagonalization iteration
      # Compute the Lanczos bidiagonal decomposition:
      # AV  = WB
      # t(A)W = VB + Ft(E)
      # with full reorthogonalization.
      # This routine updates W,V,F,B,mprod
      # ---------------------------------------------------------------------
      j <- 1
      #   Normalize starting vector:
      if (iter==1) V[,1] <- V[,1, drop=FALSE]/norm2(V[,1, drop=FALSE]) 
      else j <- k + 1
      
      #   Compute W=AV (the use of as.matrix here converts Matrix class objects)
      if(!is.null(matmul)) {
        #     User-specified matrix multiply function
        #print(matmul(A, V[,j,drop=FALSE], transpose=TRUE))
        if(interchange)
          W[,j] <- matmul(A, V[,j,drop=FALSE], transpose=TRUE)
        else
          W[,j] <- matmul(A, V[,j,drop=FALSE])
      }
      else {
        if (interchange)  W[,j] <- t (as.matrix(crossprod (V[,j,drop=FALSE], A)))
        else              W[,j] <- as.matrix(A %*% V[,j, drop=FALSE])
      }
      mprod <- mprod + 1
      
      #   Orthogonalize
      if (iter != 1) {
        W[,j] <- orthog (W[,j, drop=FALSE], W[,1:(j-1), drop=FALSE])
      }
      
      S <- norm2(W[,j, drop=FALSE])
      #   Check for linearly dependent vectors
      if ((S < SVTol) && (j==1)) stop ("Starting vector near the null space")
      if (S < SVTol) {
        W[,j] <- rnorm(nrow(W))
        W[,j] <- orthog(W[,j, drop=FALSE],W[,1:(j-1), drop=FALSE])
        W[,j] <- W[,j, drop=FALSE]/norm2(W[,j, drop=FALSE])
        S <- 0 
      }
      else W[,j] <- W[,j, drop=FALSE]/S
      
      #   Lanczos process
      while (j <= m_b) {
        if(!is.null(matmul)) {
          #       User-specified matrix multiply function
          if(interchange)
            F <- matmul(A, W[,j,drop=FALSE])
          else
            F <- matmul(A, W[,j,drop=FALSE], transpose=TRUE)
        }
        else{
          if (interchange) F <- as.matrix(A %*% W[,j, drop=FALSE])
          else F <- t(as.matrix(crossprod(W[,j,drop=FALSE],A)))
        }
        
        mprod <- mprod + 1
        F <- F - S*V[,j, drop=FALSE]
        #     Orthogonalize
        F <- orthog(F,V[,1:j, drop=FALSE])
        
        if (j+1 <= m_b) {
          R <- norm2(F)
          #       Check for linear dependence
          if (R<=SVTol) {
            F <- matrix(rnorm(dim(V)[1]),dim(V)[1],1)
            F <- orthog(F, V[,1:j, drop=FALSE])
            V[,j+1] <- F/norm2(F)
            R <- 0 
          }
          else V[,j+1] <- F/R
          
          #       Compute block diagonal matrix 
          if (is.null(B)) B <- cbind(S, R)
          else            B <- rbind(cbind(B,0),c(rep(0,j-1),S,R))
          
          if(!is.null(matmul)) {
            #         User-specified matrix multiply function
            if(interchange)
              W[,j+1] <- matmul(A, V[,j+1,drop=FALSE], transpose=TRUE)
            else
              W[,j+1] <- matmul(A, V[,j+1,drop=FALSE])
          }
          else{
            if (interchange) 
              W[,j+1] <- t (as.matrix(crossprod (V[,j+1, drop=FALSE],A)))
            else             W[,j+1] <- as.matrix(A %*% V[,j+1, drop=FALSE])
          }
          mprod <- mprod + 1
          
          #       One step of the classical Gram-Schmidt process
          W[,j+1] <- W[,j+1, drop=FALSE] - W[,j, drop=FALSE]*R
          
          #       Full reorthogonalization of W
          if (iter==1 || reorth==2)
            W[,j+1] <- orthog(W[,j+1, drop=FALSE],W[,1:j, drop=FALSE])
          S <- norm2(W[,j+1, drop=FALSE])
          if (S<=SVTol) {
            W[,j+1] <- rnorm(nrow(W))
            W[,j+1] <- orthog(W[,j+1, drop=FALSE],W[,1:j, drop=FALSE])
            W[,j+1] <- W[,j+1, drop=FALSE]/norm2(W[,j+1, drop=FALSE])
            S <- 0
          }
          else W[,j+1] <- W[,j+1, drop=FALSE]/S
        }
        else {
          #       Add a last block to matrix B
          B <- rbind(B,c(rep(0,j-1),S))
        }
        j <- j + 1
      }
      #cat ("iter = ",iter," j = ",j-1, "mprod = ",mprod,"\n",file=stderr())
      # ---------------------------------------------------------------------
      # (End of the Lanczos bidiagonalization part)
      # ---------------------------------------------------------------------
      
      Bsz <- nrow(B)
      R_F <- norm2(F)
      F <- F/R_F
      #   Compute singular triplets of B. Expect svd to return s.v.s in order
      #   from largest to smallest.
      Bsvd <- svd(B)
      
      #   Estimate ||A|| using the largest singular value over all iterations
      #   and estimate the cond(A) using approximations to the largest and 
      #   smallest singular values. If a small singular value is less than sqrteps
      #   use only Ritz vectors to augment and require two-sided reorthogonalization.
      if (iter ==1) {
        Smax <- Bsvd$d[1]
        Smin <- Bsvd$d[Bsz]
      }
      else {
        Smax <- max(Smax, Bsvd$d[1])
        Smin <- min(Smin, Bsvd$d[Bsz])
      }
      Smax <- max(eps23,Smax)
      if ((Smin/Smax < sqrteps) && reorth <2) {
        warning ("The matrix is ill-conditioned. Each basis will be reorthogonalized.")
        reorth <- 2
        aug <- "ritz"
      }
      
      #   Re-order the singular values accordingly.
      if (sigma=="ss") {
        jj <- seq (ncol (Bsvd$u), 1, by=-1)
        Bsvd$u <- Bsvd$u[,jj]
        Bsvd$d <- Bsvd$d[jj]
        Bsvd$v <- Bsvd$v[,jj]
      }
      
      #   Compute the residuals
      R <- R_F * Bsvd$u[Bsz,, drop=FALSE]
      
      #   Check for convergence
      ct <- convtests(Bsz, tol, k_org, Bsvd$u, Bsvd$d, Bsvd$v, abs(R), k, SVTol, Smax)
      k <- ct$k
      
      #   If all desired singular values converged, then exit main loop
      if (ct$converged) break
      
      if (iter>=maxit) break
      
      #   Compute the starting vectors and first block of B[1:k,1:(k+1), drop=FALSE]
      if (aug=="harm") {
        #     Update the SVD of B to be the svd of [B ||F||E_m]
        Bsvd2.d <- Bsvd$d
        Bsvd2.d <- diag(Bsvd2.d, nrow=length(Bsvd2.d))
        Bsvd2 <- svd (cbind (Bsvd2.d, t (R)))
        if (sigma=="ss") {
          jj <- seq (ncol (Bsvd2$u), 1, by=-1)
          Bsvd2$u <- Bsvd2$u[,jj]
          Bsvd2$d <- Bsvd2$d[jj]
          Bsvd2$v <- Bsvd2$v[,jj]
        }
        Bsvd$d <- Bsvd2$d
        Bsvd$u <- Bsvd$u %*% Bsvd2$u
        Bsvd$v <- cbind (rbind (Bsvd$v, rep (0,Bsz)), 
                         c (rep (0,Bsz), 1)) %*% Bsvd2$v
        V_B_last <- Bsvd$v [Bsz + 1, 1:k, drop=FALSE]
        s <- R_F * solve (B, cbind (c (rep(0,Bsz-1), 1)))
        Bsvd$v <- Bsvd$v[1:Bsz, , drop=FALSE] + s %*% Bsvd$v[Bsz+1, ,drop=FALSE]
        
        qrv <- qr (cbind ( rbind (Bsvd$v[,1:k], 0), rbind (-s, 1)))
        Bsvd$v <- qr.Q(qrv)
        R <- qr.R(qrv)
        V[,1:(k+1)] <- cbind(V, F) %*% Bsvd$v
        
        #     Update and compute the k x k+1 part of B
        UT <- t(R[1:(k+1), 1:k, drop=FALSE] + R[,k+1,drop=FALSE] %*% V_B_last)
        B <- diag(Bsvd$d[1:k],nrow=k) %*% (UT*upper.tri(UT,diag=TRUE))
      }
      else {
        #     Use the Ritz vectors
        V[,1:(k + dim(F)[2])] <- cbind(V[,1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[,1:k, drop=FALSE], F)
        B <- cbind( diag(Bsvd$d[1:k],nrow=k), R[1:k, drop=FALSE])
      }
      
      #   Update the left approximate singular vectors
      W[,1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[,1:k, drop=FALSE]
      
      iter <- iter + 1
    }
    # ---------------------------------------------------------------------
    # End of the main iteration loop
    # Output results
    # ---------------------------------------------------------------------
    d <- Bsvd$d[1:k_org]
    u <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[,1:k_org, drop=FALSE]
    v <- V[,1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[,1:k_org, drop=FALSE]
    # Adjust ordering if smallest singular values selected so that singular values
    # are reported in non-increasing order.
    if(sigma=="ss") {
      reverse <- seq(length(d),1)
      d <- d[reverse]
      u <- u[,reverse,drop=FALSE]
      v <- v[,reverse,drop=FALSE]
    }
    return (list(d=d, u=u[,1:nu,drop=FALSE], v=v[,1:nv,drop=FALSE], iter=iter,mprod=mprod))
  }


getCcaClusters = function(adjacencyMat, covariates, nBlocks,
                          method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  gilSingVec = getGilSvd(getGraphMatrix(adjacencyMat, method), covariates,
                         nBlocks)$singVec
  
  return( kmeans(gilSingVec, nBlocks, nstart = randStarts)$cluster )
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering
# ---------------------------------------------------------------------
getCascClusters = function(adjacencyMat, covariates, hTuningParam,
                           nBlocks, method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  cascSingVec = getCascSvd(getGraphMatrix(adjacencyMat, method), covariates,
                           hTuningParam, nBlocks)$singVec
  
  return( kmeans(cascSingVec, nBlocks, nstart = randStarts)$cluster )
  
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering takes graphMat
# ---------------------------------------------------------------------
getCascResults = function(graphMat, covariates, hTuningParam,
                          nBlocks) {
  
  randStarts = 10 #number of random starts for kmeans
  
  cascSvd = getCascSvd(graphMat, covariates, hTuningParam, nBlocks)
  
  ortho = getOrtho(graphMat, covariates, cascSvd$singVec, cascSvd$singVal,
                   hTuningParam, nBlocks)
  
  kmeansResults = kmeans(cascSvd$singVec, nBlocks, nstart = randStarts)
  
  return( list(cluster = kmeansResults$cluster,
               wcss = kmeansResults$tot.withinss,
               singGap = cascSvd$singVal[nBlocks] -
                 cascSvd$singVal[nBlocks + 1],
               orthoL = ortho$orthoL,
               orthoX = ortho$orthoX,
               singVecK = cascSvd$singVec[, nBlocks],
               singVecKPlus = cascSvd$singVecKPlus) )    
}

# ---------------------------------------------------------------------
# returns cluster memberships for SC based graph clustering
# ---------------------------------------------------------------------
getGraphScClusters = function(adjacencyMat, nBlocks,
                              method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  scSingVec = getGraphScSvd(getGraphMatrix(adjacencyMat, method),
                            nBlocks)$singVec
  
  return( kmeans(scSingVec, nBlocks, nstart = randStarts)$cluster )
}

# ---------------------------------------------------------------------
# returns cluster memberships for SC based covariate clustering
# ---------------------------------------------------------------------
getCovScClusters = function(covariates, nBlocks) {
  
  randStarts = 10 #number of random starts for kmeans
  
  scSingVec = getCovScSvd(covariates, nBlocks)$singVec
  
  return( kmeans(scSingVec, nBlocks, nstart = randStarts)$cluster )
}


# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns left singular vectors and values for GIL based clustering
# ---------------------------------------------------------------------
getGilSvd = function(graphMat, covariates, nBlocks) {
  
  #insure irlba internal representation is large enough
  if(nBlocks > 10) {
    internalDim = 2 * nBlocks
  }
  else {
    internalDim = 20
  }
  
  #approximate the generalized inverse of L using top eigenvectors
  svdL = irlba(graphMat, nu = 2*nBlocks, m_b = internalDim)
  graphMatP = svdL$u %*% Diagonal(1/svdL$d) %*% svdL$u^T
  
  #define a custom matrix vector multiply function
  matrixMulti = function(aList, aVector, transposeBool) {
    return( as.vector(aList$graphMat %*%
                        (aList$covariates %*% (aList$covariates^T
                                               %*% aVector))) )
  } 
  
  singDecomp = irlbaMod(list(graphMat = graphMatP, covariates = covariates,
                             hTuningParam = hTuningParam), nu = nBlocks + 1, nv = 0,
                        m_b = internalDim, matmul = matrixMulti) 
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CCA based clustering
# ---------------------------------------------------------------------
getCcaSvd = function(graphMat, covariates, nBlocks) {
  
  singDecomp = svd(graphMat %*% covariates, nu = nBlocks, nv = 0)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CASC based clustering
# ---------------------------------------------------------------------
getCascSvd = function(graphMat, covariates, hTuningParam, nBlocks) {
  
  #insure irlba internal representation is large enough
  internalDim = max(2*nBlocks, 20)
  
  #define a custom matrix vector multiply function
  matrixMulti = function(aList, aVector, transposeBool) {
    return( as.vector(aList$graphMat %*% aVector +
                        aList$hTuningParam * aList$covariates %*%
                        crossprod(aList$covariates, aVector) ))
  } 
  
  singDecomp = irlbaMod(list(graphMat = graphMat, covariates = covariates,
                             hTuningParam = hTuningParam), nu = nBlocks + 1, nv = 0,
                        m_b = internalDim, matmul = matrixMulti)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d,
               singVecKPlus = singDecomp$u[, nBlocks+1]) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for graph SC based clustering
# ---------------------------------------------------------------------
getGraphScSvd = function(graphMat, nBlocks) {
  
  #insure irlba internal representation is large enough
  if(nBlocks > 10) {
    internalDim = 2 * nBlocks
  }
  else {
    internalDim = 20
  }
  
  singDecomp = irlba(graphMat, nu = nBlocks, nv = 0, m_b = internalDim)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for covariate SC based clustering
# ---------------------------------------------------------------------
getCovScSvd = function(covMat, nBlocks) {
  
  singDecomp = svd(covMat, nu = nBlocks, nv = 0)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns the proportion of the eigenvalues due to X in the top eigenspace
# ---------------------------------------------------------------------
getOrtho <- function(graphMat, covariates, cascSvdSingVec, cascSvdSingVal,
                     h, nBlocks) {
  orthoL <- as.numeric((t(cascSvdSingVec[, nBlocks])%*%graphMat%*%
                          cascSvdSingVec[, nBlocks])/cascSvdSingVal[nBlocks])
  orthoX <- as.numeric(h*(t(cascSvdSingVec[, nBlocks])%*%covariates%*%
                            t(covariates)%*%cascSvdSingVec[, nBlocks])/
                         cascSvdSingVal[nBlocks])
  return( list(orthoL = orthoL/(orthoL + orthoX),
               orthoX = orthoX/(orthoL + orthoX)) )
}

# ---------------------------------------------------------------------
# returns the graph matrix corresponding to the given method
# ---------------------------------------------------------------------
getGraphMatrix = function(adjacencyMat, method) {
  
  if(method == "regLaplacian") {
    rSums = rowSums(adjacencyMat)
    tau = mean(rSums)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
    return(forceSymmetric(normMat %*% adjacencyMat %*% normMat))
  }
  else if(method == "laplacian") {
    rSums = rowSums(adjacencyMat)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums))
    return(forceSymmetric(normMat %*% adjacencyMat %*% normMat))
  }
  else if(method == "adjacency"){
    return(adjacencyMat)
  }
  else {
    stop("Method given not valid.")
  }
  
  return(-1)
}

getCascAutoClusters = function(adjacency, covariates, nBlocks,
                               method = "regLaplacian", nPoints = 100) {
  
  graphMat = getGraphMatrix(adjacency)
  rangehTuning = getTuningRange(graphMat, covariates, nBlocks)
  
  hTuningSeq = seq(rangehTuning[1], rangehTuning[2], nPoints)
  wcssVec = rep(0, nPoints)
  clusterMat = matrix(rep(0, nPoints*dim(graphMat)[1]), nrow = nPoints)
  
  for(i in 1:nPoints) {
    cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
                                 nBlocks, enhancedTuning)
    wcssVec[i] = cascResults$wcss
    clusterMat[i, ] = cascResults$cluster
  }
  
  minWcssIndex = match(min(wcssVec), wcssVec)
  
  return(clusterMat[minWcssIndex, ])
}

# ---------------------------------------------------------------------
# returns CASC optimal h tuning parameter SVD
# ---------------------------------------------------------------------
getCascAutoSvd = function(graphMat, covariates, nBlocks,
                          nPoints = 100, enhancedTuning = T) {
  
  # value for detecting a transition
  epsilon = .05
  
  rangehTuning = getTuningRange(graphMat, covariates, nBlocks)
  
  hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
                   length.out = nPoints)
  wcssVec = vector(length = nPoints)
  gapVec = vector(length = nPoints)
  orthoX = vector(length = nPoints)
  orthoL = vector(length = nPoints)
  
  for(i in 1:nPoints) {
    cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
                                 nBlocks)
    orthoX[i] = cascResults$orthoX
    orthoL[i] = cascResults$orthoL
    wcssVec[i] = cascResults$wcss
    gapVec[i] = cascResults$singGap
  }
  
  # get transition points of static eigenvectors
  subspaces = getSubspaces(orthoX, orthoL, nPoints, epsilon)
  nSubspaces = length(subspaces$subintervalStart)    
  
  if((enhancedTuning == T) & (nSubspaces > 1)) {
    
    subMinIndex = vector(length = nSubspaces)
    subMaxIndex = vector(length = nSubspaces)
    for(i in 1:nSubspaces) {
      subMinIndex[i] = which.min(wcssVec[
        subspaces$subintervalStart[i]:
          subspaces$subintervalEnd[i]]) +
        subspaces$subintervalStart[i] - 1
      subMaxIndex[i] = which.max(wcssVec[
        subspaces$subintervalStart[i]:
          subspaces$subintervalEnd[i]]) +
        subspaces$subintervalStart[i] - 1
    }
    
    # keep only those intervals that are not dominated in terms of wcss
    includeVec = (rowSums(outer(wcssVec[subMinIndex], wcssVec[subMaxIndex],
                                function(x, y) {x > y})) == 0)
    
    minCountSubspaces = ((1:nSubspaces)[includeVec == 1])[
      which.min(subspaces$orthoCounts[includeVec == 1])]
    
    # min WCSS on most overlapping set of subspaces
    startIndex = subspaces$subintervalStart[minCountSubspaces]
    endIndex = subspaces$subintervalEnd[minCountSubspaces]
    minInterval = unlist(apply(cbind(startIndex, endIndex), 1, function(x)
    {x[1]:x[2]}))
    minWcssSubindex = which.min(wcssVec[minInterval])
    hOpt = (hTuningSeq[minInterval])[minWcssSubindex]
  } else {
    hOpt = hTuningSeq[which.min(wcssVec)]
  }
  
  return( getCascSvd(graphMat, covariates, hOpt, nBlocks) )
}

# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------
getTuningRange = function(graphMatrix, covariates, nBlocks) {
  
  #insure irlba internal representation is large enough
  if(nBlocks > 10) {
    internalDim = 2 * nBlocks
  } else {
    internalDim = 20
  }
  
  singValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, m_b =
                         internalDim)$d
  singValCov = svd(covariates, nu = nBlocks)$d
  
  R = length(singValCov)
  if (R <= nBlocks) {
    denum = singValCov[R]^2
  } else {
    # denum = singValCov[nBlocks]^2
    denum = singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2
  }
  # hmax = singValGraph[1]/singValCov[nBlocks]^2
  hmax = singValGraph[1]/denum
  
  hmin = (singValGraph[nBlocks] - singValGraph[nBlocks + 1])/singValCov[1]^2
  
  return( list( hmax = hmax, hmin = hmin ) )
}

# ---------------------------------------------------------------------
# Finds leading subspace discontinuities.
# Returns the start and end of a continuous interval and
# the number of orthogonal components in the leading subspace
# on the interval.
# ---------------------------------------------------------------------
getSubspaces = function(orthoX, orthoL, nPoints, epsilon) {
  
  indicatorOut = vector(length = nPoints)
  indicatorIn = vector(length = nPoints)
  
  for(i in 1:(nPoints - 1)) {
    if((orthoX[i] < epsilon) & (orthoX[i+1] > epsilon)) {
      indicatorOut[i+1] = 1
    }
    else if((orthoL[i+1] < epsilon) & (orthoL[i] > epsilon)) {
      indicatorIn[i+1] = 1
    }
  }
  
  orthoCounts = cumsum(indicatorIn) - cumsum(indicatorOut) +
    max(cumsum(indicatorOut))
  subintervalStart = unique(c(which(indicatorIn == 1),
                              which(indicatorOut == 1)))
  subintervalEnd = sort(c(subintervalStart-1, nPoints))
  subintervalStart = sort(c(1, subintervalStart))
  orthoCounts = orthoCounts[subintervalStart]
  
  return( list(orthoCounts = orthoCounts,
               subintervalStart = subintervalStart,
               subintervalEnd = subintervalEnd) )
}



########################## OUR MODEL ###########################################
nsim <- 1000; burnin <- 500; thin <- 2; nchain <- 3; K <- 2
#' Implementation of CALF-SBM 
#' 
#' Using MCMC parameters and number of clusters as input, return raw MCMC output
#' @param links List of elements of the network, 
#' requires adjacency matrix A, matrix of covariates X, and distance matrix dis
#' @param nsim Total number of MCMC iterations per chain
#' @param burnin Number of iterations in each chain to be discarded
#' @param thin Post-burnin thinning parameter
#' @param nchain Number of MCMC chains to run
#' @param K Number of clusters
#' @param directed Boolean indicating whether to use directed network 
#' (default = FALSE)
#' @param offset Boolean indicating whether to use offset terms
#' @param beta_scale Prior standard deviation of beta terms
#' @return List of beta, z, and history of K
#' @export
calf_sbm_nimble <- function(links, nsim, burnin, thin, nchain, K, directed, 
                            offset = TRUE, beta_scale = 10){
  ## Inits
  const <- list(n = nrow(links$A), K = K)
  data <- list(A = links$A, x = links$dis)
  
  inits <- list(beta0 = rnorm(1, 0, 5)
                , beta = rnorm(const$K^2, 0, 5)
                , z = cluster::pam(links$X, const$K)$clustering
                , gamma = matrix(1, const$n, const$K)
  )
  if (offset){
    if (!directed){
      inits$theta <- log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
    } else {
      inits$theta_in <- log(rowSums(links$A) * const$n / sum(links$A) + 0.0001)
      inits$theta_out <- log(colSums(links$A) * const$n / sum(links$A) + 0.0001)
    }
  }
  ## Initialize betas
  group <- gen_factor(inits$z, links$A, links$dis)
  initial_beta <- update_beta(const$K, group$cluster)
  inits$beta0 <- initial_beta$beta0
  inits$beta <- c(initial_beta$beta)
  
  ## NIMBLE code
  monitors <- c('z', 'beta', 'beta0')
  if(offset){monitors <- c(monitors, 'sigma', 'theta')}
  
  code <- nimble::nimbleCode({
    ## Priors for parameter matrix
    beta0 ~ dnorm(mean = 0, sd = 10)
    for (a in 1:K^2){
      beta[a] ~ dnorm(mean = 0, sd = 10)
    }
    ## Priors for offset    
    if (offset){
      if (!directed) {
        for (i in 1:n){
          theta[i] ~ dnorm(mean = 0, var = sigma)
        }
        sigma ~ nimble::dinvgamma(1, 1)
      } else {
        for (i in 1:n){
          theta_in[i] ~ dnorm(mean = 0, var = sigma_in)
          theta_out[i] ~ dnorm(mean = 0, var = sigma_out)
        }
        sigma_in ~ nimble::dinvgamma(1, 1)
        sigma_out ~ nimble::dinvgamma(1, 1)
      }
    }
    ## Node membership
    for (i in 1:n){
      z[i] ~ nimble::dcat(alpha[i, 1:K])
      alpha[i, 1:K] ~ nimble::ddirch(gamma[i, 1:K])
    }
    ## Adjacency matrix from fitted values
    for (i in 1:n){
      ## Undirected network
      if (!directed){
        for (j in (i+1):n){
          if (offset){
            A[i, j] ~ dbin(expit(beta0 + 
                                 theta[i] + theta[j] + 
                                 beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
          } else {
            A[i, j] ~ dbin(expit(beta0 + 
                                 beta[(max(z[i], z[j]) - 1) * K + min(z[i], z[j])] * x[i, j]), 1)
          }
        }
      } else {
        for (j in 1:n){
          A[i, j] ~ dbin(expit(beta0 + 
                                 theta_in[i] + theta_out[j] + 
                                 beta[(z[i] - 1) * K + z[j]] * x[i, j]), 1)
        }
      }
    }
  })
  ## Compile model
  model <- nimble::nimbleModel(code, const, data, inits, check = FALSE)
  cmodel <- nimble::compileNimble(model)
  ## Compile MCMC sampler
  modelConf <- nimble::configureMCMC(model, monitors = monitors, enableWAIC = TRUE)
  modelMCMC <- nimble::buildMCMC(modelConf)
  cmodelMCMC <- nimble::compileNimble(modelMCMC, project = model) #1 min
  ## Multiple chains runner
  mcmcSamples <- nimble::runMCMC(cmodelMCMC, niter = nsim, nburnin = burnin, 
                                 thin = thin, nchains = nchain)
  mcmcSamples <- rbind(mcmcSamples$chain1, mcmcSamples$chain2, mcmcSamples$chain3)
  #print(head(mcmcSamples))
  ## Post-process samples using label.switching library
  mcmcSamples <- post_label_mcmc_samples(mcmcSamples, const$K, const$n, directed)
  return(mcmcSamples)
}

