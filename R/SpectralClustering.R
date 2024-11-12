#' Spectral clustering
#'
#' Implementation of the spectral clustering algorithm of Ng et al. (2001)
#'
#' @param Dt Object that contains the data.
#' @param S The similarity matrix.
#' @param Laplacian Normalization of the Laplacian matrix. Symmetric or random walk.
#' @param k Number of clusters
#' @param sigma Sigma in the Gaussian similarity measure.
#' @param sparse Indicates whether sparse matrices are used.
#' @param OptKmeans Measure to select the optimal k-means clustering solution.
#' @param LocalScaling Whether local scaling has to be applied.
#' @param kNNLocalScaling Number of nearest neighbors when using local scaling.
#' @param ... Arguments passed to kmeans.
#'
#' @references Ng, A., Jordan, M. & Weiss, Y. (2001). On spectral clustering: Analysis and an algorithm. \emph{Advances in Neural Information Processing Systems}, 14.
#'
#' @examples
#' library(kernlab)
#' data(spirals)
#' x = SpectralClusteringGridSearch(spirals, k = 2)
#' plot(x$Results)
#'
#' Dt = as.data.frame(T4cluster::genSMILEY())[, -3]
#' x  = SpectralClustering(Dt, k = 4, LocalScaling = "kernlab")
#' x
#' plot(x)
SpectralClustering <- function(Dt, S, Laplacian = c("sym", "rw"), k = 2:10, sigma = 1, sparse = F,
                               OptKmeans = c("SSw", "PseudoF", "DBindex", "DunnIndex"),
                               LocalScaling = c("no", "kernlab", "k = 7", "kNN"), kNNLocalScaling = NULL,
                               ...) {
  if((missing(Dt) & missing(S)) | (!missing(Dt) & !missing(S)))
    stop("Provide either Dt or S.")

  OptKmeans = match.arg(OptKmeans)
  s <- function(x) as(x, "sparseMatrix")
  arpack = F

  Call         = match.call()
  Laplacian    = match.arg(Laplacian)
  LocalScaling = match.arg(LocalScaling)

  if(any(k < 1))
    k = k[k > 1]
  if(missing(S)) {
    n = nrow(Dt)
    M = as.matrix(Dt)
    S =
      if(LocalScaling == "no") {
        GSM(M, sigma)
      } else {
        DistM = as.matrix(dist(Dt))
        s =
          switch(LocalScaling,
                 "kernlab" = apply(DistM, 1, function(x) median(sort(x)[1:5])),
                 "k = 7" = apply(DistM, 1, function(x) sort(x)[8]),
                 "kNN" = apply(DistM, 1, function(x) sort(x)[kNNLocalScaling + 1]))
        GSMLocalScaling(M, s)
      }
  } else {
    if(!is.matrix(S))
      S = as.matrix(S)
    if(!Matrix::isSymmetric(S))
      stop("Similarity matrix is not symmetric.")
    n = nrow(S)
  }
  if(sparse)
    S %<>% s
  I = diag(1, ncol = nrow(S), nrow = ncol(S))

  Dinv     = diag(1 / colSums(S))
  Dinvsqrt = diag(1 / sqrt(rowSums(S)))

  if(sparse) {
    I        %<>% s
    Dinvsqrt %<>% s
    Dinv     %<>% s
  }

  L = switch(Laplacian,
             sym = Dinvsqrt %*% S %*% Dinvsqrt,
             rw = I - Dinv %*% S)
  if(any(is.nan(if(inherits(L, "sparseMatrix")) L@x else L)) |
     any(is.infinite(if(inherits(L, "sparseMatrix")) L@x else L)))
    return(NULL)
  if(sparse)
    L = as(L, "dgeMatrix")
  EigenD =
    if (arpack) {
      func <- function(x, extra = NULL) as.vector(L %*% x)
      arpack(
        func,
        options = list(
          n = nrow(L),
          nev = max(k),
          ncv = ceiling(max(k) * 1.5),
          which = "LM",
          maxiter = 200
        ),
        sym = TRUE,
        complex = FALSE
      )
    } else {
      eigen(L, symmetric = T)
    }

  e = EigenD$vectors
  l = EigenD$values

  if(length(k) > 1) {
    LBartlett = diag(n) - Dinv %*% S
    k = speccalt::bartlett(eigen(LBartlett)$values, maxk = max(k))$k
  }

  ek = e[, 1:k, drop = F]
  X  = ek / sqrt(rowSums(ek^2))
  BestSeed = lapply(1:50, function(i) {
    set.seed(i)
    tryCatch(kmeans(X, k, ...), warning = function(w) w, error = function(e) e)
  })
  BestSeed = BestSeed[!sapply(BestSeed, function(x) inherits(x, "warning")) &
                        !sapply(BestSeed, function(x) inherits(x, "error"))]
  if(length(BestSeed) < 10) {
    BestSeed = lapply(50:150, function(i) {
      set.seed(i)
      tryCatch(kmeans(X, k, ...), warning = function(w) w, error = function(e) e)
    })
    BestSeed = BestSeed[!sapply(BestSeed, function(x) inherits(x, "warning")) &
                          !sapply(BestSeed, function(x) inherits(x, "error"))]
    if(length(BestSeed) < 10)
      stop("Unstable results")
  }
  SSw  = sapply(BestSeed, "[[", "tot.withinss")
  Seed = which.min(SSw)
  set.seed(Seed)
  km   = kmeans(X, k, ...)
  Results =
    list(
      Call = Call,
      k = k,
      cluster = km$cluster,
      Dt = if(!missing(Dt)) Dt else NULL,
      EigenVec = X,
      S = S,
      kmeans = km
    )
  class(Results) = "SpectralClustering"
  return(Results)
}

SpectralClusteringGridSearch <- function(Dt, Laplacian = c("sym", "rw"), k = NULL,
                                         sigma = NULL, sparse = F,
                                         Parallel = F, NrCores = parallel::detectCores() - 2,
                                         OptFn    = function(x) x$kmeans$tot.withinss,
                                         OptValue = min) {
  s <- function(x) as(x, "sparseMatrix")

  if(is.null(k))
    stop("Provide sequence of k values.")

  Call = match.call()
  Argz = as.list(Call)[-1]
  Laplacian      = match.arg(Laplacian)
  Argz$Laplacian = Laplacian

  if(any(k < 1))
    k = k[k > 1]
  n = nrow(Dt)
  M = as.matrix(Dt)
  ktmp = as.matrix(dist(M))

  if(is.null(sigma)) {
    kmax <- max(ktmp)
    kmin <- min(ktmp + diag(rep(Inf, dim(ktmp)[1])))
    kmea <- mean(ktmp)
    lsmin <- log2(kmin)
    lsmax <- log2(kmax)
    midmax <- min(c(2*kmea, kmax))
    midmin <- max(c(kmea/2,kmin))
    rtmp <- c(seq(midmin,0.9*kmea,0.05*kmea), seq(kmea,midmax,0.08*kmea))
    if ((lsmax - (Re(log2(midmax))+0.5)) < 0.5)
      step <- (lsmax - (Re(log2(midmax))+0.5))
    else
      step <- 0.5
    if (((Re(log2(midmin))-0.5)-lsmin) < 0.5 )
      stepm <-  ((Re(log2(midmin))-0.5) - lsmin)
    else
      stepm <- 0.5

    tmpsig <- c(2^(seq(lsmin, (Re(log2(midmin))-0.5), stepm)), rtmp, 2^(seq(Re(log2(midmax)) + 0.5, lsmax,step)))
    tmpsig = sort(unique(c(tmpsig, 1)))
  }
  pGrid = expand.grid(k, tmpsig)
  colnames(pGrid) = c("k", "sigma")
  pGrid = pGrid[order(pGrid$k), ]

  cat("\nPerforming grid search.\n")
  Res =
    if(Parallel){
      require(doParallel, quietly = T)
      Argz[c("Parallel", "OptFn", "OptValue")] = NULL
      Pkgs = c(Pkgs, "RcppFunctions")
      cl   = makeCluster(NrCores)
      on.exit(stopCluster(cl))
      clusterExport(cl, "Pkgs", envir = environment())
      clusterEvalQ(cl, lapply(Pkgs, library, character.only = T))
      clusterExport(cl,  c("pGrid", "Dt", "Argz", "SpectralClustering"), envir = environment())
      parLapply(cl = cl, seq_len(nrow(pGrid)), function(i) {
        Argz$k     = pGrid$k[i]
        Argz$sigma = pGrid$sigma[i]
        tryCatch(do.call("SpectralClustering", Argz), warning = function(w) NULL,
                 error = function(e) NULL)
      })
    } else{
      Argz[c("Parallel", "OptFn", "OptValue")] = NULL
      lapply(seq_len(nrow(pGrid)), function(i) {
        Argz$k     = pGrid$k[i]
        Argz$sigma = pGrid$sigma[i]
        tryCatch(do.call("SpectralClustering", Argz), warning = function(w) NULL,
                 error = function(e) NULL)
      })
    }
  pGrid  = pGrid[!sapply(Res, is.null), ]
  Res    = Res[!sapply(Res, is.null)]
  Crit   = sapply(Res, OptFn)
  OptRes = Res[[which(Crit == OptValue(Crit))]]
  FinalRes =
    list(
      Results = OptRes,
      gridSearch = pGrid,
      AllRes  = Res,
      Crit    = Crit
    )
  class(FinalRes) = "SpectralClusteringGridSearch"
  FinalRes
}

print.SpectralClusteringGridSearch <- function(x, ...) {
  print(x$Results$Call)
  cat("\nSelected number of clusters =", x$k, "\n")
  print(x$Results$kmeans)
}

plot.SpectralClusteringGridSearch <- function(x, ...) {
  plot(x$Results$Dt, col = x$Results$cluster, ...)
}

print.SpectralClustering <- function(x, ...) {
  print(x$Call)
  cat("\nSelected number of clusters =", x$k, "\n")
  print(x$kmeans)
}

plot.SpectralClustering <- function(x, ...) {
  plot(x$Dt, col = x$cluster, ...)
}

