
getPrimes = function(n) {
  
  ###################################################################################################
  # return list of prime numbers
  # e.g. getPrimes(100) = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
  ###################################################################################################
  
  # if n is not numeric, return error
  if(is.numeric(n) == FALSE) {
    stop("the input must be numeric.")
  }
  
  # if n is not an integer, transform it to the integer
  nn = as.integer(n)
  
  # if nn is less than 4, directly return list of primes. Else, implement Sieve of Eratosthenes
  if(nn < 4) {
    
    if(nn < 2) {
      return(NULL)
    } else if(nn == 2) {
      return(2)
    } else if(nn == 3) {
      return(c(2,3))
    }
    
  } else {
    
    # initialize flags
    flags = rep(TRUE, nn)
    
    # 1 is not a prime
    flags[1] = FALSE
    
    # implement Sieve of Eratosthenes
    for(k in 2:as.integer(sqrt(nn))) {
      if(flags[k] == TRUE) {
        flags[k * 2:as.integer(nn/k)] = FALSE
      }
    }
    
    # return primes vector
    return((1:nn)[flags])
    
  }
}

reldist = function(x, y) {
  
  ######################
  # relative distance
  ######################
  
  abs(x-y) / min(abs(x), abs(y))
}

absdist = function(x, y) {
  
  ######################
  # absolute distance
  ######################
  
  abs(x-y)
}

computeGradient = function(coefs, degrees, monomialPrimes) {
  
  ##################################
  # compute gradient of polynomial
  ##################################
  
  # read dimension of variables
  L = ncol(degrees)
  
  # initialize list of derivatives
  grad = vector("list", L)
  
  # fill in the derivatives
  for(j in 1:L) {
    # compute derivative
    grad[[j]] = computeDerivative_cpp(coefs, degrees, j, monomialPrimes)
    # remove empty rows
    nonzeros = (grad[[j]]$coefs != 0)
    grad[[j]] = within(grad[[j]], { coefs <- coefs[nonzeros]; degrees <- degrees[nonzeros,] })
  }
  
  # return gradient
  return(grad)
}

rankSymmetricMatrix = function(mat, ereltol=1e-03) {
  
  ####################################################################
  # compute rank of symmetric matrix by counting nonzero eigenvalues
  ####################################################################
  
  # compute eigenvalues
  eig = eigen(mat, symmetric = TRUE, only.values = TRUE)
  # read dimension
  D = length(eig$values)
  # if dimension 1, return 1, else, return rank determined by eigenvalues
  if(D == 1) {
    return(1)
  } else {
    # compute large relative drop in the value of the eigenvalue
    eigRel = with(eig, values[2:D] / values[1:(D-1)]) < ereltol
    # guess rank to be the last location before large relative drop
    return(min(c((1:(D-1))[eigRel], D)))
  }
}

columnEchelon = function(mat) {
  
  ##################################################################
  # compute column echelon form using elementary column operations
  ##################################################################
  
  # read number of columns and rows
  dCol = ncol(mat)
  dRow = nrow(mat)
  
  # initialize echelon form
  ech = mat
  
  # set leading entries to be 1
  for(j in 1:dCol) {
    # read first row where nonzero entry appears
    beginRow = match(TRUE, do.call(pmax, lapply(j:dCol, function(x) abs(ech[,x]))) > 0, nomatch = 0)
    # if everything is zero, exit. If not, continue
    if(beginRow == 0) { break }
    # sort columns according to the descending order of the absolute values
    ech[,j:dCol] = ech[,j-1+order(abs(ech[beginRow,j:dCol]), decreasing=TRUE)]
    # scale the column so that the leading entry is 1
    ech[,j] = ech[,j] / ech[beginRow,j]
    # eliminate other entries in _beginRow_
    for(jPrime in 1:dCol) {
      if(jPrime != j) {
        if(ech[beginRow,jPrime] != 0) {
          ech[,jPrime] = ech[,jPrime] - ech[beginRow,jPrime] * ech[,j]
        }
      }
    }
  }
  
  # return echelon form
  return(ech)
}

createMonomialVector = function(varDim, order) {
  
  ########################################################################################################
  # first create matrix that records lengths of monomial vectors
  #   [i,j] entry represents length of monomial vector of i-th order
  #   which is created by multiplying x_j to the (i-1)-th order monomial vectors starting with x_j^(i-1)
  ########################################################################################################
  
  # initialize matrix
  lengthMatrix = matrix(0, nrow=max(order,1), ncol=varDim)
  
  # fill first orders
  lengthMatrix[1,] = rep(1, varDim)
  
  # fill orders of 2 or higher
  if(order >= 2) {
    for(r in 2:order) {
      for(v in 1:varDim) {
        lengthMatrix[r,v] = sum(lengthMatrix[r-1, v:varDim])
      }
    }
  }
  
  ########################################################################
  # then create a vector of monomials
  #   each variable of the polynomial is identified by the prime numbers
  #   e.g. x_1 is represented by 2, and x_1^3 is represented by 8
  ########################################################################
  
  # we use prime numbers to identify monomials
  # n-th prime is approximately n * log(n). Generate more to be safe
  primes = getPrimes(as.integer(2 + 2 * varDim * log(varDim)))
  
  # initialize monomial vector. We first create a list and then collapse
  monomialVectorList = vector("list", 1+order)
  
  # identifier for the constant term
  monomialVectorList[[1]] = 1
  
  # identifier for the first order terms
  if(order >= 1) {
    monomialVectorList[[2]] = primes[1:varDim]
  }
  
  # record identifiers for the monomials using prime numbers
  if(order >= 2) {
    for(r in 2:order) {
      for(v in 1:varDim) {
        if(v == 1) {
          idVec = primes[v] * monomialVectorList[[r]]
        } else {
          idVec = primes[v] * monomialVectorList[[r]][(1+sum(lengthMatrix[r-1,1:(v-1)])):sum(lengthMatrix[r-1,])]
        }
        monomialVectorList[[1+r]] = c(monomialVectorList[[1+r]], idVec)
      }
    }
  }
  
  # collapse list to vector
  monomialVector = unlist(monomialVectorList, use.names = FALSE)
  
  # return monomial vector along with list of primes used
  return(list(vec=monomialVector, primes=primes[1:varDim], order=order))
}

createMomentMatrix = function(monomialVector) {
  
  ##################################################################
  # create moment matrix, the outer product of the monomial vector
  ##################################################################
  
  return(outer(monomialVector, monomialVector))
}

createMomentMatrixSparse = function(monomialVector) {
  
  ###############################################
  # create moment matrix in sparse triplet form
  ###############################################
  
  # create moment matrix
  momentMatrix = createMomentMatrix(monomialVector)
  
  # transform it into a sparse triplet format
  momentMatrixSparse = data.frame(k=as.vector(row(momentMatrix)), 
                                  l=as.vector(col(momentMatrix)), 
                                  v=as.vector(momentMatrix))
  
  # only keep the lower triangular part since the matrix is symmetric
  momentMatrixSparse = subset(momentMatrixSparse, k >= l)
  
  # return sparse matrix
  return(momentMatrixSparse)
}

restoreOptimalMomentMatrix = function(monomialSystem, optimalMomentVector) {
  
  ############################################################
  # restore optimal moment matrix from optimal moment vector
  ############################################################
  
  # read dimension of moment matrix
  dimMomentMatrix = length(monomialSystem$vec)
  
  # initialize optimal moment matrix
  optimalMomentMatrix = matrix(0, nrow=dimMomentMatrix, ncol=dimMomentMatrix)
  
  # fill in the optimal moment matrix
  optimalMomentMatrix[row(optimalMomentMatrix) >= col(optimalMomentMatrix)] = optimalMomentVector
  optimalMomentMatrix = optimalMomentMatrix + t(optimalMomentMatrix)
  diag(optimalMomentMatrix) = diag(optimalMomentMatrix) / 2
  
  # return optimal moment matrix
  return(optimalMomentMatrix)
}

checkFlatExtension = function(optimalMomentMatrix, varDim, largeOrder, smallOrder, ereltol) {
  
  ##################################################
  # check flat extension condition for certificate
  ##################################################
  
  # read monomial orders
  largeMonomialOrder = as.integer((largeOrder+1)/2)
  smallMonomialOrder = as.integer((smallOrder+1)/2)
  
  # create monomial vectors
  largeMonomialSystem = createMonomialVector(varDim, largeMonomialOrder)
  smallMonomialSystem = createMonomialVector(varDim, smallMonomialOrder)
  
  # read monomial lengths
  largeMonomialLen = length(largeMonomialSystem$vec) # choose(varDim+largeMonomialOrder, varDim)
  smallMonomialLen = length(smallMonomialSystem$vec) # choose(varDim+smallMonomialOrder, varDim)
  
  # extract large and small moment matrices
  largeMatrix = as.matrix(optimalMomentMatrix[1:largeMonomialLen, 1:largeMonomialLen])
  smallMatrix = as.matrix(optimalMomentMatrix[1:smallMonomialLen, 1:smallMonomialLen])
  
  # compute ranks
  largeRank = rankSymmetricMatrix(largeMatrix, ereltol)
  smallRank = rankSymmetricMatrix(smallMatrix, ereltol)
  
  return(list(largeRank=largeRank, smallRank=smallRank, certificate=(largeRank == smallRank)))
}

checkSolutionUnderUniqueness = function(coefs, degrees, objective_primal, objective_dual, optimalMomentMatrix, abstol, reltol) {
  
  #####################################################################
  # compare sdp solution and function value assuming unique optimizer
  #####################################################################
  
  # read dimension of variables
  varDim = ncol(degrees)
  P = nrow(degrees)
  
  # compute function value
  value = evaluatePolynomial_cpp(optimalMomentMatrix[2:(1+varDim),1], coefs, degrees)
  
  # check tolerance condition
  check = (min(absdist(objective_primal, value), absdist(objective_dual, value)) < abstol) + 
    (min(reldist(objective_primal, value), reldist(objective_dual, value)) < reltol)
  
  # return check status
  if(check > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}

addLabelsToSkeleton = function(mosekModel, sense, multithread) {
  
  ####################################################################
  # add model labels to Mosek skeleton object created by cppfunction
  ####################################################################
  
  mosekModel$sense = sense
  mosekModel$A = Matrix(mosekModel$A) # depends on Matrix package
  rownames(mosekModel$bc) = c("blc", "buc")
  rownames(mosekModel$bx) = c("blx", "bux")
  if(multithread == TRUE) {
    mosekModel$iparam = list(MSK_IPAR_INTPNT_MULTI_THREAD = "MSK_ON") # turn on multithreading
  } else {
    mosekModel$iparam = list(MSK_IPAR_INTPNT_MULTI_THREAD = "MSK_OFF") # turn off multithreading
  }
  return(mosekModel)
  
}

sdpmodel = function(sense, coefs, degrees, opt=NULL) {
  
  #################################################
  # creates sdp model for polynomial optimization
  #################################################
  
  # read options list
  if(is.null(opt) == TRUE) {
    options = list()
  } else {
    options = opt
  }
  
  # read radius
  if(is.null(options$radius) == TRUE) {
    constrained = 0
    radius = 0
  } else {
    constrained = 1
    radius = options$radius
  }
  
  # read max length of hierarchy
  if(is.null(options$hierarchy) == TRUE) {
    hierarchy = 1
  } else {
    hierarchy = options$hierarchy
  }
  
  # read multithreading option
  if(is.null(options$multithread) == TRUE) {
    multithread = FALSE
  } else {
    multithread = options$multithread
  }
  
  # read the dimension of the variables and the order of polynomial
  varDim = ncol(degrees)
  order = max(rowSums(degrees))
  
  # initialize list of models
  models = vector("list", hierarchy)
  
  if(constrained == 0) {
    
    ########## if unconstrained, implement gradient ideal method
    
    # return error if order is odd and the problem is unconstrained
    if(order %% 2 != 0 & constrained == 0) {
      stop("order of the polynomial must be multiples of two in unconstrained optimization")
    }
    
    # create monomial vector of half of nearest-even order,
    monomialSystem = createMonomialVector(varDim, as.integer((order+1)/2))
    
    # compute gradient
    grad = computeGradient(coefs, degrees, monomialSystem$primes)
    
    # solve hierarchy of SDPs
    for(j in 1:hierarchy) {
      
      # set order at the hierarchy
      orderHierarchy = order + 2 * (j-1)
      
      # create monomial vector and moment matrix
      #   first create a monomial vector of half of nearest-even order,
      #   which creates moment matrix of the nearest-even order
      monomialSystem = createMonomialVector(varDim, as.integer((orderHierarchy+1)/2))
      momentMatrixSparse = as.matrix(createMomentMatrixSparse(monomialSystem$vec))
      
      # create mosek model skeleton
      model = createMosekSdpModelSkeletonWithGradientIdeals_cpp(varDim, order, orderHierarchy, grad, 
                                                                monomialSystem$primes, momentMatrixSparse)
      
      # add labels to skeleton
      model = addLabelsToSkeleton(model, sense, multithread)
      
      # record objective polynomial coefficients
      model$barc = createMosekSdpCoefficientMatrixFromDegrees_cpp(coefs, degrees, monomialSystem$primes, momentMatrixSparse)
      
      # save model
      models[[j]] = list(model=model, monomialSystem=monomialSystem)
      
    }
    
  } else {
    
    ########## if constrained, implement Lasserre's localizing matrix method
    
    # solve hierarchy of SDPs
    for(j in 1:hierarchy) {
      
      # set order at the hierarchy
      orderHierarchy = order + 2 * (j-1)
      
      # create monomial vector and moment matrix
      #   first create a monomial vector of half of nearest-even order,
      #   which creates moment matrix of the nearest-even order
      monomialSystem = createMonomialVector(varDim, as.integer((orderHierarchy+1)/2))
      momentMatrixSparse = as.matrix(createMomentMatrixSparse(monomialSystem$vec))
      
      # create mosek model skeleton
      model = createMosekSdpModelSkeleton_cpp(varDim, orderHierarchy, constrained, radius, monomialSystem$primes, momentMatrixSparse)
      
      # add labels to skeleton
      model = addLabelsToSkeleton(model, sense, multithread)
      
      # record objective polynomial coefficients
      model$barc = createMosekSdpCoefficientMatrixFromDegrees_cpp(coefs, degrees, monomialSystem$primes, momentMatrixSparse)
      
      # save model
      models[[j]] = list(model=model, monomialSystem=monomialSystem)
      
    }
    
  }
  
  # return models
  return(models)
}

optpoly = function(sense, coefs, degrees, opt=NULL) {
  
  ########################################
  # function for polynomial optimization
  ########################################
  
  # read options list
  if(is.null(opt) == TRUE) {
    options = list()
  } else {
    options = opt
  }
  
  # read radius
  if(is.null(options$radius) == TRUE) {
    constrained = 0
    radius = 0
  } else {
    constrained = 1
    radius = options$radius
  }
  
  # read max length of hierarchy
  if(is.null(options$hierarchy) == TRUE) {
    hierarchy = 1
  } else {
    hierarchy = options$hierarchy
  }
  
  # read multithreading option
  if(is.null(options$multithread) == TRUE) {
    multithread = FALSE
  } else {
    multithread = options$multithread
  }
  
  # read verbose option
  if(is.null(options$verbose) == TRUE) {
    verbose = 0
  } else {
    verbose = options$verbose
  }
  
  # absolute tolerance for simple certificate under uniqueness
  if(is.null(options$fabstol) == TRUE) {
    fabstol = 1e-10
  } else {
    fabstol = options$fabstol
  }
  
  # relative tolerance for simple certificate under uniqueness
  if(is.null(options$freltol) == TRUE) {
    freltol = 1e-06
  } else {
    freltol = options$freltol
  }
  
  # read relative eigenvalue tolerance for rank check
  if(is.null(options$ereltol) == TRUE) {
    ereltol = 1e-03
  } else {
    ereltol = options$ereltol
  }
  
  # read the dimension of the variables and the order of polynomial
  varDim = ncol(degrees)
  order = max(rowSums(degrees))
  
  # create sdp models
  models = sdpmodel(sense, coefs, degrees, opt = opt)
  
  if(constrained == 0) {
    
    ########## if unconstrained, implement gradient ideal method
    
    # return error if order is odd and the problem is unconstrained
    if(order %% 2 != 0 & constrained == 0) {
      stop("order of the polynomial must be multiples of two in unconstrained optimization")
    }
    
    # solve hierarchy of SDPs
    for(j in 1:hierarchy) {
      
      # solve SDP
      mosekSol = mosek(models[[j]]$model, opts = list(verbose=verbose, soldetail=1))
      
      # check if SDP returned error
      checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
      
      # process solution according to the error status
      if(checkError) {
        
        sol = list(sdpstatus = mosekSol$response$msg,
                   varDim    = varDim,
                   order     = order,
                   hierarchy = j)
        
      } else {
        
        # construct optimal moment matrix
        optimalMomentMatrix = restoreOptimalMomentMatrix(models[[j]]$monomialSystem, mosekSol$sol$itr$barx[[1]])
        
        # check flat extension condition
        check = checkFlatExtension(optimalMomentMatrix, varDim, models[[j]]$monomialSystem$order, models[[j]]$monomialSystem$order-order, ereltol)
        
        # if condition holds or last hierarchy, record solution
        if(check$certificate > 0 | j == hierarchy) {
          sol = list(objective_primal = mosekSol$sol$itr$pobjval,
                     objective_dual   = mosekSol$sol$itr$dobjval,
                     sdpstatus        = mosekSol$response$msg,
                     solstatus        = mosekSol$sol$itr$solsta,
                     momentmatrix     = optimalMomentMatrix,
                     varDim           = varDim,
                     order            = order)
          sol = c(sol, check, hierarchy=j)
          break
        }
        
      }
      
    }
    
  } else {
    
    ########## if constrained, implement Lasserre's localizing matrix method
    
    # solve hierarchy of SDPs
    for(j in 1:hierarchy) {
      
      # solve SDP
      mosekSol = mosek(models[[j]]$model, opts = list(verbose=verbose, soldetail=1))
      
      # check if SDP returned error
      checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
      
      # process solution according to the error status
      if(checkError) {
        
        sol = list(sdpstatus = mosekSol$response$msg,
                   varDim    = varDim,
                   order     = order,
                   hierarchy = j)
        
      } else {
        
        # construct optimal moment matrix
        optimalMomentMatrix = restoreOptimalMomentMatrix(models[[j]]$monomialSystem, mosekSol$sol$itr$barx[[1]])
        
        # check flat extension condition
        check = checkFlatExtension(optimalMomentMatrix, varDim, models[[j]]$monomialSystem$order, models[[j]]$monomialSystem$order-order, ereltol)
        
        # if condition holds or last hierarchy, record solution
        if(check$certificate > 0 | j == hierarchy) {
          sol = list(objective_primal = mosekSol$sol$itr$pobjval,
                     objective_dual   = mosekSol$sol$itr$dobjval,
                     sdpstatus        = mosekSol$response$msg,
                     solstatus        = mosekSol$sol$itr$solsta,
                     momentmatrix     = optimalMomentMatrix,
                     varDim           = varDim,
                     order            = order)
          sol = c(sol, check, hierarchy=j)
          break
        }
        
      }
      
    }
    
  }
  
  # return solution
  return(sol)
}

optmanypoly = function(sense, coefMatrix, varDim, orderObj, orderMom, constrained, radius=NULL, multithread, verbose) {
  
  ######################################################################################################
  # run optimizations of many polynomials under the same SDP configuration
  # each row of coefMatrix, which has the same length as momentMatrixSparse, represent each polynomial
  # coefficients in the rows of coefMatrix must be arranged in the same order as momentMatrixSparse
  ######################################################################################################
  
  # read number of polynomials
  N = nrow(coefMatrix)
  
  # produce error if N / multithread < 2
  if( (multithread > 1) & (N / multithread < 2) ) {
    stop("each thread must process at least two polynomials when multithreading")
  }
  
  # create monomial system and moment matrix that corresponds to coefMatrix
  monomialSystemObj = createMonomialVector(varDim, as.integer((orderObj+1)/2))
  momentMatrixSparseObj = as.matrix(createMomentMatrixSparse(monomialSystemObj$vec))
  
  # create monomial system and moment matrix for the sdp
  monomialSystemMom = createMonomialVector(varDim, as.integer((orderMom+1)/2))
  momentMatrixSparseMom = as.matrix(createMomentMatrixSparse(monomialSystemMom$vec))
  
  # transform monomials to degrees
  degrees = monomialsToDegrees_cpp(momentMatrixSparseObj[,3], monomialSystemObj$primes)
  
  # partition individuals according to the number of threads
  # this array looks like: c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,...)
  partitionIndex = 1 + as.integer(0:(N-1) * (multithread / N))
  
  if(constrained == 0) {
    
    ########## if unconstrained, implement gradient ideal method
    
    # if singlethread, use plain for loop
    if(multithread == 1) {
      
      # create data frame that stores results
      res = data.frame(objective_primal = rep(NA, nrow(coefMatrix)),
                       objective_dual   = NA,
                       sdpstatus        = " ",
                       solstatus        = " ",
                       matrix(NA, nrow=nrow(coefMatrix), ncol=nrow(momentMatrixSparseMom)), stringsAsFactors = FALSE)
      
      # fill in the data frame
      for(i in 1:nrow(coefMatrix)) {
        
        # compute gradient
        grad = computeGradient(coefMatrix[i,], degrees, monomialSystemObj$primes)
        
        # configure SDP problem
        model = createMosekSdpModelSkeletonWithGradientIdeals_cpp(varDim, orderObj, orderMom, grad, 
                                                                  monomialSystemMom$primes, momentMatrixSparseMom)
        
        # add labels to skeleton
        model = addLabelsToSkeleton(model, sense, FALSE)
        
        # record objective polynomial coefficients
        model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrix[i,], momentMatrixSparseObj[,3], momentMatrixSparseObj)
        
        # solve SDP
        mosekSol = mosek(model, opts = list(verbose=verbose, soldetail=1))
        
        # check if SDP returned error
        checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
        
        # process solution according to the error status
        if(checkError) {
          
          res$sdpstatus[i]        = mosekSol$response$msg
          
        } else {
          
          # record solutions
          res$objective_primal[i] = mosekSol$sol$itr$pobjval
          res$objective_dual[i]   = mosekSol$sol$itr$dobjval
          res$sdpstatus[i]        = mosekSol$response$msg
          res$solstatus[i]        = mosekSol$sol$itr$solsta
          res[i,-c(1:4)]          = mosekSol$sol$itr$barx[[1]]
          
        }
        
      }
      
    } else {
      
      # if multithread, use foreach dopar
      res = foreach(h=1:multithread, .combine = "rbind") %dopar% {
        
        # subset coefficient matrix
        coefMatrixSub = coefMatrix[partitionIndex == h,]
        
        # create data frame that stores results
        mosekSolDF = data.frame(objective_primal = rep(NA, nrow(coefMatrixSub)),
                                objective_dual   = NA,
                                sdpstatus        = " ",
                                solstatus        = " ",
                                matrix(NA, nrow=nrow(coefMatrixSub), ncol=nrow(momentMatrixSparseMom)), stringsAsFactors = FALSE)
        
        # fill in the data frame
        for(i in 1:nrow(coefMatrixSub)) {
          
          # compute gradient
          grad = computeGradient(coefMatrixSub[i,], degrees, monomialSystemObj$primes)
          
          # configure SDP problem
          model = createMosekSdpModelSkeletonWithGradientIdeals_cpp(varDim, orderObj, orderMom, grad, 
                                                                    monomialSystemMom$primes, momentMatrixSparseMom)
          
          # add labels to skeleton
          model = addLabelsToSkeleton(model, sense, FALSE)
          
          # record objective polynomial coefficients
          model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrixSub[i,], momentMatrixSparseObj[,3], momentMatrixSparseObj)
          
          # solve SDP
          mosekSol = mosek(model, opts = list(verbose=verbose, soldetail=1))
          
          # check if SDP returned error
          checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
          
          # process solution according to the error status
          if(checkError) {
            
            mosekSolDF$sdpstatus[i]        = mosekSol$response$msg
            
          } else {
            
            # record solutions
            mosekSolDF$objective_primal[i] = mosekSol$sol$itr$pobjval
            mosekSolDF$objective_dual[i]   = mosekSol$sol$itr$dobjval
            mosekSolDF$sdpstatus[i]        = mosekSol$response$msg
            mosekSolDF$solstatus[i]        = mosekSol$sol$itr$solsta
            mosekSolDF[i,-c(1:4)]          = mosekSol$sol$itr$barx[[1]]
            
          }
          
        }
        
        # return data frame
        mosekSolDF
        
      }
      
    }
    
  } else {
    
    ########## if constrained, implement Lasserre's localizing matrix method
    
    # configure SDP problem
    model = createMosekSdpModelSkeleton_cpp(varDim, orderMom, constrained, radius, monomialSystemMom$primes, momentMatrixSparseMom)
    
    # add labels to skeleton
    model = addLabelsToSkeleton(model, sense, FALSE)
    
    # if singlethread, use plain for loop
    if(multithread == 1) {
      
      # create data frame that stores results
      res = data.frame(objective_primal = rep(NA, nrow(coefMatrix)),
                       objective_dual   = NA,
                       sdpstatus        = " ",
                       solstatus        = " ",
                       matrix(NA, nrow=nrow(coefMatrix), ncol=nrow(momentMatrixSparseMom)), stringsAsFactors = FALSE)
      
      # fill in the data frame
      for(i in 1:nrow(coefMatrix)) {
        
        # record objective polynomial coefficients
        model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrix[i,], momentMatrixSparseObj[,3], momentMatrixSparseObj)
        
        # solve SDP
        mosekSol = mosek(model, opts = list(verbose=verbose, soldetail=1))
        
        # check if SDP returned error
        checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
        
        # process solution according to the error status
        if(checkError) {
          
          res$sdpstatus[i]        = mosekSol$response$msg
          
        } else {
          
          # record solutions
          res$objective_primal[i] = mosekSol$sol$itr$pobjval
          res$objective_dual[i]   = mosekSol$sol$itr$dobjval
          res$sdpstatus[i]        = mosekSol$response$msg
          res$solstatus[i]        = mosekSol$sol$itr$solsta
          res[i,-c(1:4)]          = mosekSol$sol$itr$barx[[1]]
          
        }
        
      }
      
    } else {
      
      res = foreach(h=1:multithread, .combine = "rbind") %dopar% {
        
        # subset coefficient matrix
        coefMatrixSub = coefMatrix[partitionIndex == h,]
        
        # create data frame that stores results
        mosekSolDF = data.frame(objective_primal = rep(NA, nrow(coefMatrixSub)),
                                objective_dual   = NA,
                                sdpstatus        = " ",
                                solstatus        = " ",
                                matrix(NA, nrow=nrow(coefMatrixSub), ncol=nrow(momentMatrixSparseMom)), stringsAsFactors = FALSE)
        
        # fill in the data frame
        for(i in 1:nrow(coefMatrixSub)) {
          
          # record objective polynomial coefficients
          model$barc = createMosekSdpCoefficientMatrixFromMonomials_cpp(coefMatrixSub[i,], momentMatrixSparseObj[,3], momentMatrixSparseObj)
          
          # solve SDP
          mosekSol = mosek(model, opts = list(verbose=verbose, soldetail=1))
          
          # check if SDP returned error
          checkError = (substr(mosekSol$response$msg, 1, 11) == "MSK_RES_ERR")
          
          # process solution according to the error status
          if(checkError) {
            
            mosekSolDF$sdpstatus[i]        = mosekSol$response$msg
            
          } else {
            
            # record solutions
            mosekSolDF$objective_primal[i] = mosekSol$sol$itr$pobjval
            mosekSolDF$objective_dual[i]   = mosekSol$sol$itr$dobjval
            mosekSolDF$sdpstatus[i]        = mosekSol$response$msg
            mosekSolDF$solstatus[i]        = mosekSol$sol$itr$solsta
            mosekSolDF[i,-c(1:4)]          = mosekSol$sol$itr$barx[[1]]
            
          }
          
        }
        
        # return data frame
        mosekSolDF
        
      }
      
    }
    
  }
  
  # return results
  return(list(objective_primal = res$objective_primal,
              objective_dual   = res$objective_dual,
              sdpstatus        = factor(res$sdpstatus),
              solstatus        = factor(res$solstatus),
              varDim           = varDim,
              orderobj         = orderObj,
              orderMom         = orderMom,
              moment_vectors   = res[,-c(1:4)]))
}

extractSolution = function(sol, points=NULL, vreltol=1e-10, ereltol=1e-03) {
  
  # read varDim
  varDim = sol$varDim
  # read order of monomial vector
  monomialVectorOrder = as.integer((sol$order+1)/2)
  
  # create monomial vector
  monomialSystem = createMonomialVector(varDim, monomialVectorOrder)
  
  # read optimal moment matrix
  optimalMomentMatrix = sol$momentmatrix
  
  # read number of solutions from rank
  if(is.null(points) == TRUE) {
    nSol = sol$smallRank
  } else {
    nSol = points
  }
  
  # compute left-matrix of SVD of optimal moment matrix
  if(nSol == 1) {
    singularMatrix = with(svd(optimalMomentMatrix, nu = nSol, nv = 0), u * d[1])
  } else {
    singularMatrix = with(svd(optimalMomentMatrix, nu = nSol, nv = 0), u %*% diag(d[1:nSol]))
  }
  
  # set infinitesimal values to be zero
  singularMatrix[abs(singularMatrix) < vreltol] = 0
  
  # compute column echelon form
  echelon = columnEchelon(singularMatrix)
  
  # extract basis polynomials of the solution
  basisMonomials = rep(0, nSol)
  for(j in 1:nSol) {
    basisMonomials[j] = monomialSystem$vec[min((1:nrow(echelon))[echelon[,j] == 1])]
  }
  
  # compute multiplication matrices
  multipMatrixList = vector("list", varDim)
  for(n in 1:varDim) {
    multipMatrix = matrix(0, nrow=nSol, ncol=nSol)
    for(j in 1:nSol) {
      multipMatrix[j,] = echelon[match(basisMonomials[j] * monomialSystem$primes[n], monomialSystem$vec),]
    }
    multipMatrixList[[n]] = multipMatrix
  }
  
  # take random average of multiplication matrices
  weights = runif(varDim, min=1, max=2)
  weights = weights / sum(weights)
  aggregateMultipMatrix = matrix(0, nrow=nSol, ncol=nSol)
  for(n in 1:varDim) {
    aggregateMultipMatrix = aggregateMultipMatrix + weights[n] * multipMatrixList[[n]]
  }
  
  # compute Schur decomposition and take Q
  qMat = Schur(aggregateMultipMatrix, vectors = TRUE)$Q
  
  # record solutions
  argMat = matrix(0, nrow=nSol, ncol=varDim)
  for(n in 1:varDim) {
    for(j in 1:nSol) {
      argMat[j,n] = qMat[,j] %*% multipMatrixList[[n]] %*% qMat[,j]
    }
  }
  
  # return solution
  return(argMat)
}
