
getPrimes = function(n) {

  ###################################################################################################
  # return list of prime numbers
  # e.g. getPrimes(100) = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
  ###################################################################################################

  return(getPrimes_cpp(n))
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

  # # read number of columns and rows
  # dCol = ncol(mat)
  # dRow = nrow(mat)
  # 
  # # initialize echelon form
  # ech = mat
  # 
  # # set leading entries to be 1
  # for(j in 1:dCol) {
  #   # read first row where nonzero entry appears
  #   beginRow = match(TRUE, do.call(pmax, lapply(j:dCol, function(x) abs(ech[,x]))) > 0, nomatch = 0)
  #   # if everything is zero, exit. If not, continue
  #   if(beginRow == 0) { break }
  #   # sort columns according to the descending order of the absolute values
  #   ech[,j:dCol] = ech[,j-1+order(abs(ech[beginRow,j:dCol]), decreasing=TRUE)]
  #   # scale the column so that the leading entry is 1
  #   ech[,j] = ech[,j] / ech[beginRow,j]
  #   # eliminate other entries in _beginRow_
  #   for(jPrime in 1:dCol) {
  #     if(jPrime != j) {
  #       if(ech[beginRow,jPrime] != 0) {
  #         ech[,jPrime] = ech[,jPrime] - ech[beginRow,jPrime] * ech[,j]
  #       }
  #     }
  #   }
  # }
  # 
  # # return echelon form
  # return(ech)
  
  # Above operations are wrapped in the function columnEchelon_cpp
  return(columnEchelon_cpp(mat))
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

  # transform into sparse triplet format
  # only keep the lower triangular part since the matrix is symmetric
  momentMatrixSparse = matrixToSparseTriplet_cpp(momentMatrix, 1)

  # label columns
  colnames(momentMatrixSparse) = c("k","l","v")

  # return sparse matrix
  return(momentMatrixSparse)
}

restoreOptimalMomentMatrix = function(monomialSystem, optimalMomentVector) {

  ############################################################
  # restore optimal moment matrix from optimal moment vector
  ############################################################
  
  # transform to optimal moment matrix and return
  return(vecToMatrix_cpp(optimalMomentVector, length(monomialSystem$vec), 1))
}

checkFlatExtension = function(optimalMomentMatrix, varDim, largeOrder, smallOrder, ereltol) {

  ##################################################
  # check flat extension condition for certificate
  ##################################################

  # # create monomial vectors
  # largeMonomialSystem = createMonomialVector(varDim, largeOrder)
  # smallMonomialSystem = createMonomialVector(varDim, smallOrder)
  # # read monomial lengths
  # largeMonomialLen = length(largeMonomialSystem$vec) # choose(varDim+largeMonomialOrder, varDim)
  # smallMonomialLen = length(smallMonomialSystem$vec) # choose(varDim+smallMonomialOrder, varDim)
  
  # skip above steps and directly compute monomial length by formula
  largeMonomialLen = choose(varDim+largeOrder, varDim)
  smallMonomialLen = choose(varDim+smallOrder, varDim)
  
  # extract large and small moment matrices
  largeMatrix = optimalMomentMatrix[1:largeMonomialLen, 1:largeMonomialLen]
  smallMatrix = optimalMomentMatrix[1:smallMonomialLen, 1:smallMonomialLen]

  # compute ranks
  largeRank = rankSymmetricMatrix(largeMatrix, ereltol)
  smallRank = rankSymmetricMatrix(smallMatrix, ereltol)

  return(list(largeRank=largeRank, smallRank=smallRank, certificate=(largeRank == smallRank)))
}

checkFlatExtensions = function(optimalMomentMatrix, varDim, order, rankStep, ereltol) {

  # read order of moment matrix
  momentMatrixOrder = as.integer((order+1)/2)

  # check rank condition
  for(r in 1:(momentMatrixOrder-rankStep)) {
    check = checkFlatExtension(optimalMomentMatrix, varDim, r+rankStep, r, ereltol)
    if(check$certificate == TRUE) {
      return(list(certificate=check$certificate, rank=check$smallRank))
    }
  }

  # if rank condition not satisfied, return FALSE and return smallRank of last rank check
  return(list(certificate=FALSE, rank=check$smallRank))
}

checkEvaluation = function(coefs, degrees, objective_primal, objective_dual, optimizers, abstol, reltol) {

  #########################################################################
  # compare sdp solution and function value given the number of solutions
  #########################################################################

  # read dimension of variables
  varDim = ncol(degrees)

  # compute function values
  value = rep(0, nrow(optimizers))
  for(k in 1:nrow(optimizers)) {
    value[k] = evaluatePolynomial_cpp(optimizers[k,], coefs, degrees)
  }

  # check tolerance condition
  checks = rep(0, nrow(optimizers))
  for(k in 1:nrow(optimizers)) {
    checks[k] = (min(absdist(objective_primal, value[k]), absdist(objective_dual, value[k])) < abstol) +
      (min(reldist(objective_primal, value[k]), reldist(objective_dual, value[k])) < reltol)
  }

  # return check status
  if(sum(checks) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}

checkCertificate = function(varDim, orderMom, coefs, degrees, objective_primal, objective_dual, optimalMomentMatrix, options) {
  
  # first, perform evaluation test
  
  # extract optimizer assuming unique solution
  uniqueOptimum = optimalMomentMatrix[-1, 1]
  uniqueOptimum = matrix(uniqueOptimum, nrow=1, ncol=length(uniqueOptimum))
  
  # check evaluation
  checkEval = checkEvaluation(coefs, degrees, objective_primal, objective_dual, uniqueOptimum, options$fabstol, options$freltol)
  
  # quit if pass evaluation test
  if(checkEval == TRUE) { return(list(certificate = TRUE, rank = 1)) }
  
  # else, proceed to flat extension test
  # max order is 1 for unconstrained optimization
  # max order of g_j is also 1 for contrained optimization with the constraint radius^2 - sum_j x_j^2 >= 0
  checkFlat = checkFlatExtensions(optimalMomentMatrix, varDim, orderMom, 1, options$ereltol)
  
  # quit if pass flat extension test
  if(checkFlat$certificate == TRUE) { return(checkFlat) }
  
  # else, perform evaluation test using moment matrix rank as number of optimizers
  
  # extract optimizers
  multipleOptima = extractSolution(list(varDim=varDim, order=orderMom, momentmatrix=optimalMomentMatrix), 
                                   checkFlat$rank, options$vabstol, options$ereltol)
  
  checkEvalMultiple = checkEvaluation(coefs, degrees, objective_primal, objective_dual, multipleOptima, options$fabstol, options$freltol)
  
  # return final result
  if(checkEvalMultiple == TRUE) {
    return(list(certificate = TRUE, rank = checkFlat$rank))
  } else {
    return(checkFlat)
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

createOptionsList = function(opt=NULL) {
  
  # read options list
  if(is.null(opt) == TRUE) {
    options = list()
  } else {
    options = opt
  }
  
  # read bounds or radius
  if(is.null(options$bounds) == TRUE) {
    if(is.null(options$radius) == TRUE) {
      options$constrained = 0
      options$rectangular = 0
      options$bounds = matrix(0, nrow=1, ncol=1)
      options$radius = 0
    } else {
      options$constrained = 1
      options$rectangular = 0
      options$bounds = matrix(0, nrow=1, ncol=1)
      # options$radius = opt$radius
    }
  } else {
    options$constrained = 1
    options$rectangular = 1
    # options$bounds = opt$bounds
    if(is.null(options$radius) == TRUE) {
      options$radius = max(abs(options$bounds))
    } else {
      # options$radius = opt$radius
    }
  }
  
  # read max length of hierarchy
  if(is.null(options$hierarchy) == TRUE) {
    options$hierarchy = 1
  } else {
    # options$hierarchy = opt$hierarchy
  }
  
  # read multithreading option
  if(is.null(options$multithread) == TRUE) {
    options$multithread = FALSE
  } else {
    # multithread = options$multithread
  }
  
  # read verbose option
  if(is.null(options$verbose) == TRUE) {
    options$verbose = 0
  } else {
    # verbose = options$verbose
  }
  
  # absolute tolerance for simple certificate from function evaluations
  if(is.null(options$fabstol) == TRUE) {
    options$fabstol = 1e-10
  } else {
    # fabstol = options$fabstol
  }
  
  # relative tolerance for simple certificate from function evaluations
  if(is.null(options$freltol) == TRUE) {
    options$freltol = 1e-06
  } else {
    # freltol = options$freltol
  }
  
  # absolute tolerance for zero value
  if(is.null(options$vabstol) == TRUE) {
    options$vabstol = 1e-10
  } else {
    # vabstol = options$vabstol
  }
  
  # relative eigenvalue tolerance for rank check
  if(is.null(options$ereltol) == TRUE) {
    options$ereltol = 1e-03
  } else {
    # ereltol = options$ereltol
  }
  
  # return options list
  return(options)
}

sdpmodel = function(sense, coefs, degrees, opt=NULL) {

  #################################################
  # creates sdp model for polynomial optimization
  #################################################

  # read options list
  options = createOptionsList(opt)

  # read the dimension of the variables and the order of polynomial
  varDim = ncol(degrees)
  order = max(rowSums(degrees))

  # initialize list of models
  models = vector("list", options$hierarchy)

  if(options$constrained == 0) {

    ########## if unconstrained, implement Nie et al's gradient ideal method

    # return error if order is odd and the problem is unconstrained
    if(order %% 2 != 0 & options$constrained == 0) {
      stop("order of the polynomial must be multiples of two in unconstrained optimization")
    }

    # create monomial vector of half of nearest-even order,
    monomialSystem = createMonomialVector(varDim, as.integer((order+1)/2))

    # compute gradient
    grad = computeGradient(coefs, degrees, monomialSystem$primes)

    # solve hierarchy of SDPs
    for(j in 1:(options$hierarchy)) {

      # set order at the hierarchy
      orderHierarchy = order + 2 * (j-1)

      # create monomial vector and moment matrix
      #   first create a monomial vector of half of nearest-even order,
      #   which creates moment matrix of the nearest-even order
      monomialSystem = createMonomialVector(varDim, as.integer((orderHierarchy+1)/2))
      momentMatrixSparse = as.matrix(createMomentMatrixSparse(monomialSystem$vec))

      # create mosek model skeleton
      model = createMosekSdpModelSkeletonNieetal_cpp(varDim, order, orderHierarchy, grad,
                                                     monomialSystem$primes, momentMatrixSparse)

      # add labels to skeleton
      model = addLabelsToSkeleton(model, sense, options$multithread)

      # record objective polynomial coefficients
      model$barc = createMosekSdpCoefficientMatrixFromDegrees_cpp(coefs, degrees, monomialSystem$primes, momentMatrixSparse)

      # save model
      models[[j]] = list(model=model, monomialSystem=monomialSystem)

    }

  } else {

    ########## if constrained, implement Lasserre's localizing matrix method

    # solve hierarchy of SDPs
    for(j in 1:(options$hierarchy)) {

      # set order at the hierarchy
      orderHierarchy = order + 2 * (j-1)

      # create monomial vector and moment matrix
      #   first create a monomial vector of half of nearest-even order,
      #   which creates moment matrix of the nearest-even order
      monomialSystem = createMonomialVector(varDim, as.integer((orderHierarchy+1)/2))
      momentMatrixSparse = as.matrix(createMomentMatrixSparse(monomialSystem$vec))

      # create mosek model skeleton
      model = createMosekSdpModelSkeletonLasserre_cpp(varDim, orderHierarchy, options$constrained, options$rectangular, options$bounds,
                                                      options$radius, monomialSystem$primes, momentMatrixSparse)

      # add labels to skeleton
      model = addLabelsToSkeleton(model, sense, options$multithread)

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
  options = createOptionsList(opt)

  # read the dimension of the variables and the order of polynomial
  varDim = ncol(degrees)
  order = max(rowSums(degrees))

  # create sdp models
  models = sdpmodel(sense, coefs, degrees, opt = opt)

  if(options$constrained == 0) {

    ########## if unconstrained, implement gradient ideal method by Nie et al

    # return error if order is odd and the problem is unconstrained
    if(order %% 2 != 0) {
      stop("order of the polynomial must be multiples of two in unconstrained optimization")
    }

    # solve hierarchy of SDPs
    for(j in 1:(options$hierarchy)) {

      # solve SDP
      mosekSol = mosek(models[[j]]$model, opts = list(verbose=options$verbose, soldetail=1))

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

        # # check flat extension condition
        # # max order is 1 for unconstrained optimization
        # # max order of g_j is also 1 for contrained optimization with the constraint radius^2 - sum_j x_j^2 >= 0
        # check = checkFlatExtensions(optimalMomentMatrix, varDim, models[[j]]$monomialSystem$order*2, 1, ereltol)
        # 
        # # check (heuristic) function evaluation criterion
        # # first assume unique optimizer
        # if(check$certificate == FALSE) {
        #   check_eval = checkEvaluation(coefs, degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
        #                                extractSolution(list(varDim=varDim, order=order, momentmatrix=optimalMomentMatrix), 1, vabstol, ereltol),
        #                                fabstol, freltol)
        #   # overwrite certificate if function evaluation criterion is met
        #   if(check_eval == TRUE) { check$certificate = TRUE; check$rank = 1 }
        # }
        # # now use check$rank as number of solutions
        # if(check$certificate == FALSE) {
        #   check_eval = checkEvaluation(coefs, degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
        #                                extractSolution(list(varDim=varDim, order=order, momentmatrix=optimalMomentMatrix), check$rank, vabstol, ereltol),
        #                                fabstol, freltol)
        #   # overwrite certificate if function evaluation criterion is met
        #   if(check_eval == TRUE) { check$certificate = TRUE }
        # }
        
        # check certificate. Above operations are wrapped in the function checkCertificate
        check = checkCertificate(varDim, models[[j]]$monomialSystem$order*2, coefs, degrees, 
                                 mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval, optimalMomentMatrix, options)

        # if certificate obtained or last hierarchy, record solution
        if(check$certificate > 0 | j == options$hierarchy) {
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
    for(j in 1:(options$hierarchy)) {

      # solve SDP
      mosekSol = mosek(models[[j]]$model, opts = list(verbose=options$verbose, soldetail=1))

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

        # # check flat extension condition
        # # max order is 1 for unconstrained optimization
        # # max order of g_j is also 1 for contrained optimization with the constraint radius^2 - sum_j x_j^2 >= 0
        # check = checkFlatExtensions(optimalMomentMatrix, varDim, models[[j]]$monomialSystem$order*2, 1, ereltol)
        # 
        # # check (heuristic) function evaluation criterion
        # # first assume unique optimizer
        # if(check$certificate == FALSE) {
        #   check_eval = checkEvaluation(coefs, degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
        #                                extractSolution(list(varDim=varDim, order=order, momentmatrix=optimalMomentMatrix), 1, vabstol, ereltol),
        #                                fabstol, freltol)
        #   # overwrite certificate if function evaluation criterion is met
        #   if(check_eval == TRUE) { check$certificate = TRUE; check$rank = 1 }
        # }
        # # now use check$rank as number of solutions
        # if(check$certificate == FALSE) {
        #   check_eval = checkEvaluation(coefs, degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
        #                                extractSolution(list(varDim=varDim, order=order, momentmatrix=optimalMomentMatrix), check$rank, vabstol, ereltol),
        #                                fabstol, freltol)
        #   # overwrite certificate if function evaluation criterion is met
        #   if(check_eval == TRUE) { check$certificate = TRUE }
        # }
        
        # check certificate. Above operations are wrapped in the function checkCertificate
        check = checkCertificate(varDim, models[[j]]$monomialSystem$order*2, coefs, degrees, 
                                 mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval, optimalMomentMatrix, options)

        # if certificate obtained or last hierarchy, record solution
        if(check$certificate > 0 | j == options$hierarchy) {
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

optmanypoly = function(sense, coefMatrix, varDim, orderObj, orderMom, opt=NULL) {

  ######################################################################################################
  # run optimizations of many polynomials under the same SDP configuration
  # each row of coefMatrix, which has the same length as momentMatrixSparse, represent each polynomial
  # coefficients in the rows of coefMatrix must be arranged in the same order as momentMatrixSparse
  ######################################################################################################

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

  # read multithreading option
  if(is.null(options$multithread) == TRUE) {
    multithread = 1
  } else {
    multithread = options$multithread
  }

  # read verbose option
  if(is.null(options$verbose) == TRUE) {
    verbose = 0
  } else {
    verbose = options$verbose
  }

  # absolute tolerance for simple certificate from function evaluations
  if(is.null(options$fabstol) == TRUE) {
    fabstol = 1e-10
  } else {
    fabstol = options$fabstol
  }

  # relative tolerance for simple certificate from function evaluations
  if(is.null(options$freltol) == TRUE) {
    freltol = 1e-06
  } else {
    freltol = options$freltol
  }

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

    ########## if unconstrained, implement Nie et al's gradient ideal method

    # if singlethread, use plain for loop
    if(multithread == 1) {

      # create data frame that stores results
      res = data.frame(objective_primal = rep(NA, nrow(coefMatrix)),
                       objective_dual   = NA,
                       sdpstatus        = " ",
                       solstatus        = " ",
                       certificate      = FALSE,
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
          res$certificate[i]      = checkEvaluation(coefMatrix[i,], degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
                                                    matrix(mosekSol$sol$itr$barx[[1]][1+1:varDim], nrow=1),
                                                    fabstol, freltol)
          res[i,-c(1:5)]          = mosekSol$sol$itr$barx[[1]]

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
                                certificate      = FALSE,
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
            mosekSolDF$certificate[i]      = checkEvaluation(coefMatrixSub[i,], degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
                                                             matrix(mosekSol$sol$itr$barx[[1]][1+1:varDim], nrow=1),
                                                             fabstol, freltol)
            mosekSolDF[i,-c(1:5)]          = mosekSol$sol$itr$barx[[1]]

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
                       certificate      = FALSE,
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
          res$certificate[i]      = checkEvaluation(coefMatrix[i,], degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
                                                    matrix(mosekSol$sol$itr$barx[[1]][1+1:varDim], nrow=1),
                                                    fabstol, freltol)
          res[i,-c(1:5)]          = mosekSol$sol$itr$barx[[1]]

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
                                certificate      = FALSE,
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
            mosekSolDF$certificate[i]      = checkEvaluation(coefMatrixSub[i,], degrees, mosekSol$sol$itr$pobjval, mosekSol$sol$itr$dobjval,
                                                             matrix(mosekSol$sol$itr$barx[[1]][1+1:varDim], nrow=1),
                                                             fabstol, freltol)
            mosekSolDF[i,-c(1:5)]          = mosekSol$sol$itr$barx[[1]]

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
              certificate      = res$certificate,
              varDim           = varDim,
              orderobj         = orderObj,
              orderMom         = orderMom,
              moment_vectors   = res[,-c(1:5)]))
}

extractSolution = function(sol, points=NULL, vabstol=1e-10, ereltol=1e-03) {

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
    nSol = sol$rank
  } else {
    nSol = points
  }

  # compute left-matrix of SVD of optimal moment matrix
  if(nSol == 1) {

    ##### if nSol = 1, return closed-form solution and quit #####
    # singularMatrix = with(svd(optimalMomentMatrix, nu = nSol, nv = 0), u * d[1])
    argMat = matrix(optimalMomentMatrix[2:(1+varDim),1], nrow=1, ncol=varDim)
    return(argMat)

  } else {

    singularMatrix = with(svd(optimalMomentMatrix, nu = nSol, nv = 0), u %*% diag(d[1:nSol]))

  }

  # set infinitesimal values to be zero
  singularMatrix[abs(singularMatrix) < vabstol] = 0

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

optquad = function(sense, coefs, degrees, constraints=NULL, opt=NULL) {
  
  # read options
  if(is.null(opt) == TRUE) {
    options = list()
  } else {
    options = opt
  }
  
  if(is.null(options$multithread) == TRUE) {
    options$multithread = FALSE
  }
  
  if(is.null(options$verbose) == TRUE) {
    options$verbose = 0
  }
  
  # read constraints
  if(is.null(constraints) == TRUE) {
    constr = list()
  } else {
    constr = constraints
  }
  
  # initialize Mosek model
  monomialVector = createMonomialVector(ncol(degrees), 2)
  model = createMosekQuadraticModelSkeleton_cpp(coefs, degrees, monomialVector$primes)
  
  # add labels to model skeleton
  model$sense = sense
  if(is.null(constr$A)  == FALSE) {
    model$A  = Matrix(constr$A) # depends on Matrix package
  } else {
    model$A = Matrix(0, nrow=1, ncol=ncol(degrees))
  }
  if(is.null(constr$bc) == FALSE) {
    model$bc = constr$bc
  } else {
    model$bc = matrix(0, nrow=2, ncol=1)
  }
  rownames(model$bc) = c("blc", "buc")
  if(is.null(constr$bx) == FALSE) {
    model$bx = constr$bx
  } else {
    model$bx = matrix(c(-Inf, Inf), nrow=2, ncol=ncol(degrees))
  }
  rownames(model$bx) = c("blx", "bux")
  if(options$multithread == TRUE) {
    model$iparam = list(MSK_IPAR_INTPNT_MULTI_THREAD = "MSK_ON") # turn on multithreading
  } else {
    model$iparam = list(MSK_IPAR_INTPNT_MULTI_THREAD = "MSK_OFF") # turn off multithreading
  }
  
  # solve quadratic program
  mosekSol = mosek(model, opts = list(verbose=options$verbose, soldetail=1))
  
  # return result
  if(is.nan(mosekSol$response$code) == TRUE) {
    sol = list(qpstatus = mosekSol$response$msg) # return error message if error occurred
  } else {
    sol = list(objective_primal = mosekSol$sol$itr$pobjval,
               objective_dual   = mosekSol$sol$itr$dobjval,
               qpstatus         = mosekSol$response$msg,
               solstatus        = mosekSol$sol$itr$solsta,
               solution         = mosekSol$sol$itr$xx)
  }
  return(sol)
}

optquadGurobi = function(sense, coefs, degrees, constraints=NULL, opt=NULL) {
  
  # read options
  if(is.null(opt) == TRUE) {
    options = list()
  } else {
    options = opt
  }
  
  if(is.null(options$Threads) == TRUE) {
    options$Threads = 0
  }
  
  # if(is.null(options$OutputFlag) == TRUE) {
  #   options$OutputFlag = 1 # can omit if OutputFlag = 1
  # }
  
  # read constraints
  if(is.null(constraints) == TRUE) {
    constr = list()
  } else {
    constr = constraints
  }
  
  # initialize Gurobi model
  monomialVector = createMonomialVector(ncol(degrees), 2)
  model = createGurobiQuadraticModelSkeleton_cpp(coefs, degrees, monomialVector$primes)
  
  # add labels to model skeleton
  model$modelsense = sense
  if(is.null(constr$A)  == FALSE) {
    model$A  = constr$A
  } else {
    model$A = Matrix(0, nrow=1, ncol=ncol(degrees)) # depends on Matrix package
  }
  if(is.null(constr$rhs) == FALSE) {
    model$rhs = constr$rhs
    model$sense = constr$sense
  } else {
    model$rhs = 0
    model$sense = "="
  }
  
  if(is.null(constr$lb) == FALSE) {
    model$lb = constr$lb
  }
  if(is.null(constr$ub) == FALSE) {
    model$ub = constr$ub
  }
  
  # solve quadratic program
  gurobiSol = gurobi(model, params = options)
  
  # return result
  if(gurobiSol$status != "OPTIMAL" & gurobiSol$status != "SUBOPTIMAL") {
    sol = list(qpstatus = gurobiSol$status) # return error message if error occurred
  } else {
    sol = list(objective = gurobiSol$objval,
               qpstatus  = gurobiSol$status,
               solution  = gurobiSol$x)
  }
  return(sol)
}
