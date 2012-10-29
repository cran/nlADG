selectInd = function(indexList, selectedIndices) {
  ind = c()

  for (v in selectedIndices) {
    ind = c(ind, indexList[[v]])
  }

  ind
}

assertDsubspace = function(formulas, D, ...) {
  for (v_ in colnames(D)) {
    for (p_ in pa(v_, D)) {
      designmatrix_v = as.matrix(model.matrix(formulas[[v_]]))
      expanded_designmatrix_v = cbind(designmatrix_v, as.matrix(model.matrix(formulas[[p_]])))
      
      if (qr(designmatrix_v, ...)$rank < qr(expanded_designmatrix_v, ...)$rank) {
        stop("The regression system does not fulfill the D-linear subspace characteristics.")
      }
    }
  }
}

assertNoObservations = function(formulas, D, ...) {
  n_ = dim(model.frame(formulas[[1]]))[1]
  
  for (v_ in colnames(D)) {
    rankZ = qr(as.matrix(model.matrix(formulas[[v_]])), ...)$rank
    dimvandparents = dim(as.matrix(lm(formulas[[v_]])$coefficients))[2]
    for (p_ in pa(v_, D)) {
      dimvandparents = dimvandparents + dim(as.matrix(lm(formulas[[p_]])$coefficients))[2]
    }
    
    if (n_ < rankZ + dimvandparents) {
      stop("Not enough observations are given for the maximum likelihood estimation to exist and be unique.")
    }
  }
}

assertPDI = function(M, D, I) {
  D_ = diag(1, ncol=dim(D)[2], nrow=dim(D)[1])
  dimnames(D_) = dimnames(D)

  nd_ = list()
  for(v_ in colnames(D)) {
    nd_[[v_]] = setdiff(colnames(D), c(v_, pa(v_, D)))
  }

  for(v_ in colnames(D)) {
    D_ = D_ %*% D

    for(v__ in colnames(D)) {
      nd_[[v__]] = setdiff(nd_[[v__]], colnames(D_)[round(D_[v__,])==1])
    }
  }

  for(v_ in colnames(D)) {
    nd_v = nd_[[v_]]

    if(length(pa(v_, D)) >= 1 && length(nd_v) >= 1) {
      lM = M[selectInd(I, v_), selectInd(I, nd_v)]
      rM = M[selectInd(I, v_), selectInd(I, pa(v_, D))] %*% solve(M[selectInd(I, pa(v_, D)),selectInd(I, pa(v_, D))]) %*% M[selectInd(I, pa(v_, D)), selectInd(I, nd_v)]

      diff = abs(lM - rM)
      diff = diff * (diff > sqrt(.Machine$double.eps))

      if(max(diff) > 0) {
        stop("The matrix is not in P(D,I)")
      }
    }
  }
}


nlADG = function(formulas, ...) UseMethod("nlADG")

nlADG.default = function(formulas, D, ...) {
  formulas = as.list(formulas)
  
  if(!hasArg(D) || !is.matrix(D)) {
    D = matrix(0, ncol=length(formulas), nrow=length(formulas), dimnames=list(names(formulas), names(formulas)))
  } else {
    D = as.matrix(D)
  }
    
  D = topSort(D)
  
  assertDsubspace(formulas, D, ...)
  assertNoObservations(formulas, D, ...)
  
  DParameterEst = estimateDParameters(formulas, D, ...)
  
  modelParameterEst = reconstructDistParameters(DParameterEst$beta_hat, DParameterEst$mu_hat, DParameterEst$lambda_hat, D, DParameterEst$I, ...)
  
  model = list()
  
  model$coefficients = DParameterEst$coefficients
  model$sigma = modelParameterEst$sigma_hat
  model$xi = modelParameterEst$xi_hat
  model$fitted.values = modelParameterEst$xi_hat
  model$residuals = DParameterEst$residuals
  model$MANOVAmodels = DParameterEst$models
  model$ADG = D
  model$I = DParameterEst$I
  model$call = match.call()
  
  class(model) = "nlADG"
  
  model
}

estimateDParameters = function(formulas, D, ...) {
  n_ = dim(model.frame(formulas[[1]]))[1]
  
  mu_hat = list()
  beta_hat = list()
  lambda_hat = list()
  fittedvalues = data.frame( dummy=1:n_ )
  coefficients = list()
  residuals = data.frame( dummy=1:n_ )
  models = list()
  
  I = list()
  Icounter = 0
  
  for(v_ in colnames(D)) {
    model_v = lm(formulas[[v_]])
    
    if (length(pa(v_, D)) >= 1) {
      no_predictors = dim(model.matrix(formulas[[v_]]))[2]
      
      model_frame_expanded_v = expand.model.frame(formulas[[v_]], pa(v_, D))
      model_v = lm(model_frame_expanded_v)
      
      no_predictors_expanded = dim(as.matrix(model_v$coef))[1]
      
      beta_hat[[v_]] = matrix( as.matrix(model_v$coef)[(no_predictors+1):no_predictors_expanded,] , nrow=no_predictors_expanded-no_predictors)
      
      clearedDataframe = model_frame_expanded_v
      for(p_ in pa(v_, D)) {
        if (is.matrix(clearedDataframe[[p_]])) {
          clearedDataframe[[p_]] = matrix(0, ncol=dim(as.matrix(model_frame_expanded_v[[p_]]))[2], nrow=dim(as.matrix(model_frame_expanded_v[[p_]]))[1])
        } else {
          clearedDataframe[[p_]] = rep(0, dim(as.matrix(model_frame_expanded_v[[p_]]))[1])
        }
      }
      mu_hat[[v_]] = as.matrix(predict(model_v, clearedDataframe))
      
      fittingDataframe = model_frame_expanded_v
      for (p_ in pa(v_, D)) {
        if (is.matrix(fittingDataframe[[p_]])) {
          fittingDataframe[[p_]] = fittedvalues[[p_]]
        } else {
          fittingDataframe[[p_]] = as.vector(fittedvalues[[p_]])
        }
      }
      fittedvalues[[v_]] = as.matrix(predict(model_v, fittingDataframe))
    } else {
      mu_hat[[v_]] = as.matrix(model_v$fitted.values)
      
      beta_hat[[v_]] = matrix(0, ncol=dim(mu_hat[[v_]])[2], nrow=0)
      
      fittedvalues[[v_]] = as.matrix(model_v$fitted.values)
    }
    
    original_values_v = as.matrix(model_v$residuals) + as.matrix(model_v$fitted.values)
    
    residuals[[v_]] = original_values_v - fittedvalues[[v_]]
    
    lambda_hat[[v_]] = t(as.matrix(model_v$residuals)) %*% original_values_v / dim(mu_hat[[v_]])[1]
    
    coefficients[[v_]] = as.matrix(model_v$coefficients)
    
    if (dim(mu_hat[[v_]])[2] >= 2) {
      rv_names = paste(v_, 1:dim(mu_hat[[v_]])[2], sep=".")
    } else {
      rv_names = v_
    }
    
    dimnames(beta_hat[[v_]]) = list(c(), rv_names)
    dimnames(mu_hat[[v_]]) = list(c(), rv_names)
    dimnames(lambda_hat[[v_]]) = list(rv_names, rv_names)
    dimnames(coefficients[[v_]]) = list(rownames(coefficients[[v_]]), rv_names)

    models[[v_]] = model_v
    
    I[[v_]] = Icounter + 1:dim(mu_hat[[v_]])[2]
    Icounter = Icounter + dim(mu_hat[[v_]])[2]
  }
  
  fittedvalues[["dummy"]] = NULL
  residuals[["dummy"]] = NULL

  
  list(beta_hat=beta_hat, mu_hat=mu_hat, lambda_hat=lambda_hat, residuals=residuals, coefficients=coefficients, I=I, models=models, fittedvalues=fittedvalues)
}

reconstructDistParameters = function(beta, mu, lambda, D, I, ...) {
  V = colnames(D)
  
  n_ = dim(mu[[1]])[1]
  
  dimI = length(selectInd(I, V))
  
  sigma_hat = matrix(0, nrow=dimI, ncol=dimI)
  xi_hat = data.frame( dummy=1:n_ )
  
  cumIndices = c()
  
  columnnames = vector(mode="character", length=dimI)
  
  for(v_ in V) {
    sigma_hat[selectInd(I, v_), selectInd(I, pa(v_, D))] =  t(beta[[v_]]) %*% sigma_hat[selectInd(I, pa(v_, D)), selectInd(I, pa(v_, D))]
    sigma_hat[selectInd(I, pa(v_, D)), selectInd(I, v_)] = t(sigma_hat[selectInd(I, v_), selectInd(I, pa(v_, D))])
    sigma_hat[selectInd(I, v_), selectInd(I, v_)] = lambda[[v_]] + t(beta[[v_]]) %*% sigma_hat[selectInd(I, pa(v_, D)), selectInd(I, v_)]
    sigma_hat[selectInd(I, v_), selectInd(I, setdiff(cumIndices, pa(v_, D)))] = t(beta[[v_]]) %*% sigma_hat[selectInd(I, pa(v_, D)),selectInd(I, setdiff(cumIndices, pa(v_, D)))]
    sigma_hat[selectInd(I, setdiff(cumIndices, pa(v_, D))), selectInd(I, v_)] = t(sigma_hat[selectInd(I, v_), selectInd(I, setdiff(cumIndices, pa(v_, D)))])
    
    cumIndices = c(cumIndices, v_)
    
    xi_hat[[v_]] = mu[[v_]]
    
    if (length(pa(v_, D)) >= 1) {
      xi_hat[[v_]] = mu[[v_]] + as.matrix(xi_hat[pa(v_, D)]) %*% beta[[v_]]
    }
    
    colnames(xi_hat[[v_]]) = NULL

    columnnames[selectInd(I,v_)] = colnames(mu[[v_]])
  }
  
  dimnames(sigma_hat) = list(columnnames, columnnames)

  xi_hat[["dummy"]] = NULL
  
  list(xi_hat=xi_hat, sigma_hat=sigma_hat)
}

predict.nlADG = function(object, newdata, ...) {
  newdata = as.data.frame(newdata)
  models = object$MANOVAmodels
  D = object$ADG
  V = colnames(D)
  n_ = dim(newdata)[1]
  
  
  prediction = data.frame( dummy=1:n_ )
  
  for (v_ in V) {
    data_v = newdata
    for (p_ in pa(v_, D)) {
      if (is.matrix(model.frame(models[[v_]])[[p_]])) {
        data_v[[p_]] = as.matrix(prediction[[p_]])
      } else {
        data_v[[p_]] = as.vector(prediction[[p_]])
      }
    }
    
    prediction[[v_]] = as.matrix(predict(models[[v_]], data_v))
  }

  prediction[["dummy"]] = NULL
  
  prediction
}

print.nlADG = function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMLE xi:\n")
  print(x$xi)
  cat("\nMLE sigma:\n")
  print(x$sigma)
}

summary.nlADG = function(object, ...) {
  print(object, ...)
}

print.summary.nlADG = function(x, ...) {
  print(x, ...)
}