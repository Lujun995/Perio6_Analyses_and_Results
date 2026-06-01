library(tidyverse)
library(ordinalNet)

get.adjacency <-  function(object, whichLambda = NULL){
  # get the adjacency matrix of the periodontal measurement dependency network
  if(object$call$definite_model){
    vars <- names(object$`data=1`$`step=1`)
    adjacency <- matrix(0, nrow = length(vars), ncol = length(vars),
                        dimnames = list(vars, vars))
    for(y_name in vars){
      coefs <- ordinalNet:::coef.ordinalNet(object$`data=1`$`step=1`[[y_name]],
                                            whichLambda = whichLambda)
      intercept_index <- grepl(x = names(coefs), pattern = "\\(Intercept\\)")
      coefs_var <- coefs[!intercept_index]
      adjacency[y_name, names(coefs_var)] <- coefs_var
    }
    return(adjacency)
  }
}