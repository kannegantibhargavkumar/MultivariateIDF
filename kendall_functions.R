### this function calculates the kendall distribution for seven copulas
## clayton, frank,gumbel, bb1,bb6,bb7,bb8, joe copula
kendall_dist=function(model, t){
  
  model_name=model[["familyname"]]
  
  if (model_name=="Clayton" |model_name=="Survival Clayton"){
    
    kendall.Clayton <- function(copula, t){
      par = copula[["par"]]
      
      kt <- rep(NA,length(t))
      kt <- t + t * (1 - t^par)/par
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.Clayton(model,t)
    
  }
  else if(model_name=="Frank"){
    kendall.Frank <- function(copula, t){
      par = copula[["par"]]
      
      kt <- rep(NA,length(t))
      kt <- t + log((1 - exp(-par))/(1 - exp(-par * t))) * (1 - exp(-par * t))/(par * exp(-par * t))
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.Frank(model,t)
    
  }
  else if(model_name=="Gumbel"){
    kendall.Gumbel <- function(copula, t){
      par = copula[["par"]]
      
      kt <- rep(NA,length(t))
      kt <- t - t * log(t)/(par)
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.Gumbel(model,t)
  }
  else if (model_name=="BB1"){
    kendall.BB1 <- function(copula, t){
      theta = copula[["par"]]
      delta = copula[["par2"]]
      
      kt <- rep(NA,length(t))
      kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.BB1(model,t)
    
  } else if (model_name=="BB6"){
    kendall.BB6 <- function(copula, t){
      theta = copula[["par"]]
      delta = copula[["par2"]]
      
      kt <- rep(NA,length(t))
      kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.BB6(model,t)
    
  } else if (model_name=="BB7"){
    kendall.BB7 <- function(copula, t){
      theta = copula[["par"]]
      delta = copula[["par2"]]
      
      kt <- rep(NA,length(t))
      kt <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) -  1)/
        ((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.BB7(model,t)
    
  } else if (model_name=="BB8"){
    kendall.BB8 <- function(copula, t){
      theta = copula[["par"]]
      delta = copula[["par2"]]
      
      kt <- rep(NA,length(t))
      kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/ (theta * delta)
      kt[t==1] <- 1
      kt[t==0] <- 0
      return(kt)  
    }
    kt_e=kendall.BB8(model,t)
    
  } else if (model_name=="Bijoe"){
    kendall.Joe <- function(copula, t) kdJoe(t, copula)
    kt_e=kendall.Joe(model,t)
    
  }
  
  return(kt_e)
}




