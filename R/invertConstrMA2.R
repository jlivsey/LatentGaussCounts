# return abs of MA roots
invertConstrMA2 = function(theta){
  D = theta[1]^2-4*theta[2]
  if (D >=0){
    x1=(-theta[1]+sqrt(D))/(2*theta[2])
    x2=(-theta[1]-sqrt(D))/(2*theta[2])
    if(abs(x1)<1 | abs(x2)<1){
      return(1)
    }else{
      return(0)
    }
  }
  else{
    x = sqrt((-theta[1]/2*theta[2])^2 + (sqrt(abs(D))/(2*theta[2]))^2)
    if(abs(x)<1){
      return(1)
    }else{
      return(0)
    }
  }
}
