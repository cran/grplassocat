func_int <-
function(a){
  if(a[length(a)]!=0){a[1:(length(a)-1)]=-a[length(a)]} 
  return(a)}
