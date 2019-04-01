seqWrapper <-
function(lb, ub, by=1) {
  s = c()
  if(!ub < lb) s = seq(lb,ub, by=by)
  return(s)
}
