
# Utility function for the optimal treatment rule estimation function 

# Explanation: Inflates the value of observed event times  by the inverse censoring distribution 
zeta_t <- function(t, data, w=c(0, 1, 0), hazard, S, T){
  debug = 0
  
  if(debug) print("zeta_t")
  
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(hazard$time)*(hazard$time<=t)))
  }
  events <- list()
  # print(events)
  # for each state 
  for(j in 1:length(S)){
    
    if(!(S[j] %in% T)){
      # not terminal state 
      # 
      events[[j]] <- (data$s1==S[j] & data$t1<=t & data$t2>t)
    } else {
      # terminal state
      # BOOL: Event happens in state j for patient i  
      events[[j]] <- (data$s2==S[j] & data$t2<=t)
    }
    # print(paste("S[j] %in% T)", S[j] %in% T))
    # print(events)
  }
  
  # patient-preference weighting of states
  event <- do.call(cbind, events)%*%w
  
  if(debug) print("zeta_t: weighted events")
  if(debug) print(events)
  
  # print("zeta_t_elements")
  Tt <- (1*!(data$s2 %in% T))*t + (1*(data$s2 %in% T))*mapply(min, data$t2, t) # T^t
  elements <- sapply(Tt, element)
  
  # print("hazard")
  G_t <- exp(-hazard$hazard[elements])
  res <- event/G_t
  
  if(debug) print("res")
  if(debug) print(res)
  
  return(res)
}
