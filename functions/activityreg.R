f.activityreg <- function(act_df, sim_time, n_step){
  ## -----------------------------------------------------------------------------
  #'  Creates a dose regime  approximation function that can be called during the 
  #'  PBK simulations. 
  #'  Input: act_df; dataframe containing the following columns:
  #'            level (activity level; 'light', 'moderate', 'heavy')
  #'            n_day (number of activities per day)
  #'            start (start time [h] of the activity on any day)
  #'            dur (duration of activity [h])
  #'            int (if multiple events take place, interval [h] between the start of the activities. 
  #'                (e.g. if an activity start at 8h, with an interval of 4, then the second activity starts at 12h)  
  #'            t_stop (once this time [h] has passes within a simulation, that particular activity level should stop.)
  #'         sim_time (simulation time [h])   
  #'  Returns: activity (list with components time and level)
  #'  
  ## -----------------------------------------------------------------------------
  

  # Define timing
  times <- seq(0, sim_time, length=sim_time*n_step+1)
  
  # Define activity list (default is 'rest')
  levels <- rep('rest', length(times))
  
  for (i in 1:nrow(act_df)) {
    
    # Initial checks
    if (act_df$int[i]<=0 ) {
      if (act_df$n_day[i]>1) {
        stop("Multiple activities were provided without a proper time interval")}
      else{
        act_df$int[i] <- act_df$dur[i]+1}}
    if (act_df$t_stop[i] > sim_time) {
      warning("Activity stop time exceed the simulation time")
    }
    if (act_df$int[i]<=act_df$dur[i]) {
      warning("Interval time was provided that is smaller than the activity duration. Assumed continuous activity")
    }
    
    
    for (n in 1:act_df$n_day[i]) {
     
      levels[(times%%24 >= act_df$start[i]) & 
               (times%%24 < act_df$start[i] +(act_df$n_day[i]-1)*act_df$int[i]+act_df$dur[i]) &
               times < act_df$t_stop[i]] = act_df$level[i]
      }
  }
  return(levels)
}
