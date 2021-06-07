## assumes running from ccplinux1
## input: subjs, task, session, glmname
## output: RDS files


source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))


roimeans_hilo <- array(
  NA,
  dim = c(
    length(subjs),
    length(parcellation$key),
    length(tasks),
    length(sessions)
  ),
  dimnames = list(subj = subjs, parcel = parcellation$key, task = tasks, session = sessions)
)


## loop over tasks, subjs, runs; write 3droistats to results


for (subj.i in seq_along(subjs)) {
  
  for (task.i in seq_along(tasks)) {
    
    for (session.i in seq_along(sessions)) {
      # subj.i = 1; parcel.i = 1; task.i = 1; session.i = 1
      
      fname.i <-
        file.path(
          dir.nil.dmcc2.afni, subjs[subj.i], "SURFACE_RESULTS",  tasks[task.i],
          paste0(sessions[session.i], "_", glms[task.i]),
          paste0(
            subjs[subj.i], "_timecourses_", sessions[session.i], "_", contrs[task.i],
            "_Coef_tents_Schaefer2018_400Parcels_7Networks_order_10K"
            )
          )
      
      if (!file.exists(paste0(fname.i, "_L.txt"))) next
      
      L <- colMeans(as.matrix(fread(paste0(fname.i, "_L.txt"))[target.trs[[task.i]], -(1:2)]))
      R <- colMeans(as.matrix(fread(paste0(fname.i, "_R.txt"))[target.trs[[task.i]], -(1:2)]))
      roimeans_hilo[subj.i, 001:200, task.i, session.i] <- L
      roimeans_hilo[subj.i, 201:400, task.i, session.i] <- R
      
    }

  }

}


saveRDS(roimeans_hilo, here("in", paste0("roistats_hilo_target_2trpk_unshifted.RDS")))
