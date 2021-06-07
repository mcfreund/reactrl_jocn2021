source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))



n.iter <- length(subjs55)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)


conditions <- combo_paste(tasks, c("hi", "lo"))


## wrangle betas ----
(time.start <- Sys.time())

betas <- 
  array(
    NA,
    dim = c(n.vert, length(conditions), length(subjs55)),
    dimnames = list(vertex = NULL, condition = conditions, subj = subjs55)
  )
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 2
  
  name.glm.i <- glminfo$name.glm[glm.i]
  name.task.i <- glminfo$task[glm.i]
  
  ## read betas:
  betas.i <- readRDS(
    here("in", paste0("betas_dmcc_2trpk_", name.task.i, "_", name.glm.i,  ".RDS"))
  )
  ## remove subj with missing data, get target TRs:
  betas.i <- betas.i[, , target.trs[[glm.i]], ]
  betas.i <- aperm(betas.i, c("vertex", "subj", "reg"))
  
  betas.i.hi <- rowMeans(betas.i[, , hi[[glm.i]], drop = FALSE], dim = 2) ## ave across TR, run, relevant regs
  betas.i.lo <- rowMeans(betas.i[, , lo[[glm.i]], drop = FALSE], dim = 2)
  betas.ii <- abind(betas.i.hi, betas.i.lo, along = 0)
  
  ## reshape to match dim(betas):
  
  betas.ii <- aperm(betas.ii, c(2, 1, 3))
  conditions.i <- combo_paste(name.task.i, c("hi", "lo"))  ## run varies slower
  betas[, conditions.i, ] <- betas.ii[, , subjs55]
  
  
}
rm(betas.i)
gc()






## est simil ----


simil <- 
  array(
    NA,
    dim = c(length(conditions), length(conditions), length(parcellation$key), length(subjs55)),
    dimnames = list(.row = conditions, .col = conditions, parcel = parcellation$key, subj = subjs55)
  )


for (subj.i in seq_along(subjs55))  {
  # subj.i = 1
  
  res <- enlist(parcellation$key)
  
  
  name.subj.i <- subjs55[subj.i]
  
  # if (name.subj.i %in% c("DMCC5820265")) next
  
  betas.subj.i <- betas[, , subj.i]
  
  
  for (parcel.i in seq_along(parcellation$key)) {
    # parcel.i = 20
    
    ## mask:
    
    is.parcel <- schaefer10k == parcel.i
    betas.subj.parcel.i <- betas.subj.i[is.parcel, ]
    
    has.signal.all.conds <- !rowSums(abs(betas.subj.parcel.i) < .Machine$double.eps) > 0
    if (sum(has.signal.all.conds) <= length(conditions)) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))
    
    B <- betas.subj.parcel.i[has.signal.all.conds, ]
    
    
    ## estimate similarity matrices:
    
    simil[, , parcel.i, subj.i] <- cor(B)
    
  }
  
  pb$tick()  ## progress bar
  
  
}



## save ----

saveRDS(simil, here("out", paste0("corr-biased_unpre_jocn_2trpk.RDS")))


(time.run <- Sys.time() - time.start)
