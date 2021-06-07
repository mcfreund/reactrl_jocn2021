## input: subjs, task, glmname
## output: RDS files


source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))


dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)

# subjs[subjs == "DMCC5820265"] <- "DMCC5820265_old-delete when rerun"

for (glm.i in seq_len(nrow(glminfo))) {
  
  betas.i <- read_betas_dmcc(subjs, glminfo[glm.i]$task, glminfo[glm.i]$name.glm, dir.analysis)
  saveRDS(
    betas.i, 
    here("out", paste0("betas_dmcc_2trpk_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
    )
  
}
