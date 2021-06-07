## assumes running from ccplinux1
## input: subjs, task, glmname
## output: RDS files


source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))


dirs <- expand.grid(subj = subjs55, task = tasks, session = "baseline", stringsAsFactors = FALSE)


read_betas_dmcc <- function(
  .subjs,
  .task,
  .glm,
  .dir
) {
  # glm.i = 1; .subjs = subjs55; .task = glminfo$task[glm.i]; 
  # .glm = glminfo$name.glm[glm.i]; .dir = dir.analysis
  
  ## initialize array
  
  pick.a.file <- 
    file.path(.dir, .subjs[1], "SURFACE_RESULTS",  .task, paste0(.glm), paste0("STATS_", .subjs[1], "_REML_L.func.gii"))
  labs <- afni("3dinfo", paste0("-label ", pick.a.file))
  labs <- unlist(strsplit(labs, "\\|"))
  is.reg <- !grepl("Full|block|Tstat|Fstat", labs)
  tab <- do.call(rbind, strsplit(gsub("_Coef", "", labs[is.reg]), "#"))
  trs <- as.numeric(unique(tab[, 2])) + 1
  regs <- unique(tab[, 1])
  
  n.vertex <- 10242
  n.tr <- length(trs)
  n.reg <- length(regs)
  n.subj <- length(.subjs)
  
  betas <- array(
    NA,
    dim = c(n.vertex*2, n.reg, n.tr, n.subj),
    dimnames = list(vertex = NULL, reg = regs, tr = NULL, subj = .subjs)
  )
  
  vertex.inds <- cbind(L = 1:n.vertex, R = (n.vertex + 1):(n.vertex * 2))
  
  for (subj.i in seq_along(.subjs)) {
    # subj.i = 41; hemi.i = "L"
    
    for (hemi.i in c("L", "R")) {
      # hemi.i = "R"
      
      inds <- vertex.inds[, hemi.i]
      
      fname <- file.path(
        .dir, .subjs[subj.i], "SURFACE_RESULTS",  .task, paste0(.glm),  
        paste0("STATS_", .subjs[subj.i], "_REML_", hemi.i, ".func.gii")
      )
      
      if (!file.exists(fname)) next
      
      B <- mikeutils::read_gifti2matrix(fname)[is.reg, ]
      
      is.ok.i <- isTRUE(all.equal(dim(B), c(n.reg * n.tr, n.vertex)))
      if (!is.ok.i) stop("mismatched beta array")
      
      
      for (reg.i in seq_len(n.reg)) {
        # reg.i = 1
        
        is.reg.i <- grepl(paste0("^", regs[reg.i], "#"), labs[is.reg])
        B.reg.i <- t(B[is.reg.i, ])
        
        is.ok.ii <- isTRUE(all.equal(dim(betas[inds, reg.i, , subj.i]), dim(B.reg.i)))
        if (!is.ok.ii) stop("mismatched regressor array")
        
        betas[inds, reg.i, , subj.i] <- B.reg.i
        
      }
      
    }
    
  }
  
  betas
  
}


for (glm.i in seq_len(nrow(glminfo))) {
  
  betas.i <- read_betas_dmcc(subjs55, glminfo[glm.i]$task, glminfo[glm.i]$name.glm, dir.analysis)
  saveRDS(
    betas.i, 
    here("in", paste0("betas_dmcc_2trpk_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
    )
  
}
