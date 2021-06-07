## functions ----



get.network <- function(x) {
  gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
}


cowplot.title <- function(x) {
  ggdraw() +
    draw_label(
      x,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
}


symmat4ggplot <- function(R, var.names = c("v1", "v2"), val.name = "value") {
  
  ## make factors for row and column labels
  dn <- dimnames(R)
  if (is.null(dn)) {
    dn <- setNames(list(paste0("cell_", 1:nrow(R)), paste0("cell_", 1:ncol(R))), var.names)
  } else {
    names(dn) <- var.names  
  }
  
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  labels[[2]] <- factor(labels[[2]], levels = rev(levels(labels[[2]])))
  
  r <- c(R)
  
  cbind(labels, setNames(as.data.frame(c(R)), val.name))
  
}


matplot <- function(x) {
  
  ggplot(symmat4ggplot(x), aes(v1, v2, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno") +
    theme_minimal() +
    theme(
      axis.text = element_blank(), axis.title = element_blank(), legend.position = "none",
      panel.border = element_blank(), panel.grid = element_blank()
    )
  
}



enlist <- function(x) setNames(vector("list", length(x)), x)


combo_paste <- function(a, b, ..., sep = "_") apply(expand.grid(a, b, ...), 1, paste0, collapse = sep)


read_schaefer <- function(path.atlas) {
  
  if (missing(path.atlas)) {
    nodename <- Sys.info()["nodename"]
    if (nodename == "ccplinux1") {
      path.atlas <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
    } else if (nodename == "CCP-FREUND") {  ## mike freund's (i.e., ccp's) thinkpad
      path.atlas <- "C:/local/atlases"
    }
  }
  
  fname.atlas <- file.path(
    path.atlas,
    "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti", "Schaefer2018_400Parcels_7Networks_order_info.txt"
  )
  if (file.exists(fname.atlas)) {
    fin <- file(fname.atlas, 'rt')
    tmp <- readLines(fin);
    close(fin); unlink(fin);
    if (length(tmp) != 800) { stop("not expected Schaefer key."); }
    tmp <- tmp[seq(from=1, to=800, by=2)];   # every-other entry is a label
    atlas.key <- gsub("7Networks_", "", tmp);
  }
  atlas <- cifti::read_cifti(
    file.path(
      path.atlas,
      "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti",
      "Schaefer2018_400Parcels_7Networks_order.dscalar.nii"
    ),
    drop_data = TRUE, trans_data = FALSE
  )
  atlas <- matrix(atlas$data)
  
  list(atlas = atlas, key = atlas.key)
  
}




contrast_matrix <- function(n, condition.names) {
  # n <- 10
  # condition.names <- letters[1:n]
  
  if (n < 2) stop("you need more than 1 condition you dummy")
  
  W <- matrix(0, nrow = n^2, ncol = n)
  
  if (missing(condition.names)) {
    dimnames(W) <- list(contrast = NULL, condition = NULL)
  } else {
    dimnames(W) <- list(
      contrast = paste0(rep(condition.names, each = n), "_", rep(condition.names, n)),
      condition = condition.names
    )
  }
  
  for (condition.i in seq_len(n)) {
    # condition.i = 1
    
    row.beg <- (condition.i - 1) * n + 1
    row.end <- (condition.i - 1) * n + n
    W.i <- W[row.beg:row.end, ]  ## square matrix; the contrasts that define a column of the similarity matrix
    
    W.i[, condition.i] <- 1  ## the condition to which all others are contrasted
    diag(W.i) <- diag(W.i) - 1  ## all others
    
    W[row.beg:row.end, ] <- W.i
    
  }
  
  W
  
}

mat2vec <- function(m, full.matrix = FALSE, varnames = c(".row", ".col"), ...) {
  
  if (any(is.na(m))) stop("matrix contains NA values.")
  if (!is.array(m)) stop("m is not array.")
  if (!full.matrix) m[upper.tri(m, diag = TRUE)] <- NA
  
  reshape2::melt(m, as.is = TRUE, na.rm = TRUE, varnames = varnames, ...)
  
}


## constants ----

nodename <- Sys.info()["nodename"]

dir.nil.dmcc2.afni <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/"

nodename <- Sys.info()["nodename"]

if (nodename == "ccplinux1") {
  
  dir.atlas <- "/data/nil-external/ccp/freund/atlases"
  dir.schaefer <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
  
} else if (nodename == "CCP-FREUND") {
  ## mike freund's (i.e., ccp's) thinkpad
  ## reliant on box drive
  ## assumes box drive location at ./Users/mcf/Box
  
  dir.atlas <- "C:/local/atlases"
  dir.schaefer <- dir.atlas
  
} else if (nodename == "PUTER") {
  
  dir.atlas <- "C:/Users/mcf/Documents/atlases"
  dir.schaefer <- "C:/Users/mcf/Documents/atlases/ATLASES/"
  
}


if (nodename %in% c("ccplinux1", "PUTER")) {
  
  schaefer10k <-
    c(
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
      gifti::read_gifti(file.path(dir.schaefer, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
    )
  
}

dir.analysis <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS"

n.cores <- detectCores()

tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
glms <- c(
  "Cues_EVENTS_censored", 
  "CongruencyIncentive_EVENTS_censored",
  "ListLength_EVENTS_censored",
  "Congruency_EVENTS_censored"
)
glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_Cues_EVENTS_censored",
    "baseline_CongruencyIncentive_EVENTS_censored",
    "baseline_ListLength_EVENTS_censored",
    "baseline_Congruency_EVENTS_censored"
  )
)
glminfo <- data.table::as.data.table(glminfo)

hi <- list(
  Axcpt = "BX",
  Cuedts = c("InConInc", "InConNoInc"),
  Stern = "LL5RN",
  Stroop = c("biasInCon", "PC50InCon")
)
lo <- list(
  Axcpt = "BY",
  Cuedts = c("ConInc", "ConNoInc"),
  Stern = "LL5NN",
  Stroop = c("biasCon", "PC50Con")
)




sessions <- c("baseline", "proactive", "reactive")
sessions.short <- c("Bas", "Pro", "Rea")




n.vert <- 20484


n.trs <- c(
  Axcpt_baseline   = 1220,
  Axcpt_proactive  = 1220,
  Axcpt_reactive   = 1220,
  Cuedts_baseline  = 1300,
  Cuedts_proactive = 1300,
  Cuedts_reactive  = 1300,
  Stern_baseline   = 1200,
  Stern_proactive  = 1200,
  Stern_reactive   = 1200,
  Stroop_baseline  = 1080,
  Stroop_proactive = 1080,
  Stroop_reactive  = 1180
)

target.trs <- list(
  Axcpt = 4,
  Cuedts = 4,
  Stern = 6,
  Stroop = 2
)


contrs <- c(
  "HI_LO_conf",
  "InCon_Con",
  "RN_NN_LL5",
  "InCon_Con_bias"
)


## read atlas and atlas key ----

parcellation <- read_schaefer(path.atlas = dir.schaefer)


dmcc35 <- c(
  77, 78,  86,  87,  88, 90, 91,  99, 101, 103, 105, 110, 127, 130, 132, 139, 140, 141, 148, 172, 175, 185, 189, 290, 
  306, 311, 314, 337, 340, 346, 347, 349, 350, 353, 361
  )

subjs <- c(
  "107321", "115825", "123117", "130114", "130518", "132017", "135730","138837", "141422", "150423", "155938", "158136", "160830", "161832", 
  "165032", "171330", "173738", "178243", "178647", "178950", "182436", "197449", "203418", "204319", "233326", "250427", "300618", "317332", 
  "346945", "352738", "393550", "448347", "580650", "594156", "601127", "672756", "729254", "765864", "814649", "843151", "849971", "873968", 
  "877168", "DMCC1328342", "DMCC1596165", "DMCC1624043", "DMCC1971064", "DMCC2442951", "DMCC2609759", "DMCC2803654", "DMCC2834766", 
  "DMCC3062542", "DMCC3204738", "DMCC3963378", "DMCC4191255", "DMCC4260551", "562345", "DMCC5009144", "DMCC5195268", "DMCC5775387",
  "DMCC6371570", "DMCC6418065", "DMCC6484785", "DMCC6627478", "DMCC6661074", "DMCC6671683", "DMCC6705371", "DMCC6721369", "DMCC6904377",
  "DMCC6960387", "DMCC7297690", "DMCC7921988", "DMCC8033964", "DMCC8050964", "DMCC8078683", "DMCC8214059", "DMCC8260571", 
  "DMCC9441378", "DMCC9478705", "DMCC9953810"
  )

subjs55 <- c(
  "102008", "107321", "115825", "123117", "130114", "130518", "132017", "141422", "150423", "158136", "161832", "165032", 
"173738", "176845", "178950", "182436", "203418", "204319", "233326", "300618", "346945", "393550", "432332", "448347", 
"580650", "594156", "601127", "672756", "765864", "849971", "877168", "DMCC0472647", "DMCC1328342", "DMCC1624043", 
"DMCC3963378", "DMCC4854984", "DMCC5065053", "DMCC5195268", "DMCC5244053", "DMCC5775387", "DMCC5820265", "DMCC6418065", 
"DMCC6627478", "DMCC6671683", "DMCC6705371", "DMCC7297690", "DMCC7921988", "DMCC8050964", "DMCC8078683", "DMCC8214059", 
"DMCC8260571", "DMCC9441378", "DMCC9478705", "DMCC9850294", "DMCC9953810"
)


colors.session <- setNames(qualitative_hcl(n = 3), sessions)



## settings ----

knitr::opts_chunk$set(
  cache = TRUE, warning = FALSE, messages = FALSE,
  fig.align = 'center',
  fig.width = 11.5, fig.fullwidth = TRUE
)

set.seed(0)

theme_set(theme_minimal(base_size = 14))


