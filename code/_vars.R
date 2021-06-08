## functions ----

get.network <- function(x) gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)

enlist <- function(x) setNames(vector("list", length(x)), x)

combo_paste <- function(a, b, ..., sep = "_") apply(expand.grid(a, b, ...), 1, paste0, collapse = sep)



## constants ----


## schaefer atlas information

schaefer10k <- c(
  ## downloaded from: https://osf.io/7wk8v/
  gifti::read_gifti(here::here("in", "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
  ## downloaded from: https://osf.io/7tms8/
  gifti::read_gifti(here::here("in", "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
)
## requires internet connection:
schaeferkey <- read.csv(
  "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_400Parcels_7Networks_order_info.txt",
  header = FALSE
)[[1]]
schaeferkey <- schaeferkey[seq(from = 1, to = 800, by = 2)]  # every-other entry is a label
schaeferkey <- gsub("7Networks_", "", schaeferkey)


## strings

sessions <- c("baseline", "proactive", "reactive")
sessions.short <- c("Bas", "Pro", "Rea")
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
contrs <- c(
  "HI_LO_conf",
  "InCon_Con",
  "RN_NN_LL5",
  "InCon_Con_bias"
)

subjs80 <- c(
  "107321", "115825", "123117", "130114", "130518", "132017", "135730","138837", "141422", "150423", "155938", "158136",
  "160830", "161832", "165032", "171330", "173738", "178243", "178647", "178950", "182436", "197449", "203418", "204319",
  "233326", "250427", "300618", "317332", "346945", "352738", "393550", "448347", "580650", "594156", "601127", "672756", 
  "729254", "765864", "814649", "843151", "849971", "873968", "877168", "DMCC1328342", "DMCC1596165", "DMCC1624043", 
  "DMCC1971064", "DMCC2442951", "DMCC2609759", "DMCC2803654", "DMCC2834766", "DMCC3062542", "DMCC3204738", "DMCC3963378",
  "DMCC4191255", "DMCC4260551", "562345", "DMCC5009144", "DMCC5195268", "DMCC5775387", "DMCC6371570", "DMCC6418065",
  "DMCC6484785", "DMCC6627478", "DMCC6661074", "DMCC6671683", "DMCC6705371", "DMCC6721369", "DMCC6904377", "DMCC6960387",
  "DMCC7297690", "DMCC7921988", "DMCC8033964", "DMCC8050964", "DMCC8078683", "DMCC8214059", "DMCC8260571", "DMCC9441378",
  "DMCC9478705", "DMCC9953810"
)

subjs55 <- c(
  "102008", "107321", "115825", "123117", "130114", "130518", "132017", "141422", "150423", "158136", "161832", "165032", 
  "173738", "176845", "178950", "182436", "203418", "204319", "233326", "300618", "346945", "393550", "432332", "448347", 
  "580650", "594156", "601127", "672756", "765864", "849971", "877168", "DMCC0472647", "DMCC1328342", "DMCC1624043", 
  "DMCC3963378", "DMCC4854984", "DMCC5065053", "DMCC5195268", "DMCC5244053", "DMCC5775387", "DMCC5820265", "DMCC6418065", 
  "DMCC6627478", "DMCC6671683", "DMCC6705371", "DMCC7297690", "DMCC7921988", "DMCC8050964", "DMCC8078683", "DMCC8214059", 
  "DMCC8260571", "DMCC9441378", "DMCC9478705", "DMCC9850294", "DMCC9953810"
)


## paths to time-series data on server

if (Sys.info()["nodename"] == "ccplinux1") {
  
  dir.nil.dmcc2.afni <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/"
  dir.analysis <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS"
  
}


## nums

n.cores <- parallel::detectCores()

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

dmcc35 <- c(
  77, 78,  86,  87,  88, 90, 91,  99, 101, 103, 105, 110, 127, 130, 132, 139, 140, 141, 148, 172, 175, 185, 189, 290, 
  306, 311, 314, 337, 340, 346, 347, 349, 350, 353, 361
  )


## palettes

colors.session <- setNames(colorspace::qualitative_hcl(n = 3), sessions)




## settings ----


knitr::opts_chunk$set(
  cache = TRUE, warning = FALSE, messages = FALSE, fig.align = 'center', fig.width = 11.5, fig.fullwidth = TRUE
)

set.seed(0)

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 14))
