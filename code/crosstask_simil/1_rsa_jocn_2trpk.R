source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))


rsarray <- readRDS(here("out", "corr-biased_unpre_jocn_2trpk.RDS"))
parcels <- dimnames(rsarray)$parcel
n.mods <- length(parcels) * length(subjs)


## build models:

lab <- mat2vec(rsarray[, , 1, 1])[, c(".row", ".col")]
lab <- lab %>% 
  separate(.row, c("row_task", "row_cond")) %>%
  separate(.col, c("col_task", "col_cond"))
is.bw.task <- lab$row_task != lab$col_task
is.wn.hi <- lab$row_cond == "hi" & lab$col_cond == "hi"
is.wn.lo <- lab$row_cond == "lo" & lab$col_cond == "lo"
is.bw.lo <- lab$row_cond != lab$col_cond
X <- cbind(hihi = is.wn.hi, lolo = is.wn.lo, hilo = is.bw.lo)
X <- X[is.bw.task, ]
X <- sweep(X, 2, colSums(X), "/")  ## holds the 



## fit models ----


## prep:

vecs <- apply(rsarray, c("parcel", "subj"), function(x) x[lower.tri(x)][is.bw.task])  ## vectorize, subset by btw task only
vecs <- as.data.table(vecs)


## regress:

rsafit <- function(x, y, nms = colnames(X)) {
  
  b <- crossprod(X, y)
  
  data.frame(b = b, term = nms)
  
}

stats.subjs <-
  vecs %>%
  group_by(parcel, subj) %>%
  nest %>%
  mutate(data = map(data, ~rsafit(X, .$value))) %>%
  unnest(cols = data)


## write ----

fwrite(stats.subjs, here("out", "subjs_jocn.csv"))

