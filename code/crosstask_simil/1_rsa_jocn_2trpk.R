source(here::here("code", "_packages.R"))
source(here("code", "_vars.R"))


rsarray <- readRDS(here("out", "corr-biased_unpre_jocn_2trpk.RDS"))
parcels <- dimnames(rsarray)$parcel


## build models ----
## create weight matrix X that, when dotted with a vectorized RSM, will take average of relevant cells of matrix.

## get condition labels of cells of RSM:

## NB: first two dims of rsarray are rows and columns of representational similarity matrix (RSM)
lab <- 
  expand.grid(dimnames(rsarray)[1:2]) %>%
  separate(.row, c("row_task", "row_cond")) %>%
  separate(.col, c("col_task", "col_cond"))  ## create data.frame of labels of cells of RSM
temp_rsm <- diag(dim(rsarray)[1])  ## template representational similarity matrix
lab <- lab[which(lower.tri(temp_rsm)), ]  ## get labels for unique off-diag elements

## booleans for features:

is.bw.task <- lab$row_task != lab$col_task
is.wn.hi <- lab$row_cond == "hi" & lab$col_cond == "hi"
is.wn.lo <- lab$row_cond == "lo" & lab$col_cond == "lo"
is.bw.lo <- lab$row_cond != lab$col_cond

## combine and scale

X <- cbind(hihi = is.wn.hi, lolo = is.wn.lo, hilo = is.bw.lo) & is.bw.task
X <- sweep(X, 2, colSums(X), "/")  ## scale so will take mean



## fit models ----

## prep:

vecs <- apply(rsarray, c("parcel", "subj"), function(x) x[lower.tri(x)])  ## vectorize
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

