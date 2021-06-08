# reactrl_jocn2021

contains analyses data, code, and results for two supplementary analyses: **crosstask_simil** and **profile_analysis**

## crosstask_simil
A cross-task pattern similarity analysis:

<img src="https://github.com/mcfreund/reactrl_jocn2021/blob/main/out/fig_crosstasksimil_conditional.jpg?raw=true" alt="fig_crosstasksimil_conditional" width="500"/>

| script name (./code/crosstask_simil) | reads | writes |
| ----------- | -- | --- |
| 0_est_simil_jocn_2trpk.R      | beta coefs from fMRI 1st-level GLM: in/betas*.RDS (not in remote repo; too large) | pattern similarity matrices: out/corr*.RDS |
| 1_rsa_jocn_2trpk.R            | pattern similarity matrices: out/corr*.RDS              | subject-level similarity statistics: out/subjs_jocn.csv |
| 2_rsa_jocn_2trpk.rmd          | subject-level similarity statistics: out/subjs_jocn.csv | analysis report: code/crosstask_simil/\*.html; figure: out/fig_crosstasksimil_conditional.tiff |

* in/betas*.RDS wrangled from AFNI output files, which are not in remote repo, via code/\_wrangle/save_betas_dmcc_2trpk_surface.R


## profile analysis
A multivariate analysis of differences between sessions in profiles of activation:

<img src="https://github.com/mcfreund/reactrl_jocn2021/blob/main/out/fig_multivariate_triangle.jpg?raw=true" alt="fig_multivariate_triangle" width="500"/>

| script name (./code/profile_analysis) | reads | writes |
| ----------- | -- | --- |
| distance_jocn.rmd      | ROI-wise univariate contrasts on beta coefs from fMRI 1st-level GLM: in/roistats*.RDS | analysis report: code/profile_analysis/\*.html; figure: out/fig_multivariate_triangle.tiff |
* in/roistats*.RDS wrangled from AFNI output files via code/\_wrangle/save_roistats.R
