#!/usr/bin/env Rscript
source('scripts/00_setup_utils.R')
pkgs <- c('survival','survminer','ggplot2','dplyr','tidyr','gtsummary','broom','readr','officer','flextable','MASS')
load_or_install(pkgs)

outdir <- 'outputs/lung'
safe_mkdir(outdir)

dat <- readr::read_csv('data/lung.csv')
names(dat) <- tolower(names(dat))
required <- c('time','status')
if(!all(required %in% names(dat))) stop('lung.csv must contain time and status columns')

dat <- dat %>% mutate(sex = factor(sex))

## demographics
tbl <- gtsummary::tbl_summary(dat, by = 'sex', include = c('age','ph.ecog','ph.karno')) %>% add_p() %>% add_n()
saveRDS(tbl, file = file.path(outdir,'demographics_gtsummary.rds'))
try(export_gtsummary_docx(tbl, file = file.path(outdir,'demographics.docx')), silent = TRUE)

## KM
survobj <- Surv(time = dat$time, event = dat$status)
fitkm <- survfit(survobj ~ sex, data = dat)
km <- ggsurvplot(fitkm, data = dat, pval = TRUE, risk.table = TRUE, legend.title = 'Sex')
ggsave(file.path(outdir,'km_sex.png'), plot = km$plot, width = 8, height = 6)
ggsave(file.path(outdir,'km_sex_risktable.png'), plot = km$table, width = 8, height = 2)

## Cox
vars <- c('age','sex','ph.ecog','ph.karno')
univ <- lapply(vars, function(v) coxph(as.formula(paste0('survobj ~ ', v)), data = dat))
names(univ) <- vars
saveRDS(univ, file = file.path(outdir,'cox_univ_models.rds'))

multiv <- coxph(survobj ~ age + sex + ph.ecog + ph.karno, data = dat)
saveRDS(multiv, file = file.path(outdir,'cox_multivariable.rds'))
tbl_cox <- gtsummary::tbl_regression(multiv, exponentiate = TRUE) %>% add_global_p()
saveRDS(tbl_cox, file = file.path(outdir,'cox_table_gtsummary.rds'))
try(export_gtsummary_docx(tbl_cox, file = file.path(outdir,'cox_multivariable.docx')), silent = TRUE)

## PH test
zph <- cox.zph(multiv)
saveRDS(zph, file = file.path(outdir,'cox_zph.rds'))
png(file.path(outdir,'cox_zph_plot.png'), width=1000, height=800)
par(mfrow=c(2,2)); plot(zph); dev.off()

## interaction test (ph.ecog x age)
int_model <- update(multiv, . ~ . + ph.ecog:age)
saveRDS(int_model, file = file.path(outdir,'cox_interaction_ph.ecog_age.rds'))
saveRDS(anova(multiv, int_model), file = file.path(outdir,'cox_interaction_anova.rds'))

## stepwise AIC
step_mod <- MASS::stepAIC(multiv, direction='both', trace = FALSE)
saveRDS(step_mod, file = file.path(outdir,'cox_stepAIC.rds'))
try(export_gtsummary_docx(gtsummary::tbl_regression(step_mod, exponentiate = TRUE), file = file.path(outdir,'cox_stepAIC.docx')), silent = TRUE)

## Concordance
cindex <- survConcordance(survobj ~ predict(multiv), data = dat)
saveRDS(cindex, file = file.path(outdir,'cindex.rds'))

## time-dependent demo (split at median time)
if(max(dat$time, na.rm=TRUE) > 1){
	cutpt <- median(dat$time, na.rm=TRUE)
	td <- survSplit(dat, cut = cutpt, end='time', start='tstart', event='status')
	td$age_post <- ifelse(td$tstart >= cutpt, td$age, NA)
	cox_td <- coxph(Surv(tstart, time, status) ~ age_post + sex + ph.ecog + ph.karno, data = td)
	saveRDS(cox_td, file = file.path(outdir,'cox_time_dependent_example.rds'))
}

## forest plot
gg <- ggforest(multiv, data = dat)
ggsave(file.path(outdir,'cox_forest.png'), plot = gg, width = 8, height = 6)

cat('Lung analysis finished. Outputs in', outdir, '\n')
