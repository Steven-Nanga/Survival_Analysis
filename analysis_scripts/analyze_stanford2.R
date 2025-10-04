#!/usr/bin/env Rscript
source('scripts/00_setup_utils.R')
pkgs <- c('survival','survminer','ggplot2','dplyr','tidyr','gtsummary','broom','readr','officer','flextable','MASS')
load_or_install(pkgs)

outdir <- 'outputs/stanford2'
safe_mkdir(outdir)

dat <- readr::read_csv('data/stanford2.csv')
names(dat) <- tolower(names(dat))

required <- c('time','status')
if(!all(required %in% names(dat))) stop('stanford2.csv must contain time and status columns')

dat <- dat %>% mutate(t5 = as.numeric(t5))

tbl <- gtsummary::tbl_summary(dat, include = c('age','t5')) %>% add_n()
saveRDS(tbl, file = file.path(outdir,'demographics_gtsummary.rds'))
try(export_gtsummary_docx(tbl, file = file.path(outdir,'demographics.docx')), silent = TRUE)

survobj <- Surv(time = dat$time, event = dat$status)
dat$agecat <- factor(dat$age > 50, labels = c('<=50','>50'))
fitkm <- survfit(survobj ~ agecat, data = dat)
km <- ggsurvplot(fitkm, data = dat, pval = TRUE, risk.table = TRUE)
ggsave(file.path(outdir,'km_agecat.png'), plot = km$plot, width = 8, height = 6)
ggsave(file.path(outdir,'km_agecat_risktable.png'), plot = km$table, width = 8, height = 2)

multiv <- coxph(survobj ~ age + t5, data = dat)
saveRDS(multiv, file = file.path(outdir,'cox_multivariable.rds'))
tbl_cox <- gtsummary::tbl_regression(multiv, exponentiate = TRUE) %>% add_global_p()
saveRDS(tbl_cox, file = file.path(outdir,'cox_table_gtsummary.rds'))
try(export_gtsummary_docx(tbl_cox, file = file.path(outdir,'cox_multivariable.docx')), silent = TRUE)

zph <- cox.zph(multiv)
saveRDS(zph, file = file.path(outdir,'cox_zph.rds'))
png(file.path(outdir,'cox_zph_plot.png'), width=1000, height=800)
par(mfrow=c(2,2)); plot(zph); dev.off()

## interaction age:t5
int_model <- update(multiv, . ~ . + age:t5)
saveRDS(int_model, file = file.path(outdir,'cox_interaction_age_t5.rds'))
saveRDS(anova(multiv, int_model), file = file.path(outdir,'cox_interaction_anova.rds'))

## stepwise
step_mod <- MASS::stepAIC(multiv, direction='both', trace=FALSE)
saveRDS(step_mod, file = file.path(outdir,'cox_stepAIC.rds'))
try(export_gtsummary_docx(gtsummary::tbl_regression(step_mod, exponentiate = TRUE), file = file.path(outdir,'cox_stepAIC.docx')), silent = TRUE)

## concordance
cindex <- survConcordance(survobj ~ predict(multiv), data = dat)
saveRDS(cindex, file = file.path(outdir,'cindex.rds'))

## time-dependent demo
if(max(dat$time, na.rm=TRUE) > 1){
	cutpt <- median(dat$time, na.rm=TRUE)
	td <- survSplit(dat, cut = cutpt, end='time', start='tstart', event='status')
	td$t5_post <- ifelse(td$tstart >= cutpt, td$t5, NA)
	cox_td <- coxph(Surv(tstart, time, status) ~ t5_post + age, data = td)
	saveRDS(cox_td, file = file.path(outdir,'cox_time_dependent_example.rds'))
}

gg <- ggforest(multiv, data = dat)
ggsave(file.path(outdir,'cox_forest.png'), plot = gg, width = 8, height = 6)

cat('stanford2 analysis finished. Outputs in', outdir, '\n')
