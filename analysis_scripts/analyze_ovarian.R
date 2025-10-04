#!/usr/bin/env Rscript
source('scripts/00_setup_utils.R')
pkgs <- c('survival','survminer','ggplot2','dplyr','tidyr','gtsummary','broom','readr','officer','flextable','MASS')
load_or_install(pkgs)

outdir <- 'outputs/ovarian'
safe_mkdir(outdir)

dat <- readr::read_csv('data/ovarian.csv')
names(dat) <- tolower(names(dat))
if(!all(c('futime','fustat') %in% names(dat))) stop('ovarian.csv missing futime/fustat')

dat <- dat %>% mutate(rx = factor(rx), resid.ds = factor(resid.ds), ecog.ps = factor(ecog.ps))

tbl <- gtsummary::tbl_summary(dat, by = 'rx', include = c('age','resid.ds','ecog.ps')) %>% add_p() %>% add_n()
saveRDS(tbl, file = file.path(outdir,'demographics_gtsummary.rds'))
try(export_gtsummary_docx(tbl, file = file.path(outdir,'demographics.docx')), silent = TRUE)

survobj <- Surv(time = dat$futime, event = dat$fustat)
fitkm <- survfit(survobj ~ rx, data = dat)
km <- ggsurvplot(fitkm, data = dat, pval = TRUE, risk.table = TRUE)
ggsave(file.path(outdir,'km_rx.png'), plot = km$plot, width = 8, height = 6)
ggsave(file.path(outdir,'km_rx_risktable.png'), plot = km$table, width = 8, height = 2)

multiv <- coxph(survobj ~ age + resid.ds + rx + ecog.ps, data = dat)
saveRDS(multiv, file = file.path(outdir,'cox_multivariable.rds'))
tbl_cox <- gtsummary::tbl_regression(multiv, exponentiate = TRUE) %>% add_global_p()
saveRDS(tbl_cox, file = file.path(outdir,'cox_table_gtsummary.rds'))
try(export_gtsummary_docx(tbl_cox, file = file.path(outdir,'cox_multivariable.docx')), silent = TRUE)

zph <- cox.zph(multiv)
saveRDS(zph, file = file.path(outdir,'cox_zph.rds'))
png(file.path(outdir,'cox_zph_plot.png'), width=1000, height=800)
par(mfrow=c(2,2)); plot(zph); dev.off()

## interaction: rx x age
int_model <- update(multiv, . ~ . + rx:age)
saveRDS(int_model, file = file.path(outdir,'cox_interaction_rx_age.rds'))
saveRDS(anova(multiv, int_model), file = file.path(outdir,'cox_interaction_anova.rds'))

## stepwise selection
step_mod <- MASS::stepAIC(multiv, direction='both', trace=FALSE)
saveRDS(step_mod, file = file.path(outdir,'cox_stepAIC.rds'))
try(export_gtsummary_docx(gtsummary::tbl_regression(step_mod, exponentiate = TRUE), file = file.path(outdir,'cox_stepAIC.docx')), silent = TRUE)

## concordance
cindex <- survConcordance(survobj ~ predict(multiv), data = dat)
saveRDS(cindex, file = file.path(outdir,'cindex.rds'))

## time-dependent demo
if(max(dat$futime, na.rm=TRUE) > 1){
	cutpt <- median(dat$futime, na.rm=TRUE)
	td <- survSplit(dat, cut = cutpt, end='futime', start='tstart', event='fustat')
	td$age_post <- ifelse(td$tstart >= cutpt, td$age, NA)
	cox_td <- coxph(Surv(tstart, futime, fustat) ~ age_post + resid.ds + rx + ecog.ps, data = td)
	saveRDS(cox_td, file = file.path(outdir,'cox_time_dependent_example.rds'))
}

gg <- ggforest(multiv, data = dat)
ggsave(file.path(outdir,'cox_forest.png'), plot = gg, width = 8, height = 6)

cat('Ovarian analysis finished. Outputs in', outdir, '\n')
