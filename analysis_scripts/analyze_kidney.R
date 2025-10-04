#!/usr/bin/env Rscript
source('scripts/00_setup_utils.R')
pkgs <- c('survival','survminer','ggplot2','dplyr','tidyr','gtsummary','broom','readxl','officer','flextable','MASS')
load_or_install(pkgs)

outdir <- 'outputs/kidney'
safe_mkdir(outdir)

dat <- readxl::read_xlsx('data/kidney.xlsx')
names(dat) <- tolower(names(dat))

# normalize common names
if(!('time' %in% names(dat)) && 'futime' %in% names(dat)) names(dat)[names(dat)=='futime'] <- 'time'
if(!('status' %in% names(dat)) && 'fustat' %in% names(dat)) names(dat)[names(dat)=='fustat'] <- 'status'
required <- c('time','status')
if(!all(required %in% names(dat))) stop('kidney.xlsx must contain time and status (or futime/fustat)')

demovars <- intersect(c('age','sex','prior'), names(dat))
if(length(demovars)==0) demovars <- names(dat)[1:min(3,ncol(dat))]

tbl <- gtsummary::tbl_summary(dat, include = demovars) %>% add_n()
saveRDS(tbl, file = file.path(outdir,'demographics_gtsummary.rds'))
try(export_gtsummary_docx(tbl, file = file.path(outdir,'demographics.docx')), silent = TRUE)

survobj <- Surv(time = dat$time, event = dat$status)
fitkm <- survfit(survobj ~ 1, data = dat)
km <- ggsurvplot(fitkm, data = dat, risk.table = TRUE)
ggsave(file.path(outdir,'km_overall.png'), plot = km$plot, width = 8, height = 6)
ggsave(file.path(outdir,'km_overall_risktable.png'), plot = km$table, width = 8, height = 2)

multiv_formula <- as.formula(paste('survobj ~', paste(demovars, collapse = ' + ')))
multiv <- coxph(multiv_formula, data = dat)
saveRDS(multiv, file = file.path(outdir,'cox_multivariable.rds'))
tbl_cox <- gtsummary::tbl_regression(multiv, exponentiate = TRUE) %>% add_global_p()
saveRDS(tbl_cox, file = file.path(outdir,'cox_table_gtsummary.rds'))
try(export_gtsummary_docx(tbl_cox, file = file.path(outdir,'cox_multivariable.docx')), silent = TRUE)

zph <- cox.zph(multiv)
saveRDS(zph, file = file.path(outdir,'cox_zph.rds'))
png(file.path(outdir,'cox_zph_plot.png'), width=1000, height=800)
par(mfrow=c(2,2)); plot(zph); dev.off()

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
  # pick first demovar to show post-change
  varpost <- demovars[1]
  td[[paste0(varpost,'_post')]] <- ifelse(td$tstart >= cutpt, td[[varpost]], NA)
  form <- as.formula(paste('Surv(tstart, time, status) ~', paste0(paste0(varpost,'_post'), collapse = ' + ')))
  cox_td <- coxph(formula = form, data = td)
  saveRDS(cox_td, file = file.path(outdir,'cox_time_dependent_example.rds'))
}

gg <- ggforest(multiv, data = dat)
ggsave(file.path(outdir,'cox_forest.png'), plot = gg, width = 8, height = 6)

cat('Kidney analysis finished. Outputs in', outdir, '\n')
