#!/usr/bin/env Rscript
source('scripts/00_setup_utils.R')
pkgs <- c('survival','survminer','ggplot2','dplyr','tidyr','gtsummary','broom','readr','officer','flextable')
load_or_install(pkgs)

outdir <- 'outputs/colon'
safe_mkdir(outdir)

dat <- readr::read_csv('data/colon.csv')
names(dat) <- tolower(names(dat))

## --- Data prep and checks -------------------------------------------------
required <- c('time','status')
if(!all(required %in% names(dat))) stop('colon.csv must contain time and status columns')
dat <- dat %>% mutate(
	sex = factor(sex, levels = c(0,1), labels = c('female','male')),
	rx = factor(rx),
	differ = factor(differ)
)

## --- Demographics / baseline table (gtsummary) ----------------------------
tbl <- gtsummary::tbl_summary(
	dat,
	by = 'rx',
	include = c('age','sex','nodes','differ'),
	statistic = list(all_continuous() ~ '{mean} ({sd})', all_categorical() ~ '{n} ({p}%)')
)
tbl_t <- tbl %>% add_p() %>% add_n()
print(tbl_t)
saveRDS(tbl_t, file = file.path(outdir,'demographics_gtsummary.rds'))
try(export_gtsummary_docx(tbl_t, file = file.path(outdir,'demographics.docx')),
		silent = TRUE)

## --- Kaplan-Meier curves -------------------------------------------------
survobj <- Surv(time = dat$time, event = dat$status)
fit_km_rx <- survfit(survobj ~ rx, data = dat)
kmplot <- ggsurvplot(
	fit_km_rx, data = dat, pval = TRUE, risk.table = TRUE,
	ggtheme = theme_minimal(), palette = 'Dark2',
	legend.title = 'Treatment', legend.labs = levels(dat$rx)
)
ggsave(filename = file.path(outdir,'km_rx.png'), plot = kmplot$plot, width = 8, height = 6)
ggsave(filename = file.path(outdir,'km_rx_risktable.png'), plot = kmplot$table, width = 8, height = 2)

## --- Cox proportional hazards models -------------------------------------
# Univariate screening
vars <- c('age','sex','rx','nodes','differ')
univ_models <- lapply(vars, function(v){
	f <- as.formula(paste0('survobj ~ ', v))
	coxph(f, data = dat)
})
names(univ_models) <- vars
univ_coefs <- lapply(univ_models, broom::tidy)
saveRDS(univ_models, file = file.path(outdir,'cox_univ_models.rds'))

# Multivariable model (pre-specified)
multiv_formula <- as.formula('survobj ~ age + sex + rx + nodes + differ')
cox_multi <- coxph(multiv_formula, data = dat)
saveRDS(cox_multi, file = file.path(outdir,'cox_multivariable.rds'))

# Model summary and nice table
tbl_cox <- gtsummary::tbl_regression(cox_multi, exponentiate = TRUE)
tbl_cox <- tbl_cox %>% add_global_p()
saveRDS(tbl_cox, file = file.path(outdir,'cox_table_gtsummary.rds'))
try(export_gtsummary_docx(tbl_cox, file = file.path(outdir,'cox_multivariable.docx')),
		silent = TRUE)

## --- PH assumption checks ------------------------------------------------
zph <- cox.zph(cox_multi)
saveRDS(zph, file = file.path(outdir,'cox_zph.rds'))
png(file.path(outdir,'cox_zph_plot.png'), width=1000, height=800)
par(mfrow=c(2,2))
plot(zph)
dev.off()

## --- Interaction test example -------------------------------------------
# Test interaction between rx and age (continuous)
int_model <- update(cox_multi, . ~ . + rx:age)
int_anova <- anova(cox_multi, int_model)
saveRDS(int_model, file = file.path(outdir,'cox_interaction_rx_age.rds'))
saveRDS(int_anova, file = file.path(outdir,'cox_interaction_anova.rds'))

## --- Model selection (stepwise AIC as a simple example) ------------------
step_mod <- MASS::stepAIC(cox_multi, direction = 'both', trace = FALSE)
saveRDS(step_mod, file = file.path(outdir,'cox_stepAIC.rds'))
tbl_step <- gtsummary::tbl_regression(step_mod, exponentiate = TRUE)
try(export_gtsummary_docx(tbl_step, file = file.path(outdir,'cox_stepAIC.docx')),
		silent = TRUE)

## --- Concordance (C-index) and calibration check -------------------------
cindex <- survConcordance(survobj ~ predict(cox_multi), data = dat)
saveRDS(cindex, file = file.path(outdir,'cindex.rds'))

## --- Time-dependent covariate example (demonstration) --------------------
# Create a simple time-dependent covariate using age as if it changes at time=100
if(max(dat$time, na.rm=TRUE) > 100){
	td <- survSplit(dat, cut = 100, end = 'time', event = 'status', start = 'tstart')
	td$age_gt100 <- ifelse(td$tstart >= 100, td$age, NA)
	# fit counting-process style cox
	cox_td <- coxph(Surv(tstart, time, status) ~ age_gt100 + sex + rx + nodes + differ, data = td)
	saveRDS(cox_td, file = file.path(outdir,'cox_time_dependent_example.rds'))
}

## --- Simple forest plot --------------------------------------------------
gg <- ggforest(cox_multi, data = dat)
ggsave(filename = file.path(outdir,'cox_forest.png'), plot = gg, width = 8, height = 6)

cat('Colon analysis finished. Outputs in', outdir, '\n')
