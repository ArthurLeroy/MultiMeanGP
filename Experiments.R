library(tidyverse)
library(MagmaClustR)

require("MMGP.R")

##### Data import #####
# raw_db <- read_csv('PATH_TO_YOUR_DATA')
## The format of raw_db should be a data frame (or tibble) with 4 columns named:
## 'ProbeID', 'ID', 'Input', 'Output'.

# db_coef = read_csv("Data/PedBE_clock_94coefficients.csv") %>% rename(ProbeID = ID)
# db_coef = read_csv("Data/horvath_clock_391coefficients.csv") %>% rename(ProbeID = ID)

db = raw_db %>%
  drop_na()

##### Training #####
# trained_model = training_MM(db, common_hp = TRUE, cv_threshold = 1e-2,
#                             kern_0 = "SE + LIN", kern_i = "SE + LIN")
# write_rds(trained_model, 'Training/horvath_trained_model.rds')

## Because of size restrictions for files on GitHub, the horvath model
## has been saved into two separate files
horvath_model_part1 = read_rds('Training/horvath_trained_model_part1.rds')
horvath_model_part2 = read_rds('Training/horvath_trained_model_part2.rds')

trained_model = c(horvath_model_part1, horvath_model_part2)

##### Prediction #####

## Generate prediction graphs of methylation at 6 years for all individual-CpGs
for(i in unique(db$probeID)){
  first = TRUE
  for(j in unique(db$ID)){
    cat("Pred Probe n°", i, 'ID n°', j, '\n \n')

    ProbeID_pred = i
    ID_pred = j

    db_test = db %>%
      filter(ProbeID %in% ProbeID_pred) %>%
      filter(ID %in% ID_pred) %>%
      filter(Input < 5)

    test_point = db %>%
      filter(ProbeID %in% ProbeID_pred) %>%
      filter(ID %in% ID_pred) %>%
      filter(Input > 5)

    if(first){
      pred = pred_MM(db_test,
                     trained_model = trained_model,
                     grid_inputs = seq(0, 7, 0.1),
                     kern = 'SE',
                     get_hyperpost = TRUE,
                     plot = FALSE)

      hyperpost = pred$hyperpost

      first = FALSE

    } else {
      pred = pred_MM(db_test,
                     trained_model = trained_model,
                     grid_inputs = seq(0, 7, 0.1),
                     kern = 'SE',
                     hyperpost = hyperpost,
                     get_hyperpost = FALSE,
                     plot = FALSE)
    }

    db_train_cpg = db %>% filter(ProbeID %in% ProbeID_pred)

    gg = plot_gp(pred,
                 data = db_test,
                 data_train = db_train_cpg,
                 prior_mean = hyperpost$mean,
                 size_data = 4) +
      geom_point(data = test_point, aes(x = Input, y = Output),
                 size = 3, colour = 'red') +
      xlab("Time") + ylab("Methylation value")

    ggsave(
      paste0('Figures/Methylation/Horvath/indiv_',j,'_CpG_',i,'.pdf'),
      plot = gg, dpi = 600, height=90, width= 180, units="mm"
    )
  }
}

##### Evaluation #####

# res_mmgp = eval(db = db, mod = trained_model)

# write_csv(res, 'Evaluation/res_pred_horvath.csv')
res_mmgp = read_csv('Evaluation/res_pred_horvath.csv')

## Compute the mean (sd) metrics summary
sum_res_mmgp = sum_res %>%
  dplyr::select(- c(ID, ProbeID)) %>%
  group_by(Method) %>%
  summarise_all(list('Mean' = mean, 'SD' = sd), na.rm = TRUE) %>%
  mutate(MSE_Mean = round(MSE_Mean, 4), WCIC_Mean = round(WCIC_Mean, 4),
         MSE_SD = round(MSE_SD, 4), WCIC_SD = round(WCIC_SD, 4)) %>%
  mutate('Mean' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
         'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
  dplyr::select(c(Method, Mean, WCIC))
# write_csv(sum_res, 'Evaluation/summary_pred_age_horvath.csv')

##### Age Prediction #####

# pred_6years = pred_MM_loop(db, trained_model)
# write_csv(pred_6years, 'Prediction/pred_6years_horvath.csv')
pred_6years = read_csv('Prediction/pred_6years_horvath.csv')

intercept = unlist(db_coef[1,2])

## Generate samples of the posterior epigenetic age distribution
size_sample = 1000
pred_age_6years = pred_6years %>%
  filter(Input != 6) %>%
  group_by(ProbeID, ID) %>%
  summarise(Sample = rnorm(size_sample, Mean, sqrt(Var)),
            Mean = Mean,
            Var = Var,
            Index = 1:size_sample) %>%
  left_join(db_coef) %>%
  group_by(ID, Index) %>%
  summarise(Age = anti.trafo(sum(Sample * Coef) + intercept)) %>%
  select(-Index)

# write_csv(pred_age_6years, 'Prediction/pred_age_6years_horvath.csv')
pred_age_6years = read_csv('Prediction/pred_age_6years_horvath.csv')

## Add the true CpG values to the 6 years predictions
db_6years = db %>%
  filter(Input > 5) %>%
  mutate(Input = round(Input, 5)) %>%
  inner_join(pred_6years) %>%
  drop_na()

## Extract the vector of true ages at 6-year data collection
true_age = db %>%
  filter(Input > 5) %>%
  select(ID, Input) %>%
  unique()

## Compute the vector of true epigenetic ages at 6-year data collection
true_pred_age = db %>%
  filter(Input > 5) %>%
  inner_join(db_coef, by = 'ProbeID') %>%
  group_by(ID) %>%
  summarise(Pred_age = unlist(anti.trafo(intercept + sum(Coef * Output))) )

## Plot comparison true age, true pred, Multi-mean pred

for(i in unique(true_age$ID)){
true_age_i = true_age %>% filter(ID == i) %>% pull(Input)
true_pred_age_i = true_pred_age %>% filter(ID == i) %>% pull(Pred_age)

gg = pred_age_6years %>% filter(ID == i) %>%
  ggplot() + geom_density(aes(Age), fill = '#FA9FB5') + theme_classic() +
  geom_vline(aes(xintercept = true_pred_age_i), col = 'red') +
  geom_vline(aes(xintercept = true_age_i), linetype = 'dashed') +
  ggtittle("B")

ggsave(
  paste0('Figures/Epigenetic_Age/pred_age_horvath_ID_',i, '.pdf'),
  plot = gg, dpi = 600, height=180, width= 180, units="mm")
}

##### Plot errors predicted CpG vs true CpG values #####

ggplot(db_6years) +
  geom_point(aes(x = Output, y = Mean), size = 0.2) +
  theme_classic() + guides(col = 'none') +
  xlab('True value') + ylab('Predicted mean')

ggsave(
  paste0('Figures/Error_pred_methylation/horvath_true_pred_CpG_comparison.pdf'),
  plot = gg2, dpi = 600, height=200, width= 200, units="mm")

##### Plot of errors of CpG predictions sorted by increasing uncertainty ####

db_error = db_6years %>%
  mutate(Error = Mean - Output) %>%
  mutate(CI_inf = - qnorm(0.975) * sqrt(Var)) %>%
  mutate(CI_sup = qnorm(0.975) * sqrt(Var)) %>%
  arrange(Var)

gg_error = ggplot(db_error) +
  geom_point(aes(x = Error, y = Var), size = 0.2) +
  geom_ribbon(aes(
      y = Var,
      xmin = CI_inf,
      xmax = CI_sup
    ),
    alpha = 0.3,
    fill = "#FA9FB5"
  ) +
  geom_vline(aes(xintercept = 0), col = "#DB15C1") +
  theme_classic() + xlab("Error (Mean pred - True value)") +
  ylab("Variance of predictions")

ggsave(
  paste0('Figures/Error_uncertainty_methylation/horvath_error_pred_uncertainty_variance.pdf'),
  plot = gg_error, dpi = 600, height=180, width= 360, units="mm")

##### Plot of age prediction errors sorted by increasing uncertainty #####

db_error_age = pred_age_6years %>%
  group_by(ID) %>%
  mutate(CI_inf = quantile(Age, 0.025)) %>%
  mutate(CI_sup = quantile(Age, 0.975)) %>%
  mutate(Var = var(Age)) %>%
  summarise(across(c(Age, CI_inf, CI_sup, Var), mean)) %>%
  left_join(true_age, by = c('ID')) %>%
  left_join(true_pred_age, by = c('ID')) %>%
  mutate(Error = (Age - Input)) %>%
  mutate(Error_pred = (Age - Pred_age)) %>%
  arrange(Age) %>%
  rownames_to_column('Index')

gg_error_age = ggplot(db_error_age) +
  geom_point(aes(x = Pred_age, y = as.numeric(Index)), size = 0.8) +
  geom_point(aes(x = Age, y = as.numeric(Index)), col = 'red' , size = 0.8) +
  geom_ribbon(aes(
    x = Age,
    y = as.numeric(Index),
    xmin = CI_inf,
    xmax = CI_sup
  ),
  alpha = 0.3,
  fill = "#FA9FB5"
  ) +
  # geom_vline(aes(xintercept = 0), col = "#DB15C1") +
  theme_classic() + xlab("Predicted Epigenetic Age") +
  ylab("Subject Index")

ggsave(
  paste0('Figures/Error_uncertainty_methylation/illu_age_prediction_with_comparison.pdf'),
  plot = gg_error_age, dpi = 600, height=100, width= 150, units="mm")

##### Plot mean processes distributions #####

hpost = hyperposterior(trained_model$cg06430061, grid_inputs = seq(0, 7, 0.05))

gg = hpost$pred %>% plot_gp(data_train = db %>%
                         dplyr::filter(ProbeID == 'cg06430061'))

new_gg = gg + geom_line(data = hpost$pred, aes(x = Input, y = Mean), linetype = 'dashed') +
  xlab("Time (Years)") + ylab("Methylation value")

ggsave(
  paste0('Figures/mean_process_CpG_cg06430061.pdf'),
  plot = new_gg, dpi = 300, height=60, width= 120, units="mm")
