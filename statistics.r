library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ARTool) # non-parametric anova 
library(emmeans)

# Rev Pt 1.4  cross predictions for usv and starling  ---------------------------------------------------
rm(list = ls())
cross_predict_path <- '~/rds/projects/kozlovlab_rds2/live/jphys/starling_usv_cross_predict/cross_predictions.csv' 

cross_predict_usv_star <- read_csv(cross_predict_path) %>%
    rename(cc_train_usv_test_star = cc_train_usv_test_usv, # correct variable names 
        cc_train_usv_test_usv = cc_train_usv_test_star, 
        cc_train_star_test_star = cc_train_star_test_usv, 
        cc_train_star_test_usv = cc_train_star_test_star) %>% 
    mutate(cc_train_usv_test_usv = cc_train_usv_test_usv/cc_upper_bound_usv, 
        cc_train_usv_test_star = cc_train_usv_test_star/cc_upper_bound_starling, 
        cc_train_star_test_star = cc_train_star_test_star/cc_upper_bound_starling, 
        cc_train_star_test_usv = cc_train_star_test_usv/cc_upper_bound_usv) %>% 
    extract(file, "unit_name", "/(20.*)_MNE_model.mat") %>% 
    drop_na() %>% 
    dplyr::select(unit_name, contains("train")) %>% 
    pivot_longer(
        cols = cc_train_usv_test_usv:cc_train_star_test_usv, 
        names_to = c("train_set", "test_set"),
        names_pattern = "cc_train_(.*)_test_(.*)",
        values_to = "correlation") %>% 
    mutate(type = ifelse(train_set == test_set, "within-stim", "across-stim"), train_set = as.factor(train_set)) %>% 
    mutate(type = factor(type, levels=c("within-stim", "across-stim")), unit_name =as.factor(unit_name))

ggboxplot(cross_predict_usv_star, x="type", y="correlation",width=0.5) + 
    geom_jitter(width=0.1) + 
    facet_wrap(~train_set) + 
    ylab("normalized_correlation") + 
    theme(text=element_text(size=15),
        strip.text.x = element_text(size=15),  
        axis.title.x = element_blank())

ggsave("plots/cross_predict_usv_star.png")

cross_m <- art(correlation ~ type*train_set + Error(unit_name), data=cross_predict_usv_star)
anova(cross_m)

#Table Type: Repeated Measures Analysis of Variance Table (Type I) 
#Model: Repeated Measures (aov)
#Response: art(correlation)
#
#                 Error Df Df.res   F value     Pr(>F)    
#1 type           Withn  1    159 560.09710 < 2.22e-16 ***
#2 train_set      Withn  1    159   0.73037    0.39405    
#3 type:train_set Withn  1    159  22.74244  4.159e-06 ***
#---
#Signif. codes:   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Rev Pt 1.5 Effect of dataset size ----------------------------------------------------------------------------
rm(list = ls())
summary_files <- list.files(path = '~/rds/projects/kozlovlab_rds2/live/jphys/subset_training_data', pattern="*.csv", full.names = TRUE)

dat <- as_tibble(do.call(rbind, lapply(summary_files, read_csv, id="path", col_select = 
    c("filename", "correlation_coefficient_usv", "correlation_coefficient_starling", 
    "significant_excitatory_features_usv", "significant_inhibitory_features_usv", 
    "significant_excitatory_features_starling", "significant_inhibitory_features_starling")))) %>% 
    extract(path, "percentData", ".*summary_(.*).csv")

# effect of dataset size on correlation coefficients 
cc <- dat %>% 
    extract(filename, "unit_name", "/(20.*)_split_[0-9]+_MNE_model.mat") %>% 
    dplyr::select(unit_name, percentData, starts_with("correlation")) %>% 
    pivot_longer(
        cols = correlation_coefficient_usv:correlation_coefficient_starling,
        names_to = "stim",
        names_pattern = "correlation_coefficient_(.*)",
        values_to = "correlation") %>% 
        mutate(stim = as.factor(stim), unit_name = as.factor(unit_name), percentData = factor(percentData, levels=c("50", "60", "65", "70", "80", "90", "100")))

# compare CC across different dataset sizes 
star_cc <- cc %>% filter(stim == "starling") %>% drop_na()
star_mod <- art(correlation ~ percentData + Error(unit_name), data=star_cc)
anova(star_mod) # not significant 

usv_cc <- cc %>% filter(stim == "usv") %>% drop_na()
usv_mod <- art(correlation ~ percentData + Error(unit_name), data=usv_cc)
anova(usv_mod)

#Table Type: Repeated Measures Analysis of Variance Table (Type I) 
#Model: Repeated Measures (aov)
#Response: art(correlation)
#
#              Error Df Df.res F value     Pr(>F)    
#1 percentData Withn  6    438  19.185 < 2.22e-16 ***

ggplot(cc,aes(x=percentData, y = correlation)) + 
    geom_boxplot() + 
    geom_jitter(width=0.3) + 
    facet_wrap(~stim) + 
    theme(text=element_text(size=20))

ggsave("plots/correlation_datasetsize.png")

# increasing training data and number of features 
no_sig <- dat %>% 
    extract(filename, "unit_name", "/(20.*)_split_[0-9]+_MNE_model.mat") %>% 
    dplyr::select(unit_name, percentData, starts_with("significant")) %>% 
    pivot_longer(
        cols = significant_excitatory_features_usv:significant_inhibitory_features_starling,
        names_to = c("type", "stim"),
        names_pattern = "significant_(.*)_features_(.*)",
        values_to = "number") %>%
        mutate(stim = as.factor(stim), unit_name = as.factor(unit_name), percentData = factor(percentData, levels=c("50", "60", "65", "70", "80", "90", "100")))

# compare no of excitatory features across different dataset sizes 
usv_nf <- no_sig %>% filter(stim == "usv", type=="excitatory") %>% drop_na()
usv_m2 <- art(number ~ percentData + Error(unit_name), data=usv_nf)
anova(usv_m2)

#Table Type: Repeated Measures Analysis of Variance Table (Type I) 
#Model: Repeated Measures (aov)
#Response: art(number)
#
#              Error Df Df.res F value    Pr(>F)   
#1 percentData Withn  6    438  2.9581 0.0076805 **
#---
#Signif. codes:   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

usv.lm2 <- artlm(usv_m2, "percentData")
usv.emm2 <- emmeans(usv.lm2, ~ percentData)
contrast(usv.emm2 , method="pairwise", adjust="none")

ggplot(no_sig,aes(x=percentData, y =number, fill=type)) + 
    geom_boxplot() + 
    facet_wrap(~stim) + 
    ylim(0, 15) + 
    theme(text=element_text(size=20))

ggsave("plots/no_features_subset_data.png")

# Rev. Pt 1.6 tones and USVs  ------------------------------
rm(list = ls())
results_path <- '~/rds/projects/kozlovlab_rds2/live/jphys/tones_usv/results/csv_format/usv_tones_26nfreqs_12lags.csv' 
summary_file <- read.csv(results_path)

### correlation values 
corr <- summary_file %>% 
    mutate(correlation_coefficient_usv = correlation_coefficient_usv/cc_upper_bound_usv, correlation_coefficient_tones = correlation_coefficient_tones/cc_upper_bound_tones) %>% 
    select(filename, starts_with("correlation")) %>% 
    pivot_longer(
        cols = correlation_coefficient_usv:correlation_coefficient_tones,
        names_to = "stim",
        names_pattern = "correlation_coefficient_(.*)",
        values_to = "correlation") %>%  
        mutate(stim = as.factor(stim))

corr %>% group_by(stim) %>% summarize(mean(correlation), sd(correlation)) # summary statistics

ggboxplot(corr, x="stim", y="correlation", width=0.5) + 
    geom_jitter(width=0.1) + 
    ylab("normalized_correlation") + 
    theme(text=element_text(size=15),
        strip.text.x = element_text(size=15),  
        axis.title.x = element_blank())

m_corr <- wilcox.test(correlation ~ stim, data=corr, paired=TRUE) # p-value = 1.164e-10

ggsave("plots/tones_usv_correlations.png")

### number of features 
no_sig <- summary_file %>% 
    extract(filename, "unit_name", ".*nlags/(.*)_SPK_MNE_model.mat") %>% 
    select(unit_name, starts_with("significant")) %>% 
    pivot_longer(
    cols = significant_excitatory_features_usv:significant_inhibitory_features_tones, 
    names_to = c("type", "stim"),
    names_pattern = "significant_(.*)_features_(.*)",
    values_to = "number") %>%
    mutate(type = as.factor(type), stim = as.factor(stim), unit_name = as.factor(unit_name))

no_sig %>% group_by(stim, type) %>% summarize(mean(number), sd(number)) # summary statistics 

ggboxplot(no_sig, x="type", y="number", fill="type", width=0.5) + 
    geom_jitter(width=0.1) + 
    facet_wrap(~stim) + 
    theme(text=element_text(size=15),
        strip.text.x = element_text(size=15),  
        axis.title.x = element_blank()) + 
    xlab("feature type") +
    theme(legend.position = "none")

ggsave("plots/tones_usv_no_features.png")

m_usv <- wilcox.test(number ~ type, data=no_sig %>% filter(stim == "usv"), paired=TRUE)
m_tones <- wilcox.test(number ~ type, data=no_sig %>% filter(stim == "tones"), paired=TRUE)
m_exc <- wilcox.test(number ~ stim, data=no_sig %>% filter(type == "excitatory"), paired=TRUE)

# cross predictions for tones and usvs 
rm(list = ls())
summary_file <- read.csv('~/rds/projects/kozlovlab_rds2/live/jphys/tones_usv/results/csv_format/cross_predictions_usv_tones.csv')

cross_predict <- summary_file %>% 
    extract(file, "unit_name", ".*mne_models/(.*)_model.mat") %>%
    mutate(cc_train_usv_test_usv = cc_train_usv_test_usv/cc_upper_bound_usv, 
            cc_train_usv_test_tones = cc_train_usv_test_tones/cc_upper_bound_tones, 
            cc_train_tones_test_tones = cc_train_tones_test_tones/cc_upper_bound_tones, 
            cc_train_tones_test_usv = cc_train_tones_test_usv/cc_upper_bound_usv) %>% 
    dplyr::select(unit_name, contains("train")) %>% 
    pivot_longer(
        cols = cc_train_usv_test_usv:cc_train_tones_test_usv, 
        names_to = c("train_set", "test_set"),
        names_pattern = "cc_train_(.*)_test_(.*)",
        values_to = "correlation") %>% 
    mutate(type = ifelse(train_set == test_set, "within-stim", "across-stim"), train_set = as.factor(train_set)) %>% 
    mutate(type = factor(type, levels=c("within-stim", "across-stim")), unit_name = as.factor(unit_name), train_set = as.factor(train_set))

cross_predict %>% group_by(train_set, type) %>% summarize(mean(correlation), sd(correlation))

ggboxplot(cross_predict, x="type", y="correlation",width=0.5) + 
    geom_jitter(width=0.1) + 
    facet_wrap(~train_set) + 
    ylab("normalized_correlation") + 
    theme(text=element_text(size=15),
    strip.text.x = element_text(size=15),  
    axis.title.x = element_blank())

ggsave("plots/cross_predict_usv_tones.png")

cross_m <- art(correlation ~ type*train_set + Error(unit_name), data=cross_predict)
anova(cross_m)

#Table Type: Repeated Measures Analysis of Variance Table (Type I) 
#Model: Repeated Measures (aov)
#Response: art(correlation)

#                 Error Df Df.res  F value     Pr(>F)    
#1 type           Withn  1     99 252.1609 < 2.22e-16 ***
#2 train_set      Withn  1     99  46.2048 8.1056e-10 ***
#3 type:train_set Withn  1     99   1.1616    0.28376    
#---
#Signif. codes:   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

