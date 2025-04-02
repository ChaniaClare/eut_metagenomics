library(flopr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)
library(GauPro)
library(stringr)


setwd("/Users/chaniaclare/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Documents/PhD/Barnes-PhD/metagenomics/metadata")

data <- read.csv("Metadata_metagenomics.csv")

my_theme <-   
  theme_classic() +
  theme(legend.position = "none") +
  theme(text=element_text(color="black"),axis.text=element_text(color="black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3, size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(text = element_text(size=9)) +
  theme(axis.title = element_text(face="bold", size=8))+
  theme(plot.title = element_text(hjust = 0.5))



##
## Sampling timeline for each sample
##

sampling_dates <- data %>% 
  select(patient, post_pre_status, days_from_transplant) %>% 
  filter(!(patient == "healthy"))

sampling_timeline <- ggplot() +
  geom_vline(xintercept = 0, color="gray50", linewidth = 0.4, linetype = "dashed") +
  #geom_vline(xintercept = c(-20, 20, 40), color="gray80", linewidth = 0.2) +
  geom_point(data = subset(sampling_dates, post_pre_status == "pre"), aes(x = (days_from_transplant), y = patient), shape = 22, size = 2, color = "gray40", pch = 21, fill = "gray90",) +
  geom_point(data = subset(sampling_dates, post_pre_status == "post"), aes(x = (days_from_transplant), y = patient), shape = 21, size = 2, color = "gray40", pch = 21, fill = "gray90",) +
  #xlab("Sample collection (days)") + 
  #ylab("Patient") +
  my_theme +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  scale_y_discrete(limits = rev)

sampling_timeline
ggsave(plot = sampling_timeline, filename="sampling_timeline.pdf", width = 70, height = 60, units = "mm")


##
## Antibiotic useage heatmap for each patient
##

current_antibiotic_use <- data %>% 
  #select(2,3, 12:23) %>% 
  select(1,3, 12:23) %>% 
  select(order((colnames(.)))) %>% 
  select(sample_id, post_pre_status, everything()) %>% 
  filter(!grepl('H', sample_id))


historic_antibiotic_use <- data %>% 
  select(1,3, 29:40) %>% 
  select(order((colnames(.)))) %>% 
  select(sample_id, post_pre_status, everything()) %>% 
  filter(!grepl('H', sample_id)) %>% 
  rename_with(~str_remove(., '_tot$'))


current_antibiotic_use_melted <- current_antibiotic_use %>% 
  pivot_longer(-c(sample_id, post_pre_status), names_to = 'antibiotic', values_to = 'current_use')

historic_antibiotic_use_melted <- historic_antibiotic_use %>% 
  pivot_longer(-c(sample_id, post_pre_status), names_to = 'antibiotic', values_to = 'historic_use')


combined_data <- current_antibiotic_use_melted %>%
  full_join(historic_antibiotic_use_melted, by = c("sample_id", "antibiotic", "post_pre_status")) %>%
  mutate(status = case_when(
    current_use == 1 ~ "Currently taking",
    historic_use == 1 ~ "Previously exposed",
    TRUE ~ "Never exposed"))

antibiotic_use_heatmap <- ggplot(combined_data, aes(x = antibiotic, y = sample_id)) +
  geom_point(data = subset(combined_data, post_pre_status == "pre"), size = 2, color = "gray90", fill = "gray90", shape = 15) +
  geom_point(data = subset(combined_data, post_pre_status == "post"), size = 2, color = "gray90", fill = "gray90", shape = 16) +
  geom_point(data = subset(combined_data, status == "Previously exposed"), size = 2, color = "red", fill = "white", shape = 21, stroke = 0.5) +
  geom_point(data = subset(combined_data, post_pre_status == "pre" & status == "Currently taking"), size = 2, color = "red", fill = "red", shape = 15) +
  geom_point(data = subset(combined_data, post_pre_status == "post" & status == "Currently taking"), size = 2, color = "red", fill = "red", shape = 16) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
  scale_y_discrete(limits = rev) 
  
antibiotic_use_heatmap
#ggsave(plot = antibiotic_use_heatmap, filename="antibiotic_use_heatmap.pdf", width = 70, height = 70, units = "mm")

###
### Patient age distributions
###

patient_ages <- data %>% 
  select(patient, post_pre_status, age_months, healthy) %>% 
  filter(!grepl('pre', post_pre_status))

patient_age_boxplot <- ggplot() +
  geom_hline(yintercept = c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156), color="gray95", linewidth = 0.2) +
  geom_boxplot(data = patient_ages, aes(x = factor(healthy), y = age_months), linewidth = 0.2, alpha = 0.5, outlier.shape = NA, fill = "gray90") +
  geom_dotplot(data = patient_ages, aes(x = factor(healthy), y = age_months), binaxis='y', stackdir='center',position=position_dodge(width=0.5), 
               dotsize=1, binwidth = 2, method = "histodot", stackratio = 1.5, fill = "gray40", color = "gray40") +
  coord_flip() +
  my_theme +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) 

patient_age_boxplot
ggsave(plot = patient_age_boxplot, filename="patient_age_boxplot.pdf", width = 70, height = 40, units = "mm")


###
### Summary stats
###

summary_stats <- data %>% 
  select(patient, post_pre_status, healthy, sex_female, age_months, viraemia, viraemia_anytime, bacteraemia_anytime, gvhd_acute_grade) %>% 
  filter(!grepl('pre', post_pre_status)) %>% 
  mutate(healthy = factor(healthy, labels = c("Test", "Healthy")))


