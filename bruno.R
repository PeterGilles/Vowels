library(tidyverse)
formants_all <- readRDS("formants_all.RDS")
## GAMMs

#https://nzilbb.github.io/Covariation_monophthongs_NZE/Covariation_monophthongs_analysis.html#GAMM_modelling

vowels_all <- formants_all %>%
  mutate(participant_year_of_birth = recode(Alter, `≤ 24` = 2003, `25 bis 34` = 1993, `35 bis 44` = 1983, `45 bis 54` = 1973, `55 bis 64` = 1963, `65+` = 1953, `55+` = 1963)) %>%
  rename(Vowel = labels, F1_lobanov_2.0 = F1_lob2, F2_lobanov_2.0 = F2_lob2,
         Speaker =id,
         Word = ORT,
         Gender = Geschlecht) %>%
  mutate(Gender = recode(Gender, `Weiblech` = "F")) %>%
  mutate(Gender = recode(Gender, `Männlech` = "M")) %>%
  mutate(participant_year_of_birth = as.integer(participant_year_of_birth)) %>%
  mutate(Speaker = as.factor(Speaker)) %>%
  mutate(Word = as.factor(Word)) %>%
  mutate(Vowel = as.factor(Vowel))

#update the Gender variable to allow for conrast coding
vowels_all$Gender <- as.ordered(vowels_all$Gender)
contrasts(vowels_all$Gender) <- "contr.treatment"

vowels_all <- vowels_all %>%
  arrange(as.character(Speaker))

vowels_all <- vowels_all %>%
  select(Speaker, Gender, participant_year_of_birth, Word, Vowel, F1_lobanov_2.0, F2_lobanov_2.0) %>%
  filter(Vowel %in% c("aː", "æ", "ɑ")) 

# speakers having ALL required vowels
vowels_all <- vowels_all %>% 
  group_by(Speaker) %>% 
  mutate(tokens_per_speaker = n()) %>%
  ungroup() %>%
  group_by(Speaker, Vowel) %>% 
  mutate(vowels_per_speaker = n()) %>%
  filter(tokens_per_speaker >= 20) %>%
  ungroup()


library(mgcv)
library(itsadug)
#create a data frame to store the intercepts from the models, this will initially contain just the speaker names
gam_intercepts.tmp <- vowels_all %>%
  dplyr::select(Speaker) %>%
  distinct()

#loop through the vowels
cat(paste0("Start time:\n", format(Sys.time(), "%d %B %Y, %r\n")))

#vowels_all <- vowels_nz
for (i in levels(factor(vowels_all$Vowel))) {
  
  #F1 modelling
  
  #run the mixed-effects model on the vowel, i.e. if i = FLEECE this will model F1 for FLEECE
  gam.F1 <- bam(F1_lobanov_2.0 ~
                  s(participant_year_of_birth, k=5, bs="cr", by=Gender) +
                  s(participant_year_of_birth, k=5, bs="cr") +
                  Gender +
                  #s(Speech_rate) +
                  s(Speaker, bs="re") +
                  s(Word, bs="re"),
                data=vowels_all %>% filter(Vowel == i),
                discrete=T, nthreads=2)
  
  #extract the speaker intercepts from the model and store them in a temporary data frame
  gam.F1.intercepts.tmp <- as.data.frame(get_random(gam.F1)$`s(Speaker)`)
  
  #assign the model to an object
  assign(paste0("gam_F1_", i), gam.F1)
  
  #save the model summary
  saveRDS(gam.F1, file = paste0("gam_F1_", i, ".rds"))
  
  cat(paste0("F1_", i, ": ", format(Sys.time(), "%d %B %Y, %r"), " ✅\n")) #print the vowel the loop is up to for F1, as well as the start time for the model
  
  #F2 modelling
  
  #run the mixed-effects model on the vowel, i.e. if i = FLEECE this will model F2 for FLEECE
  gam.F2 <- bam(F2_lobanov_2.0 ~
                  s(participant_year_of_birth, k=5, bs="cr", by=Gender) +
                  s(participant_year_of_birth, k=5, bs="cr") +
                  Gender +
                  #s(Speech_rate) +
                  s(Speaker, bs="re") +
                  s(Word, bs="re"),
                data=vowels_all %>% filter(Vowel == i),
                discrete=T, nthreads=2)
  
  #extract the speaker intercepts again, storing them in a separate data frame
  gam.F2.intercepts.tmp <- as.data.frame(get_random(gam.F2)$`s(Speaker)`)
  
  #assign the model to an object
  assign(paste0("gam_F2_", i), gam.F2)
  
  #save the model summary
  saveRDS(gam.F2, file = paste0("gam_F2_", i, ".rds"))
  
  #rename the variables so it clear which one has F1/F2, i.e. this will give F1_FLEECE, F2_FLEECE
  names(gam.F1.intercepts.tmp) <- paste0("F1_", i)
  names(gam.F2.intercepts.tmp) <- paste0("F2_", i)
  
  #combine the intercepts for F1 and F2 and store them in the intercepts.tmp_stress data frame
  gam_intercepts.tmp <- cbind(gam_intercepts.tmp, gam.F1.intercepts.tmp, gam.F2.intercepts.tmp)
  
  cat(paste0("F2_", i, ": ", format(Sys.time(), "%d %B %Y, %r"), " ✅\n")) #print the vowel the loop is up to for F2 , as well as the start time for the model
}

#save the intercepts as a .csv file
write.csv(gam_intercepts.tmp, "gam_intercepts_tmp_new.csv", row.names = FALSE)


vowels_to_plot <- c("ɑ", "aː", "æ")

pred_table <- function (Vowel) {
  mod_name <- paste0("gam_F1_", Vowel)
  return(cbind(Vowel, 
               plot_smooth(get(mod_name), view="participant_year_of_birth", 
                           rm.ranef=T, rug=F, n.grid=500)$fv)
  )
}

gamm_preds_to_plot <- lapply(vowels_to_plot, pred_table) %>%
  do.call(rbind, .) %>%
  mutate(participant_year_of_birth = as.integer(participant_year_of_birth))

saveRDS(gamm_preds_to_plot, "gamm_preds_to_plot.rds")


#make a long version of the intercepts
gam_intercepts.tmp_long <- gam_intercepts.tmp %>%
  pivot_longer(F1_aː:F2_ɑ, names_to = "Vowel_formant", values_to = "Intercept") %>%
  mutate(Formant = substr(Vowel_formant, 1, 2),
         Vowel = substr(Vowel_formant, 4, max(nchar(Vowel_formant)))) %>%
  left_join(vowels_all %>%
              select(Speaker, participant_year_of_birth) %>% distinct()) %>%
  left_join(gamm_preds_to_plot %>% 
              mutate(Formant = "F1") %>% 
              select(participant_year_of_birth, Vowel, Formant, fit, ll, ul))

speakers1 <- gam_intercepts.tmp_long %>%
  filter(Vowel_formant %in% c("F1_aː", "F1_æ", "F1_ɑ")) %>%
  ungroup() %>%
  arrange(participant_year_of_birth) %>%
  filter(Speaker %in% c("1000", "1630", "421", "1004")) %>%
  mutate(letters = c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3)),
         speakers1 = paste0(letters, " (", round(Intercept, 2), ")"))

F1_speaker_intercepts <- gam_intercepts.tmp_long %>%
  filter(Vowel_formant %in% c("F1_aː", "F1_æ", "F1_ɑ"),
         Formant == "F1") %>% #choose the vowels
  mutate(Vowel = factor(Vowel, levels = c("aː", "æ", "ɑ"))) %>% #order the vowels
  ggplot(aes(x = participant_year_of_birth, y = fit, colour = Vowel, fill = Vowel)) +
  geom_point(aes(x = participant_year_of_birth, y = fit + Intercept), size = 1, alpha = 0.2) +
  geom_line() +
  geom_ribbon(aes(ymin=ll, ymax=ul), colour = NA, alpha=0.2) +
  #geom_point(data = speakers1 %>% mutate(fit = ifelse(speakers1 %in% c("D (-0.36)"), fit - 0.05, fit)), aes(x = participant_year_of_birth, y = fit + Intercept), size = 4) +
  #geom_label(data = speakers1 %>% mutate(fit = ifelse(speakers1 %in% c("B (0.23)"), fit - 0.1, fit), fit = ifelse(speakers1 %in% c("D (-0.36)"), fit - 0.1, fit)), aes(x = participant_year_of_birth - 1, y = fit + Intercept, label = speakers1), fill = "white", alpha = 0.5, hjust = 1, show.legend = FALSE) +
  scale_y_reverse() +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  xlab("Participant year of birth") + #x axis title
  ylab("Normalised F1 (Lobanov 2.0)") + #y axis title
  theme_bw() + #general aesthetics
  theme(legend.position = c(0.8, 0.1), #legend position
        legend.direction = "horizontal",
        legend.background = element_rect(colour = "black"),
        axis.text = element_text(size = 14), #text size
        axis.title = element_text(face = "bold", size = 14), #axis text aesthetics
        legend.title = element_blank(), #legend title text size
        legend.text = element_text(size = 14)) #legend text size

F1_speaker_intercepts

ggsave(plot = F1_speaker_intercepts, filename = "F1_speaker_intercepts.png", width = 10, height = 5, dpi = 400)