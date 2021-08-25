library(tidyverse)
formants_all <- readRDS("formants_all.RDS")

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
  select(Speaker, Gender, participant_year_of_birth, Word, Vowel, F1_lobanov_2.0, F2_lobanov_2.0) 


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
                  s(participant_year_of_birth, k=10, bs="ad", by=Gender) +
                  s(participant_year_of_birth, k=10, bs="ad") +
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
                  s(participant_year_of_birth, k=10, bs="ad", by=Gender) +
                  s(participant_year_of_birth, k=10, bs="ad") +
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

