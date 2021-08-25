library(emuR)
library(tidyverse)
library(tictoc)

db <- load_emuDB("/Users/peter.gilles/Documents/_Daten/emuDBs/schnëssen_segements_emuDB")
summary(db)
list_levelDefinitions(db)

# all word: KAN & ORT
sl <- query(db, "[ KAN =~ .*]") 
sl2 <- sl %>%
  mutate(ORT = requery_hier(db, seglist = sl, level = "ORT")$labels)
all_KAN <- sl2 %>%
  count(labels, ORT, sort = T)

# get frequencies of all segments: MAU
sl <- query(db, "[ MAU =~ .*]")
sl2 <- sl %>%
  mutate(KAN = requery_hier(db, seglist = sl, level = "KAN")$labels) %>%
  mutate(ORT= requery_hier(db, seglist = sl, level = "ORT")$labels)
all_MAU <- sl2 %>%
  dplyr::select(labels, KAN, ORT) %>%
  count(labels, KAN, ORT, sort = T)


#vowel <- "(i|iː|ə|aː|æ|ɑ|e|eː|o|oː|u|uː)"
vowel <- "(iː|i|eː|e|æ|aː|ɑ|oː|o|uː|u)"
vowel <- "(oː|o)"
#vowel <- c("iː","i","eː","e","æ","aː","ɑ","oː","o","uː","u")
vowel <- c("ɜɪ")

# function for stressed vowels
get_stressed_vowels <- function(vowel) {
  # al KAN segments containing stressed vowels
  query_stressed <- paste0("[ KAN =~ '.*ˈ .*", vowel, " .*' ]")
  sl_kan <- query(db, query = query_stressed)
  sl_kan <- sl_kan %>%
    distinct(labels)
  # create a regex string from the vector, i.e. word|word|...
  sl_kan <- as.vector(sl_kan$labels) %>%
    paste(collapse = "|")
  
  # all MAU vowels only for these KAN words
  query_stressed <- paste0("[ #MAU == '", vowel, "' ^ KAN =~ '(", sl_kan, ")' ]")
  sl_mau <- query(db, query = query_stressed) 
  # change labels for stressed vowel, add ˈ
  toChange <- sl_mau %>% 
    select(session, bundle, level, start_item_seq_idx, attribute, labels) %>% 
    mutate(labels = paste0("ˈ", vowel))
  # update db
  update_itemsInLevel(db, toChange, verbose = TRUE)
  # attach KAN and ORT columns
  #sl_mau %>%
  #  mutate(KAN = requery_hier(db, seglist = sl_mau, level = "KAN", collapse = TRUE)$labels) %>%
  #  mutate(ORT = requery_hier(db, seglist = sl_mau, level = "ORT", collapse = TRUE)$labels)
}

stressed_vowels <- get_stressed_vowels(vowel = "ə")

sl <- query(db, query = "[ORT == schmitt]")

vowel <- query(db, query = "[MAU == ˈə ^ ORT =~ '(schëlder|schëlter|ëmmer|gëtt|schwëster|zëmmer|drénke)' ]") 
temp <- vowel
vowel <- query(db, query = "[MAU == ə ]") 
vowel <- bind_rows(vowel, temp)

vowel2 <- vowel %>%
  mutate(next_sound = requery_seq(db, vowel, offset = 1)$labels)

vowel3 <- vowel2 %>%
  mutate(ORT = requery_hier(db, vowel, level = "ORT")$labels)

vowel <- vowel3

# # column for stopwords
stopwords <- rio::import("../stopwords.txt")
# sl_stressed <- sl2 %>%
   mutate(stopword = if_else(sl2$ORT %in% stopwords$stopword, "stopword", "content word"))

# get formants
td_vowels = get_trackdata(db, vowel3,
                          ssffTrackName = "praatFms", 
                          #                          onTheFlyFunctionName = "forest",
                          #                          onTheFlyParams = list(gender=c("m")),
                          verbose = TRUE)
# normalise the length
td_vowels_norm = normalize_length(td_vowels)
# Achtung: wenn colNames verwendet wird, dann werden zusätzliche Metadatenspalten einfach gelöscht
# td_vowels_norm = normalize_length(td_vowels, colNames = c("T1", "T2"))

# convert Hertz to Bark
# convert T1, T2, T3 to Bark
td_vowels_norm$T1 <- emuR::bark(td_vowels_norm$T1)
td_vowels_norm$T2 <- emuR::bark(td_vowels_norm$T2)
td_vowels_norm$T3 <- emuR::bark(td_vowels_norm$T3)

vowel_midpoints = td_vowels_norm %>% 
  filter(times_norm == 0.4) %>%
  filter(T1 <= 8) %>%
  filter(T2 >= 9) %>%
  group_by(labels, id = word(bundle, 3, sep ="_")) %>%
  mutate(F1 = mean(T1)) %>%
  mutate(F2 = mean(T2)) 

euclid <- vowel_midpoints %>%
  pivot_wider(names_from = c("labels"), values_from = c("F1", "F2"))

# remove outliers
# without_outliers <- vowel_midpoints %>%
#   #  filter(!stopword=="stopword") %>%
#   select(ORT, labels, T1, T2, times_rel) %>%
#   group_by(labels) %>%
#   filter(!find_outliers(T1, T2, keep = 0.95)) 

ggplot(vowel_midpoints) +
  aes(x = T2, y = T1, label = labels, col = labels) +
  geom_text() +
  scale_y_reverse(breaks = seq(from=0, to=9, by=1), minor_breaks = NULL) + 
  scale_x_reverse(breaks = seq(from=3, to=14, by=1), minor_breaks = NULL) + 
  labs(x = "F2", y = "F1") +
  facet_wrap(next_sound ~ labels)
  #ggtitle("Scatter plot of all vowels in the LOD database") +



# old
# stressed vowel
# all individual vowels first will go into a list
sl_stressed <- list()
# loop over the vowels
for (vowel in vowel) {
  # construct the query
  # KAN contains the stress marker ˈ
  query_stressed <- paste0("[#MAU =='", vowel, "' ^ KAN =~ '.*ˈ ", vowel, " .*']")
  # read the query result into a list
  sl_stressed[[vowel]] <- query(db, query = query_stressed)
}
# bind all lists together to one
sl_stressed <- do.call(rbind.data.frame, sl_stressed)
sl_stressed <- sl_stressed %>%
  # insert an additional column with the corresponding word
  mutate(ORT = requery_hier(db, seglist = sl_stressed, level = "ORT", verbose = TRUE)$labels) %>%
  mutate(KAN = requery_hier(db, sl_stressed, "KAN")$labels)

# insert new columns for ORT and KAN in segment list
sl2 <- sl_stressed
sl2$ORT = requery_hier(db, seglist = sl_stressed, level="ORT")$labels
sl2$KAN = requery_hier(db, seglist = sl_stressed, level="KAN")$labels
sl2$stress <- "stressed"

# column for syllable count
sl2$syll_count <- str_count(sl2$KAN, "\\.") + 1


# unstressed vowel
# with regex look around: not preceded by ' 
query_unstressed <- paste0("[#MAU =~'", vowel, "' ^ KAN =~ '.*(?<!ˈ) ", vowel, " .*']")
#sl <- query(db, query = " [MAU ==  'ɑ' ^ #KAN =~ '.*(?<!\') ɑ .*' ]")
sl_unstressed <- query(db, query = query_unstressed)

sl2 <- sl_unstressed
sl2$stress <- "unstressed"
sl2$ORT = requery_hier(db, seglist = sl_unstressed, level="ORT")$labels
sl2$KAN = requery_hier(db, seglist = sl_unstressed, level="KAN")$labels
# column for syllable count
sl2$syll_count <- str_count(sl2$KAN, "\\.") + 1
# column for stopwords
stopwords <- rio::import("stopwords.txt")
sl_unstressed <- sl2 %>%
  mutate(stopword = if_else(sl2$ORT %in% stopwords$stopword, "stopword", "content word"))

# combine to one sl
sl <- bind_rows(sl_stressed, sl_unstressed)

# get formants
td_vowels = get_trackdata(db, sl,
                          ssffTrackName = "praatFms", 
                          #                          onTheFlyFunctionName = "forest",
                          #                          onTheFlyParams = list(gender=c("m")),
                          verbose = TRUE)
# normalise the length
td_vowels_norm = normalize_length(td_vowels)
# achtung: wenn colNames, dann werden zusätzliche Metadatenspalten einfach gelöscht
# td_vowels_norm = normalize_length(td_vowels, colNames = c("T1", "T2"))

# convert Hertz to Bark
# convert T1, T2, T3 to Bark
td_vowels_norm$T1 <- emuR::bark(td_vowels_norm$T1)
td_vowels_norm$T2 <- emuR::bark(td_vowels_norm$T2)
td_vowels_norm$T3 <- emuR::bark(td_vowels_norm$T3)

vowel_midpoints = td_vowels_norm %>% 
  filter(times_norm == 0.4)

# remove outliers
without_outliers <- vowel_midpoints %>%
  #  filter(!stopword=="stopword") %>%
  select(ORT, labels, T1, T2, times_rel) %>%
  group_by(labels) %>%
  filter(!find_outliers(T1, T2, keep = 0.95)) 

ggplot(without_outliers) +
  aes(x = T2, y = T1, label = labels, col = labels) +
  geom_text() +
  scale_y_reverse(breaks = seq(from=0, to=9, by=1), minor_breaks = NULL) + 
  scale_x_reverse(breaks = seq(from=3, to=14, by=1), minor_breaks = NULL) + 
  labs(x = "F2 (Hz)", y = "F1 (Hz)") +
  #ggtitle("Scatter plot of all vowels in the LOD database") +
  facet_wrap(~ stopword)
