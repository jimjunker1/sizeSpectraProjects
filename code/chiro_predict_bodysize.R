# extract all chironomidae and detect best model to predict size distribution within a date
library(univariateML)
library(junkR)
library(tidyverse)
library(brms)
library(stats)
library(utils)
library(readxl)
'%ni%' <- Negate('%in%')
rstan::rstan_options(auto_write = TRUE)
theme_set(theme_minimal())

## load in Iceland data 
stream_taxon_bodysize = readRDS("./data/stream_taxon_bodysize.rds") %>% bind_rows()

#### rename taxa
streams = list("st6",
               "st9",
               "hver",
               "st14")

old_names = list(c("Tanypodinae","Midge 1","Midge 2","Midge 6","Simulium spp.","Midge 3/2", "Midge 3","ST9 Midge 4","Midge 4",
                   "Midge 5","Midge 7","Midge 8","Springtail"),
                 c("Tanypodinae","Lumb","Oligo","Oligo indet","Oligochaeta","oligochaeta","oligo", "Tubificid 2",
                   "Midge 4/7", "Midge 1","Midge 2", "Midge 3", "Tanytarsus", "Midge", "Midge 4", "Midge 5",
                   "Midge 6","Midge 7","Midge 8","Midge B","Midge C","Midge D","Midge E","Springtail"),
                 c("Lumb 1","Lumb 2", "Simulium spp.","Midge 1","Springtail", "Midge 2"),
                 c("Simulium spp.","Lumb 1", "Lumb 2", "Tubificid 1","Tubificid 2","Midge 1", 
                   "Midge 2","Midge 3","Midge 4","Midge 5","Midge 6","Springtail"))

new_names = list(c("Macropelopia","Chironomidae","Chironomidae","Chironomidae","Simulium vittatum", "Chironomidae","Chironomidae","Chironomidae",
                   "Chironomidae", "Chironomidae","Chironomidae","Chironomidae","Thysanoptera"),
                 c("Macropelopia", "Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Tub. 2",
                   "Chironomidae","Chironomidae", "Chironomidae", "Chironomidae","Chironomidae","Chironomidae", "Chironomidae",
                   "Chironomidae","Chironomidae","Thaumaleidae", "Chironomidae","Chironomidae","Chironomidae","Chironomidae",
                   "Chironomidae","Thysanoptera"),
                 c("Nais spp.","Tub. 1", "Simulium vittatum", "Chironomidae","Thysanoptera", "Chironomidae"),
                 c("Simulium vittatum", "Lumbricidae","Lumbricidae","Tub. 1","Tub. 2",
                   "Chironomidae"," Chironomidae","Chironomidae","Chironomidae","Chironomidae","Chironomidae","Thysanoptera"))

#now do the bulk transformation
#map this across 
#keyval <- future_map2(new_names, old_names, ~setNames(as.list(..1),..2))
keyval <- purrr::map2(new_names, old_names, ~setNames(as.list(..1), unlist(..2)))

hengill_bodysizes = bind_rows(stream_taxon_bodysize) %>% 
  named_group_split(SITE) %>% 
  rlist::list.subset(., c("st6",
                          "st9",
                          "hver",
                          "st14")) %>% 
  purrr::map2(., keyval, ~.x %>% 
                dplyr::mutate(TAXON = recode(TAXON, !!!.y))) %>% 
  bind_rows() 


old_names2 = c("Eukiefferiella","Micropsectra sp.","Midge 3-2","Midge indet.",
                  "Orthocladius oblidens","Orthocladius","Thienemanniella sp.",
                  "Eukiefferiella sp.","Orthocladius frigidus")
new_names2 = rep("Chironomidae", length(old_names2))
keyval2 = purrr::map2(new_names2, old_names2, ~setNames(..1, unlist(..2))) %>% flatten 

chiro_bodysizes =hengill_bodysizes %>% 
  mutate(TAXON = trimws(TAXON),
         TAXON = recode(TAXON, !!!keyval2)) %>% 
  filter(TAXON == "Chironomidae")

chiro_bodysizes %>% 
  ggplot()+
  geom_density(aes(x = log(MASS), y = after_stat(scaled), group = interaction(SITE,DATE), color = interaction(SITE,DATE)))

chiroMList = chiro_bodysizes %>% 
  named_group_split(SITE,DATE) %>% 
  map(~model_select(.x$MASS, models = c("gamma","norm","lnorm","nbinom","pareto","power"), return = 'all'))
  

## west blacktail creek
#DD87:EJ207
# read in data
production_files = list.files(path = "./data/", pattern = "*.xlsm", full.names = TRUE)

sample_names = gsub(".*/WBT (\\w*) .*Prod.xlsm","\\1", production_files)
name_corrections = c("Acentrella sp.",
                     "Ameletus spp.",
                     "Attenella spp.",
                     "Baetis spp.",
                     "Barbaetis spp.",
                     "Brachycentrus spp.",
                     "Chironomidae",
                     "Classenia spp.",
                     "Cinygmula sp.",
                     "Drunella doddsi",
                     "Elmidae (adult)",
                     "Elmidae",
                     "Epeorus spp.",
                     "Glossosomatidae",
                     "Drunella grandis",
                     "Tabanidae",
                     "Hydracarina",
                     "Hydropsychidae",
                     "Leptophlebiidae",
                     "Limno",
                     "Oligochaeta",
                     "Ostracoda",
                     "Pericoma spp.",
                     "Perlodidae",
                     "Rhithrogena spp.",
                     "Rhyacophilidae",
                     "Serratella spp.",
                     "Simuliidae",
                     "Suwallia spp.",
                     "Sweltsa spp.",
                     "Tipulidae",
                     "Zapada spp.")
# get samples abundances
sample_n = purrr::map(production_files, ~read_excel(.x, sheet = 'Cobble Gravel', range = 'DD87:EJ207') %>% 
                        pivot_longer(-c(SAMPLE, DATE), names_to = 'length_mm', values_to = 'n_m2'))

# sample_ab = purrr::map(production_files, ~read_excel(.x, sheet ='Cobble Gravel', range ='A39:A43', col_names = "value", col_types = 'numeric') %>%  .[c(1,3,5),] %>% data.frame(var_name = c('a','b','perc_ash'), value = .)) %>% setNames(., name_corrections)

# get size class to mass
# G36:AK36
sample_m = purrr::map(production_files, ~read_excel(.x, sheet = 'Cobble Gravel', range = 'G35:AK36') %>%
                        pivot_longer(everything(),names_to = 'length_mm', values_to ='mass_mg'))

sample_mn = purrr::map2(sample_n, sample_m, ~left_join(.x,.y, by ='length_mm')) %>% 
  setNames(., nm= name_corrections) %>% 
  bind_rows(.id = "taxon") %>% 
  dplyr::filter(n_m2 >0) %>% 
  dplyr::mutate(mass_mg = round(mass_mg, 6)) %>% 
  group_by(taxon, DATE, length_mm, mass_mg) %>% 
  dplyr::summarise(n_m2 = round(sum(n_m2, na.rm = TRUE))) %>% 
  ungroup %>% 
  dplyr::mutate(Site_ID = paste("Junker_WBT",DATE, sep = "_"),
                Sampling_effort = "0.96 m2",
                Sampling_efficiency = "Individuals <0.002 underrepresented",
                Body_weight_units = "mg",
                Body_length_units = "mm") %>%
  rename(Species_latinname = "taxon",
         Body_weight = "mass_mg",
         Body_length = "length_mm",
         Count = "n_m2") %>% 
  dplyr::select(Site_ID, Sampling_effort, Sampling_efficiency, Species_latinname, Body_length, Body_weight, Body_length_units, Body_weight_units, Count) %>% 
  dplyr::mutate(Observations = "Body length to mass conversion dependent upon taxon, mostly total length. Units are mg ash-free dry mass.")

write.csv(sample_mn, "./output/Junker_MT_size_info.csv", row.names = FALSE)


sample_mn %>% 
  dplyr::select(Site_ID) %>% 
  unique %>% 
  dplyr::mutate(Date = gsub("Junker_WBT_(\\d{4}-\\d{2}-\\d{2}).*","\\1", Site_ID),
                yday = lubridate::yday(as.Date(Date)),
                Sampling_Year = gsub("Junker_WBT_(\\d{4})-.*","\\1", Site_ID),
                Sampling_Month = month.abb[as.numeric(gsub("Junker_WBT_\\d{4}-(\\d{2})-.*","\\1", Site_ID))],
                Sampling_Season = case_when(between(yday, 355,365) | between(yday, 0,80) ~ 'winter',
                                            between(yday, 81,172 ) ~ 'spring',
                                            between(yday, 173,264) ~ 'summer',
                                            TRUE ~ 'fall'),
                Multiple_sampling = 1) %>% 
  dplyr::select(-Date, -yday) %>% 
  ungroup


