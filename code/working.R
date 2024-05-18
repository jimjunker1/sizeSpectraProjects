source("./code/paretocounts.R")
# devtools::install_github("jswesner/isdbayes")
# pacman::p_load_gh("jswesner/isdbayes")
# library(isdbayes)
library(tidyverse)
library(brms)
library(stats)
library(utils)
library(junkR)
library(readxl)
'%ni%' <- Negate('%in%')
rstan::rstan_options(auto_write = TRUE)

#
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

new_names = list(c("Macropelopia","Micropsectra sp.","Eukiefferiella","Orthocladius oblidens","Simulium vittatum", "Corynoneura sp.","Chaetocladius dentiforceps","Midge indet.",
                   "Chaetocladius dentiforceps", "Cricotopus sylvestris","Diamesa bertrami","Rheocricotopus effusus","Thysanoptera"),
                 c("Macropelopia", "Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Lumbricidae","Tub. 2",
                   "Orthocladius fridgidus","Micropsectra sp.", "Eukiefferiella", "Diamesa bohemani_zernyi","Micropsectra sp.","Midge indet.", "Orthocladius fridgidus",
                   "Orthocladius oblidens","Orthocladius spp.","Thaumaleidae", "Rheocricotopus effusus","Parochlus sp.","Thienemanniella sp.","Thienemanniella sp.",
                   "Orthocladius fridgidus","Thysanoptera"),
                 c("Nais spp.","Tub. 1", "Simulium vittatum", "Cricotopus sylvestris","Thysanoptera", "Midge indet."),
                 c("Simulium vittatum", "Lumbricidae","Lumbricidae","Tub. 1","Tub. 2",
                   "Micropsectra"," Eukiefferiella","Diamesa","Orthocladius fridgidus","Orthocladius oblidens","Rheocricotopus effusus","Thysanoptera"))

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

# prepare data for Ignasi

hengill_bodysizes %>% 
  # dplyr::filter(SITE %ni% c('st7','oh2')) %>% 
  dplyr::mutate(Site_ID = paste("Junker",SITE,DATE, sep = "_"),
                Sampling_effort = "0.115 m2",
                Sampling_efficiency = "Individuals <0.002 underrepresented",
                Body_length = NA,
                Body_length_units = NA,
                Body_weight_units = "mg") %>%
  rename(Species_latinname = "TAXON",
         Body_weight = "MASS") %>%
  dplyr::filter(!is.na(Body_weight)) %>% 
  dplyr::select(-HABITAT, -SITE) %>% 
  dplyr::group_by(Site_ID, Sampling_effort, Sampling_efficiency, Species_latinname, Body_length, Body_weight, Body_length_units, Body_weight_units) %>% 
  dplyr::summarise(Count = n()) %>% 
  group_by(Site_ID, Sampling_effort, Sampling_efficiency) %>% 
  dplyr::mutate(size_classes = n()) %>% 
    dplyr::filter(size_classes > 2) %>% 
  dplyr::select(Site_ID, Sampling_effort, Sampling_efficiency, Species_latinname, Body_length, Body_weight, Body_length_units, Body_weight_units, Count) %>% 
  dplyr::mutate(Observations = "Body length to mass conversion dependent upon taxon, mostly total length. Units are mg ash-free dry mass.") -> Junker_size_info

write.csv(Junker_size_info, "./output/Junker_Iceland_size_info.csv", row.names = FALSE)

obs_string = "stream mean annual temperatures in celcius:"
Junker_size_info %>% 
  ungroup %>% 
  dplyr::select(Site_ID) %>% 
  unique %>% 
  dplyr::mutate(Date = gsub("Junker_.*_(\\d{4}-\\d{2}-\\d{2}).*","\\1", Site_ID),
                yday = lubridate::yday(as.Date(Date)),
                Sampling_Year = gsub("Junker_.*_(\\d{4})-.*","\\1", Site_ID),
                Sampling_Month = month.abb[as.numeric(gsub("Junker_.*_\\d{4}-(\\d{2})-.*","\\1", Site_ID))],
                Sampling_Season = case_when(between(yday, 355,365) | between(yday, 0,80) ~ 'winter',
                                            between(yday, 81,172 ) ~ 'spring',
                                            between(yday, 173,264) ~ 'summer',
                                            TRUE ~ 'fall'),
                Multiple_sampling = 1,
                Geographical_position_3 = case_when(grepl("hver", Site_ID, ignore.case = TRUE) ~ "hver",
                                                    grepl("st6", Site_ID, ignore.case = TRUE) ~ "st6",
                                                    grepl("st9", Site_ID, ignore.case = TRUE) ~ "st9",
                                                    grepl("st14", Site_ID, ignore.case = TRUE) ~ "st14",
                                                    TRUE ~ NA_character_),
                Geographical_latitude = case_when(Geographical_position_3 == "hver" ~  64.024280,
                                                  Geographical_position_3 == "st6" ~  64.055708,
                                                  Geographical_position_3 == "st9" ~ 64.056951,
                                                  Geographical_position_3 == "st14" ~  64.059309,
                                                  TRUE ~ NA_integer_),
                Geographical_longitude= case_when(Geographical_position_3 == "hver" ~ -21.212862,
                                                  Geographical_position_3 == "st6" ~ -21.306785,
                                                  Geographical_position_3 == "st9" ~ -21.308062,
                                                  Geographical_position_3 == "st14" ~  -21.318419,
                                                  TRUE ~ NA_integer_),
                Observations = case_when(Geographical_position_3 == "hver" ~ paste(obs_string,Geographical_position_3,"= 27.2"),
                                         Geographical_position_3 == "st6" ~ paste(obs_string,Geographical_position_3,"= 17.6"),
                                         Geographical_position_3 == "st9" ~ paste(obs_string,Geographical_position_3,"= 11.9"),
                                         Geographical_position_3 == "st14" ~  paste(obs_string,Geographical_position_3,"= 5.0"),
                                         TRUE ~ NA_character_)) %>% 
  dplyr::select(-Date, -yday) %>% 
  ungroup %>% 
  bind_cols(data.frame(
  Article_ID = NA,
  `2nd_screening_by` = NA,
  Contacted = 1,
  Publication_source = "Ecology Letters",
  DOI = "10.1111/ele.13608",
  Publication_Type = "J",
  Publication_Year = "2020",
  Reference = "Junker, J. R., W. F. Cross, J. P. Benstead, A. D. Huryn, J. M. Hood, D. Nelson, G. M. Gíslason, and J. S. Ólafsson. 2020. Resource supply governs the apparent temperature dependence of animal production in stream ecosystems. Ecology Letters 23:1809–1819.",
  Geographical_position_1 = "Europe",
  Geographical_position_2 = "Iceland",
  Hemisphere = "North",
  Ecosystem = "Freshwater",
  Ecosystem_specification = "Stream",
  Ecosystem_specification_bis = NA,
  Environmental_context = "Temperature",
  Environmental_context_specification = "Gradient",
  Number.trophic.levels = 1,
  Biological_organisation = "Community",
  Species_Type = "Macroinvertebrate",
  Species_Type_specification = NA,
  Species_trophic_group = "Multiple",
  Taxonomic_richness = NA,
  Sampling_methodology = "Surber sampler",
  check.names = FALSE
)) %>% dplyr::select(Article_ID,`2nd_screening_by`,Contacted,Publication_source,DOI, Publication_Type,Publication_Year,Reference,Geographical_position_1,Geographical_position_2,Geographical_position_3,Geographical_latitude,Geographical_longitude,Hemisphere,Sampling_Year,Sampling_Month,Sampling_Season,Multiple_sampling,Ecosystem,Ecosystem_specification,                Ecosystem_specification_bis,Environmental_context,                Environmental_context_specification, Number.trophic.levels,                Biological_organisation,Species_Type,Species_Type_specification,      Species_trophic_group,Taxonomic_richness,Sampling_methodology,
               Observations) -> site_metadata

write.csv(site_metadata, "./output/Junker_Iceland_site_info.csv", row.names = FALSE)

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
sample_m = purrr::map(production_files, ~read_excel(.x, sheet = 'Cobble Gravel', range = 'G35:AK36') %>% pivot_longer(everything(),names_to = 'length_mm', values_to ='mass_mg'))

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
  ungroup %>% 
  bind_cols(data.frame(
  Article_ID = NA,
  `2nd_screening_by` = NA,
  Contacted = 1,
  Publication_source = "Limnology and Oceanography",
  DOI = "10.4319/lo.2014.59.2.0507",
  Publication_Type = "J",
  Publication_Year = "2014",
  Reference = "Junker, J. R. and W. F. Cross. 2014. Seasonality in the trophic basis of a temperate stream invertebrate assemblage: Importance of temperature and food quality. Limnology and Oceanography 59(2):507–518.
",
  Geographical_position_1 = "North_America",
  Geographical_position_2 = "Wyoming, USA",
  Geographical_position_3 = "West Blacktail Deer Creek",
  Geographical_latitude = 44.953720,
  Geographical_longitude= -110.589641,
  Hemisphere = "North",
  Ecosystem = "Freshwater",
  Ecosystem_specification = "Stream",
  Ecosystem_specification_bis = NA,
  Environmental_context = "energy channeling",
  Environmental_context_specification = "aquatic vs terrestrial energy basis",
  Number.trophic.levels = 1,
  Biological_organisation = "Community",
  Species_Type = "Macroinvertebrate",
  Species_Type_specification = NA,
  Species_trophic_group = "Multiple",
  Taxonomic_richness = NA,
  Sampling_methodology = "Surber sampler",
  Observations = "",
  check.names = FALSE
)) %>% 
  dplyr::select(Article_ID,`2nd_screening_by`,Contacted,Publication_source,DOI, Publication_Type,Publication_Year,Reference,Geographical_position_1,Geographical_position_2,Geographical_position_3,Geographical_latitude,Geographical_longitude,Hemisphere,Sampling_Year,Sampling_Month,Sampling_Season,Multiple_sampling,Ecosystem,Ecosystem_specification,                Ecosystem_specification_bis,Environmental_context,                Environmental_context_specification, Number.trophic.levels,                Biological_organisation,Species_Type,Species_Type_specification,      Species_trophic_group,Taxonomic_richness,Sampling_methodology,
                Observations )-> MT_site_metadata

write.csv(MT_site_metadata, "./output/Junker_MT_site_info.csv", row.names = FALSE)

## combine the two projects together
Junker_site_info <- bind_rows(site_metadata, MT_site_metadata)
write.csv(Junker_site_info, "./output/Junker_site_info.csv", row.names = FALSE)

Junker_size_info <- bind_rows(Junker_size_info, sample_mn)
write.csv(Junker_size_info, "./output/Junker_size_info.csv", row.names = FALSE)

### old work for Iceland data. including models.
#####    Create taxon rename lists   #####

old_names <- list(c("Csylvestris"),
                  c("Radix", "RadixBalthica","Radix.Balthica","Radix.balthica"),
                  c("Eukie","Eminor","Eukiefferiella sp."),
                  c("Ooblidens"),
                  c("Tany"),
                  c("Theinem","Theinemanniella"),
                  c("Simuliumvittatum", "Svitt", "Simulium.vittatum"),
                  c("Oligochaeta.A","Oligochaeta.B","Oligo", "Lumb 3","Lumb"), 
                  c("Ofridgidus", "Orthocladius fridgidus","Orthocladius.fridgidus"),
                  c("Rhecricotopus","Rheocricotpus","Rheocricotopus.effusus","Rheocricotopus"),
                  c("Orthoclad"),
                  c("Ceratopegonid","Ceratopogonidae.A","Ceratopogonidae.B"),
                  c("Chaetocladius dentiforceps","Odentiformes"),
                  c("Diamesa.bohemani_zernyi"),
                  c("Diamesa.bertrami"),
                  c("Eukiefferiella.claripennis"),
                  c("Eukiefferiella.minor","Eukiefferiella.mino"),
                  c("Micropsectra", "Micropsectra.sp."),
                  c("Parochlus"),
                  c("Potomaphylax.cingulatus"),
                  c("Prosimulium","Prosimulium.ursinum"),
                  c("Simulium.vernum"),
                  c("Sperchon"),
                  c("Oligochaeta.A"),
                  c("Clinocera"))

new_names <- list("Cricotopus sylvestris",
                  "Radix balthica",
                  "Eukiefferiella sp.",
                  "Orthocladius oblidens",
                  "Macropelopia",
                  "Thienemanniella sp.",
                  "Simulium vittatum",
                  "Lumbricidae",
                  "Orthocladius frigidus",
                  "Rheocricotopus effusus",
                  "Orthocladius spp.",
                  "Ceratopogonid",
                  "Chaetocladius dentiforceps",
                  "Diamesa bohemani_zernyi",
                  "Diamesa bertrami",
                  "Eukiefferiella sp.",
                  "Eukiefferiella sp.",
                  "Micropsectra sp.",
                  "Parochlus sp.",
                  "Potamophylax cingulatus",
                  "Prosimulium ursinum", 
                  "Simulium vernum",
                  "Sperchon glandulosus",
                  "Lumbricidae",
                  "Clinocera stagnalis")

taxon_name_keyval = setNames(rep(new_names, lengths(old_names)), unlist(old_names))



x1 = rparetocounts(mu = -1.8) # `mu` is required wording from brms. in this case it really means the lambda exponent of the ISD
x2 = rparetocounts(mu = -1.5)
x3 = rparetocounts(mu = -1.2)

isd_data = tibble(x1 = x1,
                  x2 = x2,
                  x3 = x3) %>% 
  pivot_longer(cols = everything(), names_to = "group", values_to = "x") %>% 
  group_by(group) %>% 
  mutate(xmin = min(x),
         xmax = max(x)) %>% 
  group_by(group, x) %>% 
  add_count(name = "counts")


sim_fit1 = brm(x | vreal(counts, xmin, xmax) ~ group, 
               data = isd_data,
               stanvars = stanvars,
               family = paretocounts(),
               chains = 1, iter = 1000)



###
isd_df = hengill_bodysizes %>%
  dplyr::filter(MASS >= 0.002) %>%
  select(SITE, DATE, dw = "MASS") %>%
  group_by(SITE, DATE, dw) %>% 
  dplyr::summarise(counts = n()) %>% 
  group_by(SITE, DATE) %>% 
  dplyr::mutate(xmax = max(dw, na.rm = TRUE),
                xmin = min(dw, na.rm = TRUE)) %>% 
  unite("group",SITE,DATE, remove = FALSE) %>% 
  na.omit %>% 
  group_by(group) %>% 
  dplyr::mutate(n = n()) %>% 
  filter(n > 2)
# dplyr::filter(group == "oh2_2011-07-19") %>%
# mutate(counts = counts) %>% 
# named_group_split(group)


fit1 = brm(dw | vreal(counts, xmin, xmax) ~ SITE + (group|DATE/SITE), 
           data = isd_df,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000)

conditional_effects(fit1)




y = isd_df[1:10] %>% 
  purrr::map(~brm(dw | vreal(counts, xmin, xmax) ~ 1, 
                  data = .x,
                  # prior = isd_p,
                  stanvars = stanvars,
                  family = paretocounts(),
                  chains = 1, iter = 500))

y2 = isd_df[11:20] %>% 
  purrr::map(~brm(dw | vreal(counts, xmin, xmax) ~ 1, 
                  data = .x,
                  # prior = isd_p,
                  stanvars = stanvars,
                  family = paretocounts(),
                  chains = 1, iter = 500))


y3 = isd_df[21:30] %>% 
  purrr::map(~brm(dw | vreal(counts, xmin, xmax) ~ 1, 
                  data = .x,
                  # prior = isd_p,
                  stanvars = stanvars,
                  family = paretocounts(),
                  chains = 1, iter = 500))




isd_df[[28]] %>% View()


isd_priors = get_prior(dw | vreal(counts, xmin, xmax) ~ group, 
                       data = isd_df,
                       stanvars = stanvars,
                       prior = isd_p,
                       family = paretocounts())


isd_p = c(prior(normal(0,0.25), class = b),
          prior(normal(-1.5,1), class = Intercept))



fit2 = brm(dw | vreal(counts, xmin, xmax) ~ group, 
           data = isd_df,
           prior = isd_p,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000)

fit3 = update(fit1, newdata = isd_df,
              formula = dw | vreal(counts, xmin, xmax) ~ 1)