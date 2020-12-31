library(jsonlite)
library(tidyverse)
source('functions.R')

state_url <- "https://opendata.arcgis.com/datasets/859602b5d427456f821de3830f89301c_11.geojson"
county_url <- "https://opendata.arcgis.com/datasets/5374188992374b318d3e2305216ee413_12.geojson"



df1 <- JSON_fread_state(state_url)

df1 <- df1 %>% rename(datex = Date_reported)


#################check problems with days that have zero:  Oct 17 2020 followed by large numbers on Oct 18 and 19--special cause
vars_common <- c("GEOID","GEO","NAME", "datex")

hosp_subset <- "HOSP_YES"

death_subset <- "DEATHS"

cases_pos_subset <- "POSITIVE"

ICU_subset <- "IC_YES"

HC_subset <- "POS_HC_Y"

subset_list <- c(death_subset,ICU_subset,hosp_subset,HC_subset,cases_pos_subset)

skip_start <- c(ICU_subset,hosp_subset,HC_subset)

df_list <- lapply(1:length(subset_list), make_df)

#Note that ICU counts are NA before 4/1, hosp_N are NA before 3/29, HC_pos_N are NA before 4/13.  This means the first day 
#in the differenced series is NA (which is reasonable, because the first day in the original series typically represents a bolus)
df_all <- do.call(bind_rows,df_list)

