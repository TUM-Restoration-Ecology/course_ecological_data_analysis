The nine csv files of this publication correspond to obtained data from two field experiments testing the effects of management designed to control a native invader in wet grasslands to several aspects of plant community changes.

The data provided correspond to vegetation data obtained from 20 pre-alpine grassland sites and to environmental site conditions of the sites.
Data of plant traits were obtained from several trait databases using the R-package tr8. As these data come from public sources they are only available upon request.

File 1:

Abundances_withoutSenecio_verylow_low_prod
Associated experiment: field experiment on very low and low productive grassland sites conducted to assess suppression effects under different management treatments
Treatments were categorized according to their calculated management intensity

Field experiment conducted at 13 wet grassland sites in S Germany (seven low productive sites; six very low productive sites)
In 2018 treatments consisting of mowing regimes with different management intenisities were implemented to control the native invader; vegetation data were taken at the start (2018) and end (2021) of the experiment

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name
treatment - treatment name used in the paper
treatment_weightedMonth - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
Achi_mill to Vici_crac- abundance occurring plant species; each species was abbreviated with the first four letters of the genus and species name respectively (e.g. Achi_mill refers to Achillea_millefolium)


File 2:

Abundances_withoutSenecio_mod_prod
Associated experiment: field experiment on moderately productive grassland sites conducted to assess suppression effects under different management treatments
Treatments were categorized according to their management intensity as high, medium and low

Field experiment conducted at 7 wet grassland sites in S Germany with moderate productivity
In 2017 treatments consisting of different praxis-oriented combinations of mowing and fertilization were implemented to control the native invader; vegetation data were taken at the start (2017) and end (2020) of the experiment

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name
treatment - treatment name used in the paper
treatment_weightedMonth - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
Agro_stol to Vero_serp- abundance occurring plant species; each species was abbreviated with the first four letters of the genus and species name respectively (e.g. Agro_stol refers to Agrostis_stolonifera)


File 3:

vegetation_corr_verylow_low_prod
Associated experiment: field experiment on very low and low productive grassland sites conducted to assess suppression effects under different management treatments

Field experiment conducted at 13 wet grassland sites in S Germany (seven low productive sites; six very low productive sites)
In 2018 treatments consisting of mowing regimes with different management intenisities were implemented to control the native invader; vegetation data were taken at the start (2018) and end (2021) of the experiment

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name
treatment_intensity - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
categorization_ONLYMT - different categorization idea which was not used in the final paper
presence_WKK - presence absence data of Jacobaea aquatica
Jaco_aqua_cover - estimated cover of Jacobaea aquatica
number_spec_withoutWKK - total number of species per plot without Jacobaea aquatica
number_species - total number of species per plot


File 4:

vegetation_corr_mod_prod
Associated experiment: field experiment on moderately productive grassland sites conducted to assess suppression effects under different management treatments

Field experiment conducted at 7 wet grassland sites in S Germany with moderate productivity
In 2017 treatments consisting of different praxis-oriented combinations of mowing and fertilization were implemented to control the native invader; vegetation data were taken at the start (2017) and end (2020) of the experiment

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name
treatment_intensity - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
categorization_ONLYMT - different categorization idea which was not used in the final paper
presence_WKK - presence absence data of Jacobaea aquatica
Jaco_aqua_cover - estimated cover of Jacobaea aquatica
number_spec_withoutWKK - total number of species per plot without Jacobaea aquatica
number_species - total number of species per plot


File 5:

Change_Abundances
Associated experiment: combined data of both field experiments
Abundances of Jacobaea aquatica Monocots Dicots and total plant number at the start and the end of the experiment were used to calculate change percentages and standardized change values

Variable description
productivity - productivity of the site with (v) very low; (l) low and (m) moderate
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
site - ID of the grassland site
treatment_original - original treatment name expanded by site indication
treatment - treatment name used in the paper
monocot_start - abundance of monocotyledoneous plant species per plot at the start of the experiment in June 2017 2018 respectively
monocot_end - abundance of monocotyledoneous plant species per plot at the end of the experiment in June 2020 2021 respectively
percent_change_monocot - percentage of change in the abundance of monocots from start to end of the experiment
change_monocot - calculated response ratio of change in the abundance of monocots from start to end of the experiment
dicot_start - abundance of dicotyledoneous plant species per plot at the start of the experiment in June 2017 2018 respectively
dicot_end - abundance of dicotyledoneous plant species per plot at the end of the experiment in June 2020 2021 respectively
percent_change_dicot - percentage of change in the abundance of dicots from start to end of the experiment
change_dicot - calculated response ratio of change in the abundance of dicots from start to end of the experiment
dicot_without_WKK_start - abundance of dicotyledoneous plant species without Jacobaea aquatica per plot at the start of the experiment in June 2017 2018 respectively
dicot_without_WKK_end - abundance of dicotyledoneous plant species without Jacobaea aquatica per plot at the end of the experiment in June 2020 2021 respectively
percent_change_dicot_without_WKK - percentage of change in the abundance of dicots without Jacobaea aquatica from start to end of the experiment
change_dicot_without_WKK - calculated response ratio of change in the abundance of dicots without Jacobaea aquatica from start to end of the experiment
WKK_cover_start - abundance Jacobaea aquatica at the start of the experiment in June 2017 2018 respectively
WKK_cover_end - abundance Jacobaea aquatica at the end of the experiment in June 2020 2021 respectively
percent_change_WKK_cover - percentage of change in the abundance of Jacobaea aquatica from start to end of the experiment
change_WKK_cover - calculated response ratio of change in the abundance of Jacobaea aquatica from start to end of the experiment
nr_ssp_start - total number of species per plot at the start of the experiment in June 2017 2018 respectively
nr_ssp_end - total number of species per plot at the end of the experiment in June 2020 2021 respectively
percent_change_nr_ssp - percentage of change in the number of species from start to end of the experiment
change_nr_ssp - calculated response ratio of change in the number of species from start to end of the experiment
nr_ssp_start_without_WKK - number of species without Jacobaea aquatica per plot at the start of the experiment in June 2017 2018 respectively
nr_ssp_end_without_WKK - number of species without Jacobaea aquatica per plot at the end of the experiment in June 2020 2021 respectively
percent_change_nr_ssp_without_WKK - percentage of change in the number of species without Jacobaea aquatica from start to end of the experiment
change_nr_ssp_without_WKK - calculated response ratio of change in the number of species without Jacobaea aquatica from start to end of the experiment


File 6:

Change_Fdis_Fred_allvariations
Associated experiment: combined data of both field experiments
Indices of functional dispersion and functional redundancy at the start and the end of the experiment were used to calculate change percentages and standardized change values; values were calculated for the whole set of plant traits, a reduced and a minimum set of plant traits for comparison

Variable description
productivity - productivity of the site with (v) very low; (l) low and (m) moderate
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
site - ID of the grassland site
treatment_original - original treatment name expanded by site indication
treatment - treatment name used in the paper
percent_change_Fdis_all - percentage of change in functional dispersion from start to end of the experiment using the full set of plant traits
change_Fdis_all - calculated response ratio of change in functional dispersion from start to end of the experiment using the full set of plant traits
percent_change_Fred_all - percentage of change in functional redundancy from start to end of the experiment using the full set of plant traits
change_Fred_all - calculated response ratio of change in functional redundancy from start to end of the experiment using the full set of plant traits
percent_change_Fdis_reduced - percentage of change in functional dispersion from start to end of the experiment using a reduced set of plant traits
change_Fdis_reduced - calculated response ratio of change in functional dispersion from start to end of the experiment using a reduced set of plant traits
percent_change_Fred_reduced - percentage of change in functional redundancy from start to end of the experiment using a reduced set of plant traits
change_Fred_reduced - calculated response ratio of change in functional redundancy from start to end of the experiment using a reduced set of plant traits
percent_change_Fdis_min - percentage of change in functional dispersion from start to end of the experiment using a minimum set of plant traits
change_Fdis_min - calculated response ratio of change in functional dispersion from start to end of the experiment using a minimum set of plant traits
percent_change_Fred_min - percentage of change in functional redundancy from start to end of the experiment using a minimum set of plant traits
change_Fred_min - calculated response ratio of change in functional redundancy from start to end of the experiment using a minimum set of plant traits


File 7:

environ-verylow
Associated experiment: data from very low productive sites
Obtained environmental variables used for NMDS analysis; climate data were taken from the nearest available climate station of the respective site; species abundances and presence data were taken and calculated using the file F1

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
region - ID of regional cluster of the sites
site - ID of the grassland site
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name expanded by site indication
treatment_weightedMonth - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
mean_temp_3month - mean temperature of the three month April, May and June of the respective year
mean_precip_3month - mean precipitation of the three month April, May and June of the respective year
mean_temp_year_DWD - mean temperature of the respective year
mean_precip_year_DWD - mean precipitation of the respective year
temp_longterm_cdc - mean annual longterm temperature
precip_longterm_cdc - mean annual longterm precipitation
altitude - altitude in meter above sea level of the grassland site
grassland_number - potential productivity of the grasslands
soil_type - soil type of the respective site
soil_type_agrolab - soil type classification given from laboratory analysis done by agrolab
pH - pH of the soil obtained by soil probes taken at the start of the experiment
phosphorus - phosphorus in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
potassium - potassium in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
organic_substance - organic substance of the soil in percent
N_total_percent - total percentage of nitrogen in the soil
CN_relation - relation of carbon and nitrogen


File 8:

environ-low
Associated experiment: data from low productive sites
Obtained environmental variables used for NMDS analysis; climate data were taken from the nearest available climate station of the respective site; species abundances and presence data were taken and calculated using the file F1

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
region - ID of regional cluster of the sites
site - ID of the grassland site
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name expanded by site indication
treatment_weightedMonth - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
mean_temp_3month - mean temperature of the three month April, May and June of the respective year
mean_precip_3month - mean precipitation of the three month April, May and June of the respective year
mean_temp_year_DWD - mean temperature of the respective year
mean_precip_year_DWD - mean precipitation of the respective year
temp_longterm_cdc - mean annual longterm temperature
precip_longterm_cdc - mean annual longterm precipitation
altitude - altitude in meter above sea level of the grassland site
grassland_number - potential productivity of the grasslands
soil_type - soil type of the respective site
soil_type_agrolab - soil type classification given from laboratory analysis done by agrolab
pH - pH of the soil obtained by soil probes taken at the start of the experiment
phosphorus - phosphorus in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
potassium - potassium in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
organic_substance - organic substance of the soil in percent
N_total_percent - total percentage of nitrogen in the soil
CN_relation - relation of carbon and nitrogen


File 9:

environ-moderate
Associated experiment: data from moderately productive sites
Obtained environmental variables used for NMDS analysis; climate data were taken from the nearest available climate station of the respective site; species abundances and presence data were taken and calculated using the file F2

Variable description
year_site_plot - ID using a combination of year, experimental site and plot given to distinguish the datalines
year - year in which the data were taken
region - ID of regional cluster of the sites
site - ID of the grassland site
plot - ID of the treatment plot which is a combination of the site and the plot number in the field
treatment_original - original treatment name expanded by site indication
treatment_weightedMonth - treatment intensity calculated using the number of mowing events and their time during the vegetation period; details provided in the appendix of the paper
mean_temp_3month - mean temperature of the three month April, May and June of the respective year
mean_precip_3month - mean precipitation of the three month April, May and June of the respective year
mean_temp_year_DWD - mean temperature of the respective year
mean_precip_year_DWD - mean precipitation of the respective year
temp_longterm - mean annual longterm temperature
precip_longterm - mean annual longterm precipitation
altitude - altitude in meter above sea level of the grassland site
grassland_number - potential productivity of the grasslands
soil_type - soil type of the respective site
soil_type_agrolab - soil type classification given from laboratory analysis done by agrolab
pH - pH of the soil obtained by soil probes taken at the start of the experiment
phosphorus - phosphorus in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
potassium - potassium in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
magnesium - magnesium in mg per 100g or mg per 100ml for peat soil of the soil obtained by soil probes taken at the start of the experiment
organic_substance - organic substance of the soil in percent
N_total_percent - total percentage of nitrogen in the soil
CN_relation - relation of carbon and nitrogen


