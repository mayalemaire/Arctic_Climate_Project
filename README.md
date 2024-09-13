**Trends-in-Extreme-Climate-Events-in-the-Arctic-and-Their-Ecological-Consequences**

This repository contains all the code used in the MBiol dissertation title "Trends in Extreme Climate Events in the Arctic and Their Ecological Consequences".

List of the supplementary materials provided:

**Meta-Analysis**

*R code files:*
- meta_analysis_full.R - an R script containing the code used to convert effect sizes to SMDH and build all meta-analytical models
- Kappa_calc.R - an R script used to calculate Cohen's Kappa
- meta_baselines_16_day.R - an R script used to calculate the 90th percentile for 'days' used to identify extreme years in study 16

*Raw data files:*
- Meta_Analysis_Screening_Data.CSV - CSV file containing all the extracted data from the studies included in the meta-analysis
- Meta_Analysis_Screening_Data.xlsx - Excel file containing all the different stages of the data extraction

**Climate Analysis**

*R code files:*
(These are presented in the order in which they are run to complete the analysis)

- 01_nc_to_csv_t2m.R
- 02_nc_to_csv_precipitation.R
- 03_t2m_precip_merger.R - R script used to merge precipitation and temperature data
- 04_ROS_df_1t2m_5mm. R - R script use to detect ROS events
- 05_ROS_grid_plot_fig4.R - R script to produce figure 4 and paired t-tests
- 06_WW_0.99_percentile.R - R script used to determine 0.99 threshold for each grid cell
- 07_WW_gridded_sum_t2m_exceedance.R - R script used to calculated annual cumulative winter warming exceedance
- 08_WW_grid_plot_fig1.R - R script to produce figure 1 and paired t-tests
- 09_precipitation_press_aggregation.R - R script for aggregating precipition to winter seasons data
- 10_ROS_intensity_plot.R - R script to produce figures 5 & 6 as well as linear regressions
- 11_WW_gridded_t2m.R - R script for detecting winter warming events
- 12_t2m_press_aggregation.R - R script for aggregating temperature to winter seasons data
- 13_WW_intensity_plot.R - R script to produce figures 2 & 3 as linear regressions
- 14_Adapted_polar_map.R - R script to produce Supplementary Figure 1


*Raw data files:*
- evidence-map-scope - shapefile of the Arctic regions as established by Martin et al., 2022
- ne_10m_ocean - shapefile of the Ocean used to filter out data points in the Ocean
- ne_10m_coastline - shapefile of coastlines to add onto maps
- ne_10d_graticule - shapefile of 10 degrees graticule for plotting
- All ERA5 data used is available online at https://cds.climate.copernicus.eu/cdsapp#!/home

Thank you for your time!
