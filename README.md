# Overview

Sessile organisms often rely entirely on their environment to facilitate the dispersal of their propagules and expand their populations. Unfortunately, dispersal and propagule fate are poorly understood because they are exceptionally challenging and resource-intensive to measure. This is especially true for the many species whose dispersal involves multiple pathways, where each one contributes individually to its total dispersal, or aggregation of all potential dispersal pathways. Here, we address these challenges for the widely distributed invasive thistle *Carduus nutans*. We create an individual-based model of population spread which, when paired with sensitivity analyses, allows us to examine the relative importance of multiple dispersal pathways even when some dispersal and propagule fate parameters are not fully known. We then used this model to assess the contributions of both dispersal pathways in this species (primary dispersal via wind, and secondary dispersal via ants) to simulated spread rates, including exploring how each one is impacted by climate warming. Results indicate that an approximately 0.6°C higher growing temperature approximately doubled simulated spread rates, regardless of whether ant-driven seed movement and/or predation occurred after wind dispersal. Simulated spread rates were much higher when ant-driven dynamics were removed from the model, likely due to the exclusion of predation rather than the exclusion of movement. Parameters used to model ant-driven seed movement, such as ant nest density or maximum foraging distance, made little contribution to spread relative to those used to model wind dispersal or ant-driven predation; given limited resources, future research should thus be focused on quantifying the latter factors. Our results demonstrate that modelling-based approaches can be effectively used to develop guidance on where limited resources such as money, time, and labor should be spent, ensuring accurate data are collected on dispersal pathways that make material contributions to species spread.

<br/>

# Files

## Data

**ThistleData** *(.xlsx)* - Data including flower heights and stages from the experimental thistles.

**SeedDropData** *(.csv)* - Data on seed terminal velocities, used to model seed dispersal (old).

**SeedDropData2** *(.txt)* - Data on seed terminal velocities, used to model seed dispersal (current).

**Weather1** *(.csv)* - Weather data from the study site, with wind speeds used to model seed dispersal.

**Weather2** *(.csv)* -  Weather data from the study site, with wind speeds used to model seed dispersal.

**wv_base** *(.csv)* - Median simulated wavespeeds from the base cases.

**wv_base_aw** *(.csv)* - Summary statistics on all simulated wind and dispersal distances from the base cases.

**wv_base_d** *(.csv)* - Summary statistics on simulated wind and dispersal distance contributions to the wavefront from the base cases.

**wv_elas** *(.csv)* - Median simulated wavespeeds used in the sensitivity analysis for elasticity calculations.

**wv_pert** *(.csv)* - Median simulated wavespeeds used in the sensitivity analysis for perturbation calculations.

## Figures

**Figure 1** *(.pptx)* - Schematic for the algorithm used in the individual-based model of population spread.

**Figure 2** *(.tif)* - Elasticity of median simulated spread rates for various dispersal-related parameters.

**Figure 3** *(.tif)* - Perturbations of median simulated spread rates for various dispersal-related parameters.

**Figure S1** *(.tif)* - An example of simulated wave movement, using snapshots over time.

**Figure S2** *(.gif)* - An example of simulated wave movement, using bars to visualise density (not featured in manuscript).

**Figure S3** *(.gif)* - An example of simulated wave movement, using lines to visualise density (not featured in manuscript).

## Scripts

**01_Setup** *(.R)* - Code used to define functions and process data.

**02_Stats** *(.R)* - Code used to fit statistical models and estimate parameters.

**03_DemoDisp** *(.R)* - Code used to define functions used that model demographic and dispersal processes.

**04_Wavespeeds** *(.R)* - Code used to simulate an invasion wave.

**05_Scenarios** *(.R)* - Code used to define scenarios for base runs and sensitivity analyses.

**06_Plots** *(.R)* - Code used for plotting figures.

**S1_Extras** *(.R)* - Supplementary code not used in the main analyses.

## Other

**CombinedDispersalMS_v6_Ecology** *(.docx)* - Latest version of the manuscript for this research, submitted to Ecology.

**CombinedDispersalMS_v6_Ecology_Appendix_S1** *(.docx)* - Supplemental material containing a detailed description of methodology used to simulate C. nutans spread rates, and tables for select dispersal statistics.
