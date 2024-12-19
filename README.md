# Overview

Sessile organisms often rely on their environment to facilitate the dispersal of their seeds or propagules. For plants in particular, species may have multiple vectors of dispersal that operate in series and/or parallel, with each vector contributing individually to their total dispersal, or aggregation of all potential dispersal pathways. Here, we used an individual-based model of population spread to examine the contributions of two serial dispersal vectors (wind as a primary vector, and ants as a secondary vector) to total dispersal in *Carduus nutans*, a widely distributed invasive thistle. We simulated population spread and assessed the relative contributions of both dispersal modes to spread rates, including exploring their interactions with climate warming. Results indicate that an ~0.6 °C higher growing temperature approximately doubled simulated spread rates, regardless of whether ant-driven seed movement and/or predation occurred after wind dispersal. Simulated spread rates were much higher when ant-driven dynamics were removed from the model, however, this was driven more by the exclusion of predation than the exclusion of movement. Individual parameters related to ant-driven seed movement, such as ant nest density or maximum foraging distance, made little contribution to spread rates relative to parameters used to model wind dispersal or ant-driven predation. Our results demonstrate that spread models can be used to examine the relative importance of dispersal modes and related factors even when all parameters are not fully known, providing a way to develop guidance on where limited resources such as money, time, and labor should be spent to ensure accurate data are collected on dispersal pathways that make material contributions to species spread. In our case, information about ant-driven movement of seeds, which is challenging to measure in the field, provided little additional benefit to predictions of plant spatial spread and thus need not be a major focus of future research on spread in this species.

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

**Figure 1** *(.tif)* - Schematic for the algorithm used in the individual-based model of population spread.

**Figure 2** *(.tif)* - Diagram of how the individual-based movel for population spread works.

**Figure 3** *(.tif)* - 

**Figure S1** *(.tif)* - An example of wave movement, using snapshots over time.

**Figure S2** *(.gif)* - An example of wave movement, using bars to visualise density (not featured in manuscript).

**Figure S3** *(.gif)* - An example of wave movement, using lines to visualise density (not featured in manuscript).

## Scripts

**01_Setup** *(.R)* - Code used to define functions and process data.

**02_Stats** *(.R)* - Code used to fit statistical models and estimate parameters.

**03_DemoDisp** *(.R)* - Code used to define functions used that model demographic and dispersal processes.

**04_Wavespeeds** *(.R)* - Code used to simulate an invasion wave.

**05_Scenarios** *(.R)* - Code used to define scenarios for base runs and sensitivity analyses.

**06_Plots** *(.R)* - Code used for plotting figures.

**S1_Extras** *(.R)* - Supplementary code not used in the main analyses.

## Other

**CombinedDispersalMS_v2_Ecology** *(.docx)* - Latest version of the manuscript for this research, submitted to Ecology.

**CombinedDispersalMS_v2_Ecology_Appendix_S1** *(.docx)* - Supplemental material containing a detailed description of methodology used to simulate C. nutans spread rates, and tables for select dispersal statistics.
