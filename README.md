# Overview

An model combining primary and secondary dispersal in the invasive thistle *Carduus nutans*, as well as assessing the relative contribution of model components to the spread of this species in a spatial projection model.

<br/>

# Files

## Data

**ThistleData** *(.xlsx)* - Data including flower heights and stages from the experimental thistles.

**SeedDropData** *(.csv)* - Data on seed terminal velocities, used to model seed dispersal (old).

**SeedDropData2** *(.txt)* - Data on seed terminal velocities, used to model seed dispersal (current).

**Weather1** *(.csv)* - Weather data from the study site, with wind speeds used to model seed dispersal.

**Weather2** *(.csv)* -  Weather data from the study site, with wind speeds used to model seed dispersal.

## Scripts

**01_Setup** *(.R)* - Code used to read in the data and fit certain demographic and dispsersal parameters.

**02_DemoDisp** *(.R)* - Code used to define functions used to model demographic and dispersal processes.

**03_WaveElas** *(.R)* - Code used to simulate an invasion wave and calculate wavespeed elasticity.

**S1_Extras** *(.R)* - Extra code for two-dimensional spread that will likely not be implemented.

## Figures

**Spread1** *(.gif)* - An example of wave movement, using bars to visualise density.

**Spread2** *(.gif)* - An example of wave movement, using lines to visualise density.

## Other

**CombinedDispersalMS_v1** *(.docx)* - Manuscript for this research.
