# STEDI_python
Hydrological model for small runoff dams

This repository includes code for the hydrological tool STEDI (Spatial tool for estimating dam impacts). For further information about the purpose of this code, please refer to the following two papers:

> Morden R, Horne A, Nathan R, Bond N R (2024), XXXXXXXXXXXXXXXX

> Fowler, K., Morden, R., Lowe, L., Nathan, R., 2015. Advances in assessing the impact of hillside farm dams on streamflow. Australian Journal of Water Resources 19. https://doi.org/10.1080/13241583.2015.1116182

## Installation
I am a complete GitHub newbie as of August 2022. There is no prepared package to install, this is simply a repository of my code. If it looks useful, please download it and reference this repository or my paper Morden et al (2024).

## Overview
STEDI is a tool which was originally developed as a windows executable program in the early 2000's. It estimates the impact of small dams on downstream flows. This code represents a somewhat simplified version of the algorithm.

STEDI works by calculating a water balance around each dam in the catchment, including inflows, rainfall and evaporation on the surface, on-site consumption, and downstream spills. It is the difference between the inflow and downstream spills which represents the flow impact of each dam.

## Inputs

The code requires the following input data:
* A list of catchments where dam impacts are to be calculated. (See example data `STEDI_catlist.csv`)
* A list of dams in each catchment. All dams in all catchments can be listed in a single file. (See example data `STEDI_dams.csv`)
* A timeseries of estimated catchment runoff. This is equivalent to downstream gauged flow assuming the dams were not present. (See example data `STEDI_natflow.csv`)
* A timeseries of rainfall on the surface of the dams. (See example data `STEDI_rain.csv`)
* A timeseries of evaporation from the surface of the dams. (See example data `STEDI_shallowlakeevap.csv`)
* A pattern of on-site demands for irrigation, stock, or domestic use, but NOT including evaporation. (See example data `STEDI_demandpattern.csv`)

Each of these inputs is discussed in more detail below.

### List of catchments
A list of catchments should include the following fields:

* `catID`     a unique alphanumeric ID for each catchment (should match the field name for the respective flow, climate, and demand pattern)
* `area_km2`  area in square kilometres
                    
### Database of dams
The database of dams must have the following index and fields:
                    
* `damID`         unique id for each dam (alphanumeric, used as dataframe index)
* `catID`         name or code of the catchment where the dam is located (should match the field name for the respective flow, climate, and demand pattern, refer 'catlist' above)
* `capacityML`    maximum capacity of the dam in ML (Litres x 10^6)
* `surfarea_m2`   maximum surface area of the dam in square metres
* `catcharea_km2` upstream catchment area of the dam in square kilometres
* `demandfactor`  proportion of the dam extracted or consumed on-site in an average year (between 0 and 1)
                   
A typical demand factor for a small stock dam might be 0.3 to 0.5, whereas a dam built for irrigation might have higher annual consumption nearer 0.7 to 1, or sometimes more if inflows are plentiful all year round.

### Estimated catchment runoff
A file listing the monthly flows for each catchment assuming NO DAMS. In other words, this is the 'unimpacted' flow if no dams were present.

Must be monthly timestep, in units of ML (megalitres, or litres x 10^6). Format should be as follows, but you can always change it to whatever you want:

* `field 1`  date (month, year) in some sort of python interpretable format
* `field 2`  site 1 flow (field name should match the unique id in catlist, see above)
* `field 3`  site 2 flow

etc...
    
### Rain
Rainfall on the surface of each dam in the catchment in mm per month.

File format should be the same as estimated catchment runoff.
                
### Evaporation
Evaporation from the surface of the dam in mm per month.

File format should be the same as estimated catchment runoff.

Note that this is NOT evapo**TRANSPI**ration! This is direct evaporation from the surface of each waterbody in the catchment. An appropriate choice is Mortons shallow lake evaporation where the energy storage within the water is not significant, but where some allowance for cool humid air above the waterbody is needed. Refer Mcmahon et al (2013) equation 23 or Morton (1983b) equation 11.

> McMahon, T.A., Peel, M.C., Lowe, L., Srikanthan, R., McVicar, T.R., 2013. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences 17, 1331–1363. https://doi.org/10.5194/hess-17-1331-2013

> Morton, F. I., 1983. Operational estimates of lake evaporation. Journal of Hydrology 66, 77–100. https://doi.org/10.1016/0022-1694(83)90178-6.
                    
### Pattern of on-site demands
Pattern of on-site demands extracted from the dam each month. Note that this includes water pumped out of the dam for irrigation or other use, and includes dam consumed by stock or wildlife. It does NOT include evaporation which is counted separately. This demand pattern should be a timeseries similar to flow, rainfall, or evaporation.

The exact units are not important as the code pre-processes the timeseries so that the sum of demand in the average year is 1 (ie. the mean of the annual sum of demand = 1). This adjusted pattern is then multiplied by the demand factor (see 'dams' input below) and then by the capacity of the dam. As an equation:

_Demand from dam (ML/month) = dam capacity (ML) x demand factor (no units) x (demand pattern/mean annual demand)_

One method to estimate demand is to calculate net evapotranspiration (eg. short crop et - rain) on a daily timestep, calculate a 2 week moving average, then aggregate to monthly. This gives a first order estimate of crop irrigation demand. Or you can use whatever you like to estimate monthly demand in a given situation.
    
