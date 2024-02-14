# -*- coding: utf-8 -*-
"""
STEDI algorithm for python (Spatial Tool for Estimating Dam Impacts)

Robert Morden
University of Melbourne
February 2024

This code has been published as part of the following journal paper:
    Morden R, Horne A, Nathan R, Bond N R (2024), 
    XXXXXXXXXXXXXXXX
    

"""
# imports
import pandas as pd
import numpy as np
from numba import jit                                                          # if you are not using NUMBA, just comment this entire line out
    
# =======================================================================================================================
def STEDI_wrapper():

    """
    This subroutine is a wrapper to demonstrate the usage of STEDI.
    
    It loads some basic information for 5 hypothetical catchments and calls the STEDI engine. 
    It requires basic information about dams, on-farm demands, climate, and upstream flow.
    
    Note that this is a SIMPLIFIED VERSION of the algorithm, with no allowance for dams to be
    upstream or downstream of each other, or to have variable surface areas. This code could be
    modified to do that, but it hasn't been necessary so far.
    
    Please refer to Fowler et al (2015) for a more detailed discussion of the algorithm.
    
    Fowler, K., Morden, R., Lowe, L., Nathan, R., 2015. Advances in assessing the impact of
    hillside farm dams on streamflow. Australian Journal of Water Resources 19.
    https://doi.org/10.1080/13241583.2015.1116182

    
    Input files
    -----------
    
    catlist     A list of catchments including the following fields:
                    
                    'catID'    : a unique alphanumeric ID for each catchment (should match the
                                 field name for the respective flow, climate, and demand pattern)
                    'area_km2' : area in km2
                    
    natflow     A file listing the monthly flows for each catchment assuming NO DAMS.
                In other words, this is the 'unimpacted' flow if no dams were present.
                Must be monthly timestep, in units of ML (megalitres, or litres x 10^6).
                Format should be as follows, but you can always change it to whatever you want.
                    
                    field 1  : date (month, year) in some sort of python interpretable format
                    field 2  : site 1 flow (field name should match the unique id in catlist, see above)
                    field 3  : site 2 flow
                    etc...
    
    rain        Rainfall on the surface of each dam in the catchment in mm per month.
                File format should be the same as natflow.
                
    evap        Evaporation from the surface of the dam in mm per month.
                Note that this is NOT evapoTRANSPIration! This is direct evaporation from the
                surface of each waterbody in the catchment. An appropriate choice is
                Mortons shallow lake evaporation where the energy storage within the water is
                not significant, but where some allowance for cool humid air above the waterbody
                is needed. Refer Mcmahon et al (2013) equation 23 or Morton (1983b) equation 11.
                
                McMahon, T.A., Peel, M.C., Lowe, L., Srikanthan, R., McVicar, T.R., 2013.
                    Estimating actual, potential, reference crop and pan evaporation using
                    standard meteorological data: a pragmatic synthesis. Hydrology and Earth
                    System Sciences 17, 1331–1363. https://doi.org/10.5194/hess-17-1331-2013

                Morton, F. I., 1983. Operational estimates of lake evaporation. Journal of
                    Hydrology 66, 77–100. https://doi.org/10.1016/0022-1694(83)90178-6.
                    
    dem         Pattern of on-site demands extracted from the dam each month. Note that this
                includes water pumped out of the dam for irrigation or other use, and includes dam
                consumed by stock or wildlife. It does NOT include evaporation which is counted  
                separately. This demand pattern should be a timeseries similar to flow, rainfall,
                or evaporation. The exact units are not important as the code pre-processes the
                timeseries so that the sum of demand in the average year is 1 (ie. the mean of the
                annual sum of demand = 1). This adjusted pattern is then multiplied by the demand
                factor (see 'dams' input below) and then by the capacity of the dam.
                As an equation:
                    
                    Demand from dam (ML/month) = dam capacity (ML) x demand factor (no units) x (demand pattern/mean annual demand)
                
                One method to estimate demand is to calculate net evapotranspiration (eg. short crop et - rain) 
                on a daily timestep, calculate a 2 week moving average, then aggregate to monthly.
                This gives a first order estimate of crop irrigation demand. Or you can use whatever
                you like to estimate monthly demand in a given situation.
    
    dams        Database of dams in the catchment.
                Must have the following index and fields:
                    
                    'damID'         : unique id for each dam (alphanumeric, used as dataframe index)
                    'catID'         : name or code of the catchment where the dam is located
                                      (should match the field name for the respective flow, climate,
                                      and demand pattern, refer 'catlist' above)
                    'capacityML'    : maximum capacity of the dam in ML
                    'surfarea_m2'   : maximum surface area of the dam in square metres
                    'catcharea_km2' : upstream catchment area of the dam in square kilometres
                    'demandfactor'  : proportion of the dam extracted or consumed on-site in an
                                      average year (between 0 and 1)
                    
                A typical demand factor for a small stock dam might be 0.3 to 0.5, whereas a dam
                built for irrigation might have higher annual consumption nearer 0.7 to 1, or sometimes
                more if inflows are plentiful all year round.
                
    """

    # load basic data
    catlist = pd.read_csv('STEDI_catlist.csv',dtype={'catID':str})                 # database, string index called catID
    catlist.set_index('catID',inplace=True)
    natflow = pd.read_csv('STEDI_natflow.csv',index_col=0,parse_dates=True)        # timeseries, date index in first column
    rain = pd.read_csv('STEDI_rain.csv',index_col=0,parse_dates=True)              # timeseries, date index in first column
    evap = pd.read_csv('STEDI_shallowlakeevap.csv',index_col=0,parse_dates=True)   # timeseries, date index in first column
    dem = pd.read_csv('STEDI_demandpattern.csv',index_col=0,parse_dates=True)      # timeseries, date index in first column
    dams = pd.read_csv('STEDI_dams.csv',dtype={'damID':str,'catID':str})           # database, string index called damID, link to catID
    dams.set_index('damID',inplace=True)
    
    # run stedi
    flow_impacted = STEDI_setup(catlist,natflow,rain,evap,dem,dams)
    
    return flow_impacted

# =======================================================================================================================
def STEDI_setup(catlist,flowin,rain,evap,dem,dams): 
    """
    
    This routine loops through inputs one catchment at a time, then runs
    STEDI separately for each one.
    
    All inputs are panda dataframes. This routine will convert inputs from
    pandas to numpy so that the STEDI engine will work smoothly with NUMBA.
    
    It returns the output to a pandas series at the end.
    
    """
    
    flowout = pd.DataFrame(index=flowin.index,columns=flowin.columns)
    
    # for each catchment
    for irow in range(len(catlist)):
        
        cat_name = catlist.index[irow]
        cat_area = catlist.loc[cat_name,'area_km2']

        # get list of dams for selected catchment
        catdams = dams[dams['catID'] == cat_name].copy()
        
        # determine sum of us areas for selected catchment, adjust if required
        # (ie. the dams cannot impound an area bigger than the total catchment itself, that wouldn't make sense)
        pc_impound = catdams['catcharea_km2'].sum() / cat_area
        if pc_impound > 0.9:
            catdams['catcharea_km2'] = catdams['catcharea_km2'] * (0.9/pc_impound)
        
        # create numpy array
        catdams_array = catdams.loc[:,['capacityML','surfarea_m2','catcharea_km2','demandfactor']].to_numpy()
        
        # adjust irrigation to math STEDI algorithm - average annual pattern should equal 1
        dem_series_ann = dem[cat_name].resample('YS').sum()
        dem_series_adj = dem[cat_name] / dem_series_ann.mean()
        
        # convert timeseries to numpy arrays
        flow_array = flowin[cat_name].to_numpy()
        rain_array = rain[cat_name].to_numpy()
        evap_array = evap[cat_name].to_numpy()
        demd_array = dem_series_adj.to_numpy()
        
        # get start and end of flow series
        
        dtstart = flowin[cat_name].first_valid_index()
        dtend = flowin[cat_name].last_valid_index()
        istart = flowin[cat_name].index.get_loc(dtstart)
        iend = flowin[cat_name].index.get_loc(dtend)
        
        # run monthly calculations
        flowout_array = STEDI_engine(flow_array, rain_array, evap_array, demd_array, catdams_array, cat_area, istart, iend)
        
        # convert back to panda series and add to DataFrame
        flowout[cat_name] = pd.Series(flowout_array,index=flowin.index)
    
    return flowout

# =======================================================================================================================

@jit(nopython=True)                                                            # if you are not using NUMBA, just comment this entire line out
def STEDI_engine(flow,rain,evap,demd,dams,catcharea,istart,iend):
        
    """
    Calculate farm dam impacts on flow (the JIT decorator should get it up to warp speed 9)
    
    Takes the following inputs:
        flow / rain / evap / demd  (1d numpy array, one row per timestep, all units ML or mm)
        dams                       (2d numpy array, 1 row per dam, cols = capacity(ML)/surfarea(m2)/catcharea(km2)/demfactor)
        catcharea                  (scalar, size of catchment in km2)
        istart/iend                (integers, index of first and last useful flow entries in flow array)
    
    Note that all inputs are numpy-based so that NUMBA can make it work super fast
    """
    
    # set up flowout array same as inflows but with all nans
    flowout = np.empty_like(flow)                                           
    flowout[:] = np.nan
    numdams = len(dams)
    
    # set up ENDstore array, set ENDstore = capacity
    estore = np.zeros(numdams)                                              
    estore[:] = dams[:,0]
    
    # f=open('./debug.txt','x', newline='\n')                                  # DEBUG
    # print('numdams = ' + str(numdams),file=f)                                # DEBUG
    # #f.write('capc,surf,ctch,dfac')                                          # DEBUG
    # #np.ndarray.tofile(f,sep=',')                                            # DEBUG
    # print('imonth, netevap, demand, endstore, impact',file=f)                # DEBUG

    for imonth in range(istart,iend+1):                                        # for each month in the useful range
                
        impact_sum = 0.0                                                       # reset impact sum for the month
        # estore_sum = 0.0                                                     # DEBUG
        # demd_sum = 0.0                                                       # DEBUG
        # netevap_sum = 0.0                                                    # DEBUG
        
        for idam in range(0,numdams):                                          # for each dam in the catchment
            
            # assign dam characteristics
            capc = dams[idam,0]                                                
            surf = dams[idam,1]
            ctch = dams[idam,2]
            dfac = dams[idam,3]
           
            sstore = estore[idam]                                              # STARTstore = ENDstore
            
            # calc inflow, rain, evap, demand, outflows, impact
            inflow = ctch / catcharea * flow[imonth]                           # inflow = damarea / catcharea x flow        
            netevap = surf * (evap[imonth]-rain[imonth]) / 1e6                 # netevap = surface area x netevap / 1e6
            demand = capc * dfac * demd[imonth]                                # demand = capacity x demfactor x pattern
            bal = (  sstore + inflow - netevap - demand)                       # raw water balance
            
            estore[idam] = max(0.0, min(bal,capc) )                            # storage must be between 0 and capacity
            outflow = bal - estore[idam]                                       # outflow = spills
            outflow = max(0.0, outflow)                                        # outflow must be positive
            impact = inflow - outflow                                          # damimpact = damin - damout
            impact_sum = impact_sum + impact                                   # aggregate impact sum for each dam for the timestep
            
            # estore_sum = estore_sum + estore[idam]                           # DEBUG
            # demd_sum = demd_sum + demand                                     # DEBUG
            # netevap_sum = netevap_sum + netevap                              # DEBUG
        
        flowout[imonth] = flow[imonth] - impact_sum                            # flow out = flow in - impact
        # print(imonth, netevap_sum, demd_sum, estore_sum, impact_sum,file=f)  # DEBUG
    
    # f.close()                                                                # DEBUG
    return flowout

# =======================================================================================================================

outflow = STEDI_wrapper()