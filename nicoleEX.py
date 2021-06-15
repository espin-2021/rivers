# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 19:51:44 2021

@author: Sam
"""

from landlab.components import OverlandFlow, DetachmentLtdErosion
from landlab.io import read_esri_ascii
from landlab import imshow_grid
import numpy as np
import copy
from matplotlib import pyplot as plt


###basin_flag = 'Square' # 'Square' or Long'
###storm_flag = 'Base' # 'Base' or'HigherIntensity' or 'LongerDuration'

### If the basin flag matches one of the two select basins, the filename and
### the outlet link to sample discharge values from for plotting.


###if basin_flag == 'Square':
    ###watershed_dem = 'Square_TestBasin.asc'
    ###link_to_sample = 299
    ###node_to_sample = 300
###else:
    ###watershed_dem = 'Long_TestBasin.asc'
    ###link_to_sample = 149
    ###node_to_sample = 150
    
link_to_sample = 299
node_to_sample = 300

### Reading in the ESRI ASCII DEM given the filename from above, and setting
### the topographic__elevation field.

(rmg, z) = read_esri_ascii("steady_state_topography.txt", name="topographic__elevation")

### remembering the initial watershed elevations for plotting purposes later
z_initial = copy.deepcopy(z) 

#rmg.node_at_link_head[149] #NMG REMOVE REMOVE REMOVE!!!

rmg.set_watershed_boundary_condition(z)

of = OverlandFlow(rmg, steep_slopes=True) # instantiate OverlandFlow object
dle = DetachmentLtdErosion(rmg, K_sp = 1.259162261 * (10**-7)) 
#instantiate DetachmentLtdErosion object


"""
if storm_flag == 'Base':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 7200.
elif storm_flag == 'HigherIntensity':
    starting_precip_mmhr = 10.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 7200.
elif storm_flag == 'LongerDuration':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 14400.
"""
imshow_grid(rmg, z)  # plot the DEM

uplift_rate = 3.170979 * (10**-10) # m/s
starting_precip_ms = 2.77778 * 10 ** -7
storm_duration = 14400.
elapsed_time = 1.0
model_run_time = 43200.0 # s

## Lists for saving data
discharge_at_outlet = []
hydrograph_time = []
incision_at_outlet = []

## Setting initial fields...
rmg['node']['surface_water__discharge'] = np.zeros(rmg.number_of_nodes)
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)

while elapsed_time < model_run_time:
    # Setting the adaptive time step
    of.dt = of.calc_time_step()

    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    if elapsed_time < (storm_duration):
        of.rainfall_intensity =  starting_precip_ms   
    else: # elapsed time exceeds the storm duration, rainfall ceases.
        of.rainfall_intensity = 0.0

    of.overland_flow() # Generating overland flow based on the deAlmeida solution.

    ## Mapping water discharge from links (m^2/s) to nodes (m^3/s) for use
    ## in the DetachmentLtdErosion component.
    node_slope = (of.water_surface_slope[rmg.links_at_node] * rmg.active_link_dirs_at_node)
    incision_Q = np.abs(of.q * rmg.dx)[rmg.links_at_node]
    rmg['node']['surface_water__discharge'] = (incision_Q[np.arange(len(
                                node_slope)), np.argmax(node_slope, axis=1)])
    
    ## Calculating water surface slope from the OverlandFlow component.
    node_slope = node_slope.max(axis=1)
    rmg['node']['water_surface__slope'] = node_slope

    ## Eroding topographic__elevation using DetachmentLtdErosion component.
    dle.erode(of.dt, slope='water_surface__slope')

    ## Updating topographic__elevation after surface was eroded in
    ## DetachmentLtdErosion component.
    rmg['node']['topographic__elevation'] += uplift_rate * of.dt

    ## Append time and discharge to their lists to save data and for plotting.
    hydrograph_time.append(elapsed_time)
    discharge_at_outlet.append(np.abs(of.q[link_to_sample]) * rmg.dx)
    incision_at_outlet.append(dle.I[node_to_sample])

    ## If one wanted to see the time evolve, uncomment below, but this 
    ## produces A LOT of output, as time steps are very short.
    # print(elapsed_time)

    ## Updating elapsed_time
    elapsed_time += of.dt
    
    ## Plotting the hydrograph at the outlet
plt.figure(1)
plt.plot(hydrograph_time, discharge_at_outlet, color='mediumblue')
plt.ylabel('Discharge (cms)')
plt.xlabel('Time (seconds)')
plt.title('Hydrograph at Watershed Outlet')

## Below is if you want to save a figure.
# files=plt.savefig('HydrographAtBasinOutlet.png')

## Plotting the incision rate through time at the outlet
plt.figure(2)
plt.plot(hydrograph_time, incision_at_outlet, color='darkred')
plt.ylabel('incision rate (m/s)')
plt.xlabel('Time (seconds)')
plt.title('Incision at Watershed Outlet')

## Below is if you want to save a figure.
# files=plt.savefig('HydrographAtBasinOutlet.png')

z_diff = z - z_initial # the difference in elevation from start to finish
imshow_grid(rmg, z_diff, limits=(-0.00004, np.max(z_diff)))


### this is just fluvial incision, which will be negative
z_diff = z - z_initial - uplift_rate * model_run_time 
imshow_grid(rmg, z_diff, limits=(-0.00004, 0.0))