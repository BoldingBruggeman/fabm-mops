import cftime
import os
import datetime

import fabmos
import fabmos.transport.tmm
import fabmos.input.riverlist

import pandas as pd

tm_config_dir = "/Users/vsauerland/MITgcm_2.8deg"  # directory with a TM configuration from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
calendar = "360_day"  # any valid calendar recognized by cftime, see https://cfconventions.org/cf-conventions/cf-conventions.html#calendar

script_dir = os.path.dirname(__file__)
fabm_yaml = os.path.join(
    script_dir, "fabm_with_runoff.yaml"
)

domain = fabmos.transport.tmm.create_domain(tm_config_dir)

sim = fabmos.transport.tmm.Simulator(domain, calendar=calendar, fabm=fabm_yaml)

sim.fabm.get_dependency("mole_fraction_of_carbon_dioxide_in_air").set(280.0)

# VS: Jorns suggestion to have same annual air pressure forcing as with atmosp_ files
#     in the PETSc-TMM-MOPS version (see Jorns mail at August 29, 2024)
airp = fabmos.input.from_nc("/Users/vsauerland/tmm_matlab_code/OCMIP_Utils/Data/gasx_ocmip2.nc", "P", decode_times=False)
airp.coords["MONTH_REG"] = fabmos.transport.tmm.climatology_times(calendar)
sim.fabm.get_dependency("surface_air_pressure").set(airp * 101325.0,climatology=True)

# VS: deactivating dilution/concentration of tracers by precipitation/evaporation
#     this keeps uniformly distributet BGC tracers uniformly distributed
#     if the corresponding BGC processes are deactivated
#     (see Jorns mail at August 24, 2024)
# VS: but we have to switch it off for carbon
#     to yield the effect of the term co2emp in BGC_MODEL.F
#     Jorn mailed me how to do it (November 7, 2024, second approach)
for tracer in sim.tracers:
   if tracer.name != 'carbon_c':
      tracer.precipitation_follows_target_cell = True

# Load rivers from https://doi.org/10.1029/96JD00932
rivers = pd.read_fwf(
    os.path.join(script_dir, "rivrstat.txt"),
    skiprows=7,
    usecols=(1, 2, 3, 4, 6),
    names=("name", "country", "lat", "lon", "flow"),
)

# Exclude Arctic rivers as in https://doi.org/10.5194/bg-10-8401-2013
rivers = rivers.loc[rivers.lat <= 60.0]

# Convert river list to gridded [2D] river inputs (level increase in m s-1)
runoff = fabmos.input.riverlist.map_to_grid(sim.T, rivers.lon, rivers.lat, rivers.flow)

# Reintroduce buried phosphorus, weighted by river runoff
# (https://doi.org/10.5194/bg-10-8401-2013)
sim.request_redistribution(
    "det/burial", "runoff/source", runoff, datetime.timedelta(days=360)
)

out = sim.output_manager.add_netcdf_file(
    "output.nc", interval=12, interval_units=fabmos.TimeUnit.HOURS, save_initial=True
)

out.request("runoff_source", "temp", "salt", "wind", "ice", "surface_air_pressure", *sim.fabm.default_outputs, time_average=False) # daily snapshot output (without averaging) 

start = cftime.datetime(2000, 1, 1, calendar=calendar)
stop = cftime.datetime(2001, 1, 1, calendar=calendar)
sim.start(start, timestep = 12 * 3600, transport_timestep = 12 * 3600 )
while sim.time < stop:
    sim.advance()
sim.finish()
