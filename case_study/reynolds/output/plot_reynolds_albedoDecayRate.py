'''Basic plot to compare the impact of albedo decay rate parametrization to observations at Reynolds Mountain East site.'''

# modules
import xarray as xr
import matplotlib.pyplot as plt

# specify the files relative to the 'output' folder
file_constant = 'reynolds_constantDecayRate_timestep.nc'
file_variable = 'reynolds_variableDecayRate_timestep.nc'
file_observed = '../evaluation/reynolds_evalData.nc'

# load the simulations
sim_constant = xr.open_dataset(file_constant)
sim_variable = xr.open_dataset(file_variable)
observed = xr.open_dataset(file_observed)

# specify the time period of interest
times = ['2005-09-01','2006-09-01']

# store the observations in temporary variables for clarity
plt_time = observed['time'].sel(time=slice(times[0],times[1]))
plt_snow = observed['zs_sheltered'].sel(time=slice(times[0],times[1])) / 100 # convert cm > m

# make the figure with observations as shaded area and simulations as colored lines
plt.figure(figsize=(20,6))
plt.fill_between(plt_time,plt_snow,color=[0.7,0.7,0.7],label='observations');
sim_variable['scalarSnowDepth'].sel(time=slice(times[0],times[1])).plot(label='Variable albedo decay',color='b'); 
sim_constant['scalarSnowDepth'].sel(time=slice(times[0],times[1])).plot(label='Constant albedo decay',color='r');
plt.legend()
plt.title('Impact of albedo parametrization on snow depth simulations at Reynolds Mountain East')
plt.ylabel('Snow depth [m]');
plt.savefig('reynolds_albedoDecayParametrizationImpact.png', bbox_inches='tight')
plt.close()