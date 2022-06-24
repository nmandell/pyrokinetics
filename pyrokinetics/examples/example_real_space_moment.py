from pyrokinetics import Pyro
import matplotlib.pyplot as plt
import numpy as np
import xrft
import xarray as xr

# Load in file
pyro = Pyro(gk_file="input.cgyro")

ky = pyro.numerics.ky
shat = pyro.local_geometry.shat
Delta =  1 / (ky * shat)
length = pyro.cgyro_input['BOX_SIZE'] * Delta
a_lte = pyro.local_species.electron.a_lt


# Load in GS2 output data
pyro.gk_code.load_grids(pyro)
pyro.gk_code.load_moments(pyro)

print('Got moments')
data = pyro.gk_output.data
energy = data['moments'].sel(species='electron', drop=True).sel(moment='energy')
energy = energy.where(energy.time > 500, drop=True).mean(dim=['time', 'theta']).sel(kx=energy.kx[1:])


# Does inverse FFT (c2c along kx, c2r along ky)
rs_energy = xrft.ifft(energy, real_dim='ky')

nx = len(rs_energy.freq_kx.data)
x = np.linspace(-1/2, 1/2, nx) * length

rs_energy = rs_energy.assign_coords(freq_kx=x)

rs_energy.plot(x='freq_kx')
plt.xlabel('$x$')
plt.ylabel(r'$y$')
plt.title(r'$\delta E_e$')
plt.grid()
plt.show()

rs_energy.mean(dim='freq_ky').plot(x='freq_kx')
plt.xlabel('$x$')
plt.title(r'$\delta E_e$')
plt.grid()
plt.show()


gradT = energy * energy.kx * 1j
rs_gradT = a_lte - xrft.ifft(gradT, real_dim='ky').mean(dim='freq_ky')

rs_gradT = rs_gradT.assign_coords(freq_kx=x)
rs_gradT.plot()
plt.vlines([-Delta, 0.0, Delta], ymin=min(rs_gradT.data), ymax=max(rs_gradT.data))
plt.xlabel('$x$')
plt.ylabel(r'$\omega_{Te}^{eff}$')
plt.grid()
plt.show()

