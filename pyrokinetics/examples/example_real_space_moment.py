from pyrokinetics import Pyro
import matplotlib.pyplot as plt
import numpy as np
import xrft
import xarray as xr

# Load in file
pyro = Pyro(gk_file="input.cgyro")


a_lte = pyro.local_species.electron.a_lt


# Load in GS2 output data
pyro.gk_code.load_grids(pyro)
pyro.gk_code.load_moments(pyro)

print('Got moments')
data = pyro.gk_output.data
energy = data['moments'].sel(species='electron', drop=True).sel(moment='energy')
energy = energy.where(energy.time > 500, drop=True).mean(dim=['time', 'theta']).sel(kx=energy.kx[1:])

rs_energy = np.real(xrft.fft(energy))
rs_energy.plot(x='freq_kx')
plt.xlabel('$x$')
plt.ylabel(r'$y$')
plt.title(r'$\delta E_e$')
plt.grid()
plt.show()


gradT = energy * energy.kx * 1j
rs_gradT = a_lte - np.real(xrft.fft(gradT)).mean(dim='freq_ky')
rs_gradT.plot()
plt.xlabel('$x$')
plt.ylabel(r'$\omega_{Te}^{eff}$')
plt.grid()
plt.show()


# Include -ky in FFT
# Flip in kx, set ky -> -ky and take conjugate
energy_flip = energy.copy()
energy_flip = energy_flip.assign_coords(ky=-energy_flip.ky).assign_coords(kx=-energy_flip.kx)
energy_flip = energy_flip.reindex(ky=energy_flip.ky[::-1], kx=energy_flip.kx[::-1]).where(energy_flip.ky != -0.0, drop=True)
energy_full = xr.concat([np.conj(energy_flip), energy], dim='ky')


rs_energy = np.real(xrft.fft(energy_full)) * 0.5
rs_energy.plot(x='freq_kx')
plt.show()


gradT = energy_full * energy_full.kx * 1j
rs_gradT = a_lte - np.real(xrft.fft(gradT).mean(dim='freq_ky'))
rs_gradT.plot()
plt.xlabel('$x$')
plt.ylabel(r'$\omega_{Te}^{eff}$')
plt.grid()
plt.show()

np.real(energy_full).plot()
plt.show()

