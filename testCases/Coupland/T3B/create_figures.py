import matplotlib.pyplot as plt
import numpy as np
import os

# Get the case name and model name form the base directory. The base directory is assumed to be in the format: /path/to/case_name/model_name/
base_dir = os.getcwd()
case_name = base_dir.split('/')[-2]
model_name = base_dir.split('/')[-1]

# Read wall shear stress data from the postProcessing directory
wss = np.loadtxt("postProcessing/surfaces/5000/wallShearStress_plate.raw")
wx = wss[:, 0]
tau_xy = wss[:, 3]

# Sort the data by x-coordinate
sorted_indices = np.argsort(wx)
wx = wx[sorted_indices]
tau_xy = tau_xy[sorted_indices]

# Read velocity from the postProcessing directory
if "kkLOmega" in model_name:
    velocity = np.loadtxt("postProcessing/line_z5/5000/line_kl_kt_omega_U.xy")
    tx = velocity[:, 0]
    k  = velocity[:, 1] + velocity[:, 2]
    u  = velocity[:, 4]
else:
    velocity = np.loadtxt("postProcessing/line_z5/5000/line_k_omega_U.xy")
    tx = velocity[:, 0]
    k  = velocity[:, 1]
    u  = velocity[:, 3]
    
# Sort the velocity data by x-coordinate
sorted_indices = np.argsort(tx)
tx = tx[sorted_indices]
k = k[sorted_indices]
u = u[sorted_indices]

# Interpolate k to match the x-coordinates of wall shear stress
k_interp = np.interp(wx, tx, k)
u_interp = np.interp(wx, tx, u)

# Calculate skin friction coefficient Cf
cf = - 2 * tau_xy / (u_interp**2)

Re_x = u_interp * wx / 1.5e-5  # Assuming kinematic viscosity of air at room temperature

Tu = 100 * np.sqrt(2*k_interp/3) / u_interp  # Turbulence intensity in percentage

# Load experimental data for comparison
exp_data = np.loadtxt("../"+case_name.lower()+"y.dat", skiprows=18, usecols=(1,2,3,4,8))

# Create a cf plot
plt.figure(figsize=(10, 6))
plt.plot(Re_x, cf, label='Skin Friction Coefficient $C_f$', color='blue')
plt.scatter(exp_data[:, 1], exp_data[:, 3], label='Experimental Data', color='red', marker='o', s=50)
plt.xlabel('Reynolds Number $Re_x$')
plt.ylabel('Skin Friction Coefficient $C_f$')
plt.title(f'Skin Friction Coefficient vs Reynolds Number\nCase: {case_name}, Model: {model_name}')
plt.ylim(0,0.01)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'cf_plot_{case_name}_{model_name}.png', dpi=300)

# Create a turbulence intensity plot
plt.clf()
plt.plot(Re_x, Tu, label='Turbulence Intensity $Tu$', color='green')
plt.scatter(exp_data[:, 1], exp_data[:, 4], label='Experimental Data', color='red', marker='o', s=50)
plt.xlabel('Reynolds Number $Re_x$')
plt.ylabel('Turbulence Intensity $Tu$ (%)')
plt.title(f'Turbulence Intensity vs Reynolds Number\nCase: {case_name}, Model: {model_name}')
plt.ylim(0,10)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'tu_plot_{case_name}_{model_name}.png', dpi=300)
