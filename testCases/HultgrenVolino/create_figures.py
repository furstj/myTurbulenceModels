#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import os

# Get the case name and model name form the base directory. The base directory is assumed to be in the format: /path/to/case_name/model_name/
base_dir = os.getcwd()
csnm = base_dir.split('/')[-2]
model_name = base_dir.split('/')[-1]

if csnm.split('_')[0] == "Re50k":
    case_name = "Re = 50 000, "
    Uref = 0.75
elif csnm.split('_')[0] == "Re300k":
    case_name = "Re = 300 000, "
    Uref = 4.5
else:
    raise Exception("Unknown Reynolds numer")
    

if csnm.split('_')[1] == "lowTI":
    case_name = case_name + "low turbulence intensity"
    TuPlot_ylim = 1
elif csnm.split('_')[1] == "highTI":
    case_name = case_name + "high turbulence intensity"
    TuPlot_ylim = 5
else:
    raise Exception("Unknown turbulence level")

    
# Read wall shear stress data from the postProcessing directory
wss = np.loadtxt("postProcessing/surfaces/5000/wallShearStress_plate.raw")
wx = wss[:, 0]
tau_xy = wss[:, 3]

# Sort the data by x-coordinate
sorted_indices = np.argsort(wx)
wx = wx[sorted_indices]
tau_xy = tau_xy[sorted_indices]

# Read wall pressure 
pws = np.loadtxt("postProcessing/surfaces/5000/p_plate.raw")
pw  = pws[:, 3]
pw  = pw[sorted_indices]

# Read velocity from the postProcessing directory
# Create a cp plot
uu = np.loadtxt("postProcessing/line_Ue/5000/line_k_omega_p_U.xy")
px = uu[:, 0]
k  = uu[:, 1]
p  = uu[:, 3]
Ue = uu[:, 4]
Ue = np.sqrt(uu[:, 4]**2 + uu[:, 5]**2)

sorted_indices = np.argsort(px)
px = px[sorted_indices]
k  = k[sorted_indices]
p  = p[sorted_indices]
Ue = Ue[sorted_indices]

# Interpolate Ue to match the x-coordinates of wall shear stress
u_interp = np.interp(wx, px, Ue)

# Calculate skin friction coefficient Cf
cf = - 2 * tau_xy / (u_interp**2)

## Load experimental data for comparison
exp_data = np.loadtxt("../data/Cf_%s.txt" % csnm, skiprows=1)
# Create a cf plot
plt.figure(figsize=(10, 6))
plt.plot(wx, cf, label=model_name, color='blue')
plt.scatter(exp_data[:, 0], exp_data[:, 1], label='Experimental Data', color='red', marker='o', s=50)
plt.xlabel('$x/L$')
plt.ylabel('Skin Friction Coefficient $C_f$')
plt.title(f'Skin Friction Coefficient vs $x/L$\nCase: {case_name}')
plt.ylim(-0.005,0.015)
plt.xlim(0.2, 1.4)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'cf_plot_{csnm}_{model_name}.png')

# Pressure coefficient
#p1 = np.interp(1.0, px, p)
#cp = 2*(p - p1) / Uref**2
#cpw = 2*(pw - p1) / Uref**2

cp = 1 - (Ue/Uref)**2

exp = np.loadtxt("../data/Cp_%s.txt" % csnm)
plt.clf()
plt.plot(px, cp, label=model_name, color='magenta')
#plt.plot(wx, cpw, label='$Cp$ at wall', color='magenta', linestyle='dotted')
plt.scatter(exp[:, 1], exp[:, 2], label='Experimental Data', color='red', marker='o', s=50)
plt.xlabel('$x/L$')
plt.ylabel('Pressure Coefficient $C_p$')
plt.title(f'Pressure Coefficient vs $x/L$\nCase: {case_name}')
plt.ylim(-1,0.5)
plt.xlim(0.2, 1.4)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'cp_plot_{csnm}_{model_name}.png')

# External velocity
exp = np.loadtxt("../data/Ue_%s.txt" % csnm)
plt.clf()
plt.plot(px, Ue/Uref, label=model_name, color='green')
plt.scatter(exp[:, 0], exp[:, 1], label='Experimental Data', color='red', marker='o', s=50)
plt.xlabel('$x/L$')
plt.ylabel('Velocity $U_e/U_{ref}$')
plt.title(f'Velocity vs $x/L$\nCase: {case_name}')
plt.ylim(0.5,1.5)
plt.xlim(0.2, 1.4)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'Ue_plot_{csnm}_{model_name}.png')

# Turbulence intensity
Tu = np.sqrt(2*k/3) / Ue * 100
plt.clf()
plt.plot(px, Tu, label=model_name, color='green')
plt.xlabel('$x/L$')
plt.ylabel('Turbulence Intensity $Tu$ (%)')
plt.title(f'Turbulence Intensity vs $x/L$\nCase: {case_name}')
plt.ylim(0,TuPlot_ylim)
plt.xlim(0.2, 1.4)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'Tu_plot_{csnm}_{model_name}.png')


# Velocity profiles
xL = ("0.28", "0.33", "0.39", "0.45", "0.51", "0.57", "0.63",
      "0.69", "0.75", "0.81", "0.88", "0.94", "1.00", "1.06")
data = np.loadtxt("../data/U_%s.txt" % csnm, skiprows=1)
y = data[:, 0]
velocities = data[:, 1:]
Re = Uref * 1 / 1.5e-5
figure, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
for panel, axis in enumerate(axes):
    start = panel * 7
    for local_index in range(7):
        index = start + local_index
        axis.plot(
            velocities[:, index] + 1.5 * local_index,
            y,
            "o",
            color="red",
            markersize=4,
        )
        label = f"x/L={xL[index]}" if local_index == 0 else xL[index]
        axis.text(
            1.5 * local_index + 0.75,
            18.0,
            label,
            ha="center",
            va="center",
            fontsize=14,
        )
        
        results = np.loadtxt("postProcessing/line_st%02i/5000/line_gammaInt_k_U.xy" % (index+1))
        axis.plot(results[:,3] / Uref + 1.5 * local_index, results[:,0] * np.sqrt(Re), color="green", label=model_name)


    axis.set_xlim(-.3, 10.5)
    axis.set_ylim(0, 20.0)
    axis.set_ylabel(r"$y/L\ \mathrm{Re}^{1/2}$", fontsize=16)
    axis.xaxis.set_major_locator(MultipleLocator(1.5))
    axis.xaxis.set_minor_locator(MultipleLocator(0.3))
    axis.yaxis.set_major_locator(MultipleLocator(10))
    axis.yaxis.set_minor_locator(MultipleLocator(1))
    axis.tick_params(axis="both", which="major", labelsize=12, length=7)
    axis.tick_params(axis="both", which="minor", length=4)
    axis.grid(which="major", color="0.75", linewidth=0.8)
    axis.grid(which="minor", color="0.9", linewidth=0.5)

axes[1].set_xlabel(r"$U/U_{\mathrm{ref}}$", fontsize=16)
axes[0].tick_params(axis="x", labelbottom=True)
figure.tight_layout()

plt.savefig(f'U_{csnm}_{model_name}.png')
