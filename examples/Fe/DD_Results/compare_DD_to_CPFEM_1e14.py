import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Create as Parser to Store the Location of Excel File
parser = argparse.ArgumentParser(description="Argument Parser for Parsing the Location of Excel file")

# Add Arguments
parser.add_argument('--excel_loc', type=str, required=True, help="Path to Excel Result File")
parser.add_argument('--sheet_name', type=str, default="400K_SR10", help="Sheet Name")
parser.add_argument('--fem_csv', type=str, default="NULL", help="CSV File Output From FEM")
args = parser.parse_args()

# Read Data of Dislocation Density
dd_data = pd.read_excel(args.excel_loc, sheet_name=args.sheet_name)

# Read Data of Stress
stress_sheet = args.sheet_name + "_Stress"
stress_data = pd.read_excel(args.excel_loc, sheet_name=stress_sheet)

# Collect C1, C2 and C3 Data
stress_data_c1 = stress_data.iloc[:,:2].to_numpy()
stress_data_c2 = stress_data.iloc[:,3:5].to_numpy()
stress_data_c3 = stress_data.iloc[:,6:8].to_numpy()


dd_data_c1 = dd_data.iloc[0:33, :12].to_numpy()
dd_data_c2 = dd_data.iloc[0:33, 13:25].to_numpy()
dd_data_c3 = dd_data.iloc[0:33, 26:38].to_numpy()

# Read FEM Data if Available
if args.fem_csv != "NULL":
    fem_data = pd.read_csv(args.fem_csv)

# Plot Stress on Subplot 1
plt.subplot(1,2,1)

x = stress_data_c1[:,0]

s1 = stress_data_c1[:,1]/1e6
s2 = stress_data_c2[:,1]/1e6
s3 = stress_data_c3[:,1]/1e6

mean_stress = np.mean([s1, s2, s3], axis=0)
max_stress = np.max([s1, s2, s3], axis=0)
min_stress = np.min([s1, s2, s3], axis=0)

yerr_lower = mean_stress - min_stress
yerr_upper = max_stress - mean_stress

# Stack the yerr data
yerr = np.vstack([yerr_lower, yerr_upper])

plt.errorbar(x, mean_stress, yerr=yerr, color='blue', ecolor='blue', fmt='-o',
             elinewidth=1, capsize=3, label="DDD")

if args.fem_csv != "NULL":
    plt.plot(fem_data.iloc[:,1], fem_data.iloc[:,7], color='red', linestyle="--", linewidth=2, label="DDCP")


plt.subplot(1,2,2)
# Calculate the Errors of Rhom
rhom_C1 = (dd_data_c1[:,2] - dd_data_c1[:,-1])/1e12
rhom_C2 = (dd_data_c2[:,2] - dd_data_c2[:,-1])/1e12
rhom_C3 = (dd_data_c3[:,2] - dd_data_c3[:,-1])/1e12

mean_rhom = np.mean([rhom_C1, rhom_C2, rhom_C3], axis=0)
max_rhom = np.max([rhom_C1, rhom_C2, rhom_C3], axis=0)
min_rhom = np.min([rhom_C1, rhom_C2, rhom_C3], axis=0)

rhom_err_upper = max_rhom - mean_rhom
rhom_err_lower = mean_rhom - min_rhom
rhom_err = np.vstack([rhom_err_lower, rhom_err_lower])



# Calculate the Errors of RhoI
rhoi_C1 = dd_data_c1[:,-1]/1e12
rhoi_C2 = dd_data_c2[:,-1]/1e12
rhoi_C3 = dd_data_c3[:,-1]/1e12

mean_rhoi = np.mean([rhoi_C1, rhoi_C2, rhoi_C3], axis=0)
max_rhoi = np.max([rhoi_C1, rhoi_C2, rhoi_C3], axis=0)
min_rhoi = np.min([rhoi_C1, rhoi_C2, rhoi_C3], axis=0)
rhoi_err_upper = max_rhoi - mean_rhoi
rhoi_err_lower = mean_rhoi - min_rhoi
rhoi_err = np.vstack([rhoi_err_lower, rhoi_err_lower])


# Plot Data
x = dd_data_c1[0:33,0]

plt.errorbar(x, mean_rhom, yerr=rhom_err, color='blue', ecolor='blue', fmt='-o',
             elinewidth=1, capsize=3, label=r"DDD $\rho_m$")

plt.errorbar(x, mean_rhoi, yerr=rhoi_err, color='green', ecolor='green', fmt='-d',
             elinewidth=1, capsize=3, label=r"DDD $\rho_i$")

if args.fem_csv != "NULL":
    plt.plot(fem_data.iloc[1:,1], fem_data.iloc[1:,5]*(24/1e6), color='red', linestyle="--", linewidth=2, label=r"DDCP $\rho_m$")
    plt.plot(fem_data.iloc[1:,1], fem_data.iloc[1:,4]*(24/1e6), color='red', linestyle="-.", linewidth=2, label=r"DDCP $\rho_i$")


# plt.plot(dd_data_c1[:,0], (dd_data_c1[:,2] - dd_data_c1[:,-1])/1e12, color='red',     linestyle="--", linewidth=1, label=args.sheet_name+"_C1" )


# Final Plot Styling
fsize=16
plt.gcf().set_size_inches(16, 6)

plt.subplot(1,2,1)
plt.ylabel(r"$\sigma$, MPa", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.03)
plt.legend(fontsize=fsize, loc="lower right", frameon=True)
plt.tick_params(labelsize=fsize)


plt.subplot(1,2,2)
plt.ylabel(r"$\rho$, $1 \times 10^{12}\ \mathrm{m}^{-2}$", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.03)
plt.ylim(0, 150)
plt.legend(fontsize=fsize, loc="upper right", frameon=True, ncol=2)
plt.tick_params(labelsize=fsize)
plt.show()
