import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Create as Parser to Store the Location of Excel File
parser = argparse.ArgumentParser(description="Argument Parser for Parsing the Location of Excel file")

# Add Arguments
parser.add_argument('--fem_csv', type=str, default="NULL", help="CSV File Output From FEM")
args = parser.parse_args()


# Read FEM Data if Available
if args.fem_csv != "NULL":
    fem_data = pd.read_csv(args.fem_csv)

# Plot Stress on Subplot 1
plt.subplot(1,2,1)


if args.fem_csv != "NULL":
    plt.plot(fem_data.iloc[:,1], fem_data.iloc[:,7], color='red', linestyle="--", linewidth=2, label="DDCP")


plt.subplot(1,2,2)

if args.fem_csv != "NULL":
    plt.plot(fem_data.iloc[1:,1], fem_data.iloc[1:,5]*(24/1e6), color='red', linestyle="--", linewidth=2, label=r"DDCP $\rho_m$")
    plt.plot(fem_data.iloc[1:,1], fem_data.iloc[1:,4]*(24/1e6), color='blue', linestyle="-.", linewidth=2, label=r"DDCP $\rho_i$")


# plt.plot(dd_data_c1[:,0], (dd_data_c1[:,2] - dd_data_c1[:,-1])/1e12, color='red',     linestyle="--", linewidth=1, label=args.sheet_name+"_C1" )

max_rhoi = fem_data['rho_i'].max()
max_rhom = fem_data['rho_m'].max()
max_dd = max(max_rhoi, max_rhom)*(24/1e6) + 100

# Final Plot Styling
fsize=16
plt.gcf().set_size_inches(16, 6)

plt.subplot(1,2,1)
plt.ylabel(r"$\sigma$, MPa", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.10)
plt.legend(fontsize=fsize, loc="lower right", frameon=True)
plt.tick_params(labelsize=fsize)


plt.subplot(1,2,2)
plt.ylabel(r"$\rho$, $1 \times 10^{12}\ \mathrm{m}^{-2}$", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.10)
plt.ylim(0, max_dd)
plt.legend(fontsize=fsize, loc="lower right", frameon=True, ncol=1)
plt.tick_params(labelsize=fsize)
plt.tight_layout()
plt.show()
