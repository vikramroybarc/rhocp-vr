import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Create as Parser to Store the Location of Excel File
parser = argparse.ArgumentParser(description="Argument Parser for Parsing the Location of Excel file")

# Add Arguments
parser.add_argument('--fem_csv', type=str, default=r"NULL", help="CSV File Output From FEM")
parser.add_argument('--exp_file', type=str, default="Exp_results.xlsx", help="Excel File Containing Experimental Results")
parser.add_argument('--exportfig', type=bool, default=False, help="Whether to Export Fig")
args = parser.parse_args()

# Read Experimental File
# data_2021Soares = pd.read_excel(args.exp_file, sheet_name="300KSR0p01_2021Soares")
data_2021Soares = pd.read_excel(args.exp_file, sheet_name="300KSR1_2021Soares")
data_2003Uneishi = pd.read_excel(args.exp_file, sheet_name="300KSR1_2003Unishi")
data_2006Kuroda = pd.read_excel(args.exp_file, sheet_name="300KSR1_2006Kuroda")

# Read FEM Data if Available
if args.fem_csv != "NULL":
    fem_data1 = pd.read_csv(args.fem_csv)


plt.subplot(1,2,1)
plt.plot(data_2021Soares['PStrain2'], data_2021Soares['Stress'], color='blue', linestyle="--", linewidth=2, label="Soares et al. 2021")

if args.fem_csv != "NULL":
    rho_total = (fem_data1['rho_i'][1] + fem_data1['rho_m'][1])*(24/1e8)
    plt.plot(fem_data1['Ep_eff'], fem_data1['stress_xx'], color='red', linestyle=":", linewidth=2, 
             label=rf"DDCP $\rho_t = {rho_total:0.2f} \times 10^{{14}}\ \mathrm{{m}}^2$")


# Plot Stress on Subplot 1
plt.subplot(1,2,2)

plt.plot(data_2021Soares['PStrain'], data_2021Soares['HRate'], color='blue', linestyle="--", linewidth=2, label="Soares et al. 2021")
#plt.plot(data_2003Uneishi['PStrain'], data_2003Uneishi['HRate'], color='magenta', linestyle="-.", linewidth=2, label="Uneishi et al. 2003")
# plt.plot(data_2006Kuroda['PStrain'], data_2006Kuroda['HRate'], color='green', linestyle=":", linewidth=2, label="Kuroda et al. 2006")



if args.fem_csv != "NULL":
    fem_data1['slope'] = fem_data1['stress_xx'].diff()/fem_data1['strain_xx'].diff()
    rho_total = (fem_data1['rho_i'][1] + fem_data1['rho_m'][1])*(24/1e8)
    plt.plot(fem_data1['Ep_eff'], fem_data1['slope'], color='red', linestyle=":", 
             linewidth=2, label=rf"DDCP $\rho_t = {rho_total:0.2f} \times 10^{{14}}\ \mathrm{{m}}^2$")



# plt.plot(dd_data_c1[:,0], (dd_data_c1[:,2] - dd_data_c1[:,-1])/1e12, color='red',     linestyle="--", linewidth=1, label=args.sheet_name+"_C1" )


# Final Plot Styling
fsize=16
plt.gcf().set_size_inches(16, 6)




plt.subplot(1,2,1)
plt.ylabel(r"$\mathrm{Stress\ }$, MPa", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.1)
plt.ylim(0, 700)
plt.legend(fontsize=fsize, loc="upper right", frameon=True)
plt.title(r"(a) Stress vs $\varepsilon_p$", fontsize=fsize)
plt.tick_params(labelsize=fsize)



plt.subplot(1,2,2)
plt.ylabel(r"$\mathrm{Hardening\ Rate}$, MPa", fontsize=fsize)
plt.xlabel(r"$\varepsilon_p$", fontsize=fsize)
plt.xlim(0, 0.10)
plt.ylim(0, 5000)
plt.legend(fontsize=fsize, loc="upper right", frameon=True)
plt.title(r"(b) Hardening Rate vs $\varepsilon_p$", fontsize=fsize)
plt.tick_params(labelsize=fsize)

if args.exportfig == True:
    plt.savefig("300K_HRate_SR1", bbox_inches="tight")


plt.show()
