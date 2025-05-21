import uproot as up
import warnings
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
plt.style.use(["science", "notebook", "grid"])
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})

def Log(string):
    print("\033[94m[LOG]\033[0m :: ", string)

def Warn(string):
    print("\033[93m[WARNING]\033[0m :: ", string)

def Err(string):
    print("\033[91m[ERROR]\033[0m :: ", string)

def plot_flattree_diff_xsec(filename: str, label: str, color: str, bin_edges, linestyle = ''):
    if(linestyle == ''):
        Log("No linestyle given. Using default")
        linestyle = '-'
    infile = up.open(filename)
    mode = infile["FlatTree_VARS;1"]["cc"].array()
    Emiss_v = infile["FlatTree_VARS;1"]["Emiss"].array()
    Elep = 1000*infile["FlatTree_VARS;1"]["ELep"].array()
    Enu = 235.5
    Mp = 938.27 

    Ep, pdf = [], []
    if("noFSI" in filename):
        Ep = 1000*infile["FlatTree_VARS;1"]["E_vert"].array()
        pdg = infile["FlatTree_VARS;1"]["pdg_vert"].array()

    else:
        Ep = 1000*infile["FlatTree_VARS;1"]["E"].array()
        pdg = infile["FlatTree_VARS;1"]["pdg"].array()

    
    # XSec_scale_factor = max(infile["FlatTree_VARS;1"]["fScaleFactor"].array()) # Don't need this if we are normalising
    Emiss = []
    for index, val in enumerate(mode):
        if(mode[index] == 1):
            for particle in range(0, len(pdg[index])):
                if(pdg[index][particle] == 2212):
                    Tp = Ep[index][particle] - Mp
                    Emiss_val = Enu - Elep[index] - Tp
                    Emiss.append(Emiss_val)
            # Emiss.append(1000*Emiss_v[index]) # This is the nuisflat Emiss def


    # Compute histogram and integral
    counts, _ = np.histogram(Emiss, bins=bin_edges)
    bin_widths = np.diff(bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    integral = np.sum(counts * bin_widths)

    # This is the area of the data histogram (target scale)
    desired_total = 2*2.488923844830424

    # Compute scale factor
    scale = desired_total / integral

    weights = np.full_like(Emiss, scale)

    plt.hist(Emiss, bins=bin_edges, histtype='step', weights=weights, color=color,linewidth=1.8, label = label, linestyle=linestyle)
    plt.legend(loc = "upper right", fontsize=15)

    scaled_counts, _ = np.histogram(Emiss, bins=bin_edges, weights=weights)
    return bin_centers, scaled_counts

def plot_data_hist(filename: str):
    with up.open(filename) as data:
        data_hist = data["g_xsec"]

    # bin_edges = data_hist.axis().edges()
    bin_edges = data_hist.axis().edges()
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = (bin_edges[1:] - bin_edges[:-1]) / 2

    bin_contents = data_hist.values()
    bin_errors = np.sqrt(data_hist.variances())  # assumes Poisson errors
    integral = np.sum(bin_contents * bin_widths)

    # Plot with error bars
    plt.errorbar(
        bin_centers, bin_contents,
        xerr=bin_widths / 2, yerr=bin_errors,
        fmt='o', color='black', capsize=2, label="Data"
    )

    return bin_centers, bin_contents, bin_errors, bin_edges
    

def format_axis(ax, ax_ratio, x_title: str, xunits: str, y_title: str, yunits: str, title: str):
    ax_ratio.set_xlabel(f"{x_title} [" + xunits + "]", fontsize=20)
    ax.set_ylabel(f"{y_title} [" + yunits + "]", fontsize=20)
    ax.set_title(f"{title}", fontsize=20)

def makeMCplot(ax, MC_sample, MC_label, color, data_points, bin_edges):
    mc_centers, mc_contents = plot_flattree_diff_xsec(MC_sample, MC_label, color, bin_edges)
    
    ratio = -1
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ratio = mc_contents / data_points
        for warning in w:
            if "divide by zero encountered in divide" in str(warning.message):
                Warn("Divide by zero encountered in ratio plot")    
    ax.step(mc_centers, ratio, color=color, linestyle='-',where='mid')




# ---------------------------------- Set up axes ----------------------------------

fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)

ax_main = fig.add_subplot(gs[0])
ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)

format_axis(ax_main, ax_ratio, r"$E_{m}$", "MeV", r"Normalised $\frac{1}{\sigma}\frac{\text{d}\sigma}{\text{d}E_{m}}$", r"1 / MeV proton", "")


# ---------------------------------- Main plot ----------------------------------

plt.sca(ax_main)

data_centers, data_contents, data_errors, bin_edges = plot_data_hist("kdar_analysis/kdar_emiss_hist.root")

mc_samples = ["neutout/nuisflat_neut_5.4.0_JSPS_C.root", 
              "neutout/nuisflat_neut_5.4.0_JSPS_C_noFSI.root",
              "neutout/EDRMF/nuisflat_NEUT6_EDRMF_nominal_noFSI.root",
              "neutout/EDRMF/nuisflat_NEUT6_EDRMF_nominal.root",]
              # "neutout/EDRMF/nuisflat_NEUT_6.0.2_EDRMF_paperSFshells.root",]
            #   "neutout/RPWIA/nuisflat_NEUT_6.0.2_RPWIA_nominal_noFSI.root",
            #   "neutout/RPWIA/nuisflat_NEUT_6.0.2_RPWIA_nominal.root",
            #   "neutout/EDAIC/nuisflat_NEUT_6.0.2_EDAIC_nominal_noFSI.root",
            #   "neutout/EDAIC/nuisflat_NEUT_6.0.2_EDAIC_nominal.root",]
mc_names   = ["NEUT", "NEUT noFSI", "NEUT6.0.2 EDRMF noFSI", "NEUT6.0.2 EDRMF",]#"NEUT6.0.2 EDRMF paper shells", "NEUT6.0.2 RPWIA noFSI", "NEUT6.0.2 RPWIA", "NEUT6.0.2 EDAIC noFSI", "NEUT6.0.2 EDAIC"]
mc_colors  = ["red", "blue", "orange", "purple"]#,"green", "lime", "pink", "green", "brown"]

for index, sample in enumerate(mc_samples):
    Log(f"Reading in: {sample}")
    makeMCplot(ax_ratio, sample, mc_names[index], mc_colors[index], data_contents, bin_edges)

# ---------------------------------- Ratio plot ----------------------------------

## First two data points are 0, causes issues for ratio plot. 
## Below is a way to give a small value to the first two data points 
## so we can get a ratio.
# data_contents[0] = 1E-3
# data_contents[1] = 1E-3
ratio_errors = data_errors / data_contents  # relative error on data

ax_ratio.axhline(1, color='black', linestyle='--')
ax_ratio.set_ylabel("MC/Data", fontsize=13)
ax_ratio.set_ylim(0,2)


plt.setp(ax_main.get_xticklabels(), visible=False) 


##################################
# C12_Emiss_profile = np.loadtxt("rho_C12_neutron.txt")
# C12_Emiss_profile = np.loadtxt("rho_C12_neutron_NEUTsf_EmShells.txt")
# Em_n = C12_Emiss_profile[:,0]
# rho_n = C12_Emiss_profile[:,1]
# integral = np.trapz(rho_n, Em_n)
# rho_n_scaled = 2*2.488923844830424 * rho_n / integral

# ax_main.plot(Em_n, rho_n_scaled, color = 'red', label = r'Neutron total $\rho(E_{m})$', linestyle='dashed')
# ax_main.set_xlim(right=120)
###################################

plt.show()
