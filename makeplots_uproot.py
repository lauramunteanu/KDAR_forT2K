import uproot as up
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
plt.style.use(["science", "notebook", "grid"])
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})


def plot_flattree_diff_xsec(filename: str, label: str, color: str, bin_edges, linestyle = ''):
    if(linestyle == ''):
        print("no linestyle give, using default")
        linestyle = '-'
    infile = up.open(filename)
    mode = infile["FlatTree_VARS;1"]["cc"].array()
    pdg = infile["FlatTree_VARS;1"]["pdg"].array()
    Elep = 1000*infile["FlatTree_VARS;1"]["ELep"].array()
    Enu = 235.5
    Ep = 1000*infile["FlatTree_VARS;1"]["E"].array()
    Mp = 938.27 
    
    # XSec_scale_factor = max(infile["FlatTree_VARS;1"]["fScaleFactor"].array()) # Don't need this if we are normalising
    Emiss = []
    for index, val in enumerate(mode):
        if(mode[index] == 1):
            for particle in range(0, len(pdg[index])):
                if(pdg[index][particle] == 2212):
                    Tp = Ep[index][particle] - Mp
                    Emiss_val = Enu - Elep[index] - Tp
                    Emiss.append(Emiss_val)
                    # Emiss.append(1000*Emiss[index]) # This is the nuisflat Emiss def


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

    plt.hist(Emiss, bins=bin_edges, histtype='step', weights=weights, color=color,linewidth=1.5, label = label, linestyle=linestyle)
    plt.legend(loc = "upper right")

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



# ---------------------------------- Set up axes ----------------------------------

fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)

ax_main = fig.add_subplot(gs[0])
ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)

format_axis(ax_main, ax_ratio, r"$E_{m}$", "MeV", r"$\text{d}\sigma/\text{d}E_{m}$", r"cm$^{2}$ / MeV$^{2}$ proton", "")


# ---------------------------------- Main plot ----------------------------------

plt.sca(ax_main)

data_centers, data_contents, data_errors, bin_edges = plot_data_hist("kdar_analysis/kdar_emiss_hist.root")
mc_centers, mc_contents = plot_flattree_diff_xsec("neutout/nuisflat_neut_5.4.0_JSPS_C.root", "NEUT", "red", bin_edges)
mc_centers2, mc_contents2 = plot_flattree_diff_xsec("neutout/nuisflat_neut_5.4.0_JSPS_C_noFSI.root", "NEUT_noFSI", "blue", bin_edges)
mc_centers3, mc_contents3 = plot_flattree_diff_xsec("neutout/nuisflat_NEUT_6.0.2_JSPS_EDRMF_nominal_shells_Nocascade.root", "NEUT 6.0.2 EDRMF noFSI", "orange", bin_edges)
mc_centers4, mc_contents4 = plot_flattree_diff_xsec("neutout/nuisflat_NEUT_6.0.2_JSPS_EDRMF_nominal_shells_cascade.root", "NEUT 6.0.2 EDRMF", "purple", bin_edges)

# plot_flattree_diff_xsec("neutout/nuisflat_NEUT_6.0.2_SF.root", "Emiss", "NEUT 6.0.2 SF $M_{A}^{QE} = 1.03$", "red", "--")
# plot_flattree_diff_xsec("neutout/nuisflat_NEUT_6.0.2_SF_noFSI.root", "Emiss", "NEUT 6.0.2 SF noFSI $M_{A}^{QE} = 1.03$", "blue", "--")

# ---------------------------------- Ratio plot ----------------------------------

## First two data points are 0, causes issues for ratio plot. 
## Below is a way to give a small value to the first two data points 
## so we can get a ratio.
# data_contents[0] = 1E-3
# data_contents[1] = 1E-3

ratio = mc_contents / data_contents
ratio2 = mc_contents2 / data_contents
ratio3 = mc_contents3 / data_contents
ratio4 = mc_contents4 / data_contents
ratio_errors = data_errors / data_contents  # relative error on data


ax_ratio.step(data_centers, ratio, color='red', linestyle='-',where='mid')
ax_ratio.step(data_centers, ratio2, color='blue', linestyle='-', where='mid')
ax_ratio.step(data_centers, ratio3, color='orange', linestyle='-', where='mid')
ax_ratio.step(data_centers, ratio4, color='purple', linestyle='-', where='mid')
ax_ratio.axhline(1, color='black', linestyle='--')
ax_ratio.set_ylabel("MC/Data", fontsize=13)
# ax_ratio.set_ylim(0.5,1.5)


plt.setp(ax_main.get_xticklabels(), visible=False) 
plt.show()


##################################
# C12_Emiss_profile = np.loadtxt("rho_C12_neutron.txt")
# Em_n = C12_Emiss_profile[:,0]
# rho_n = C12_Emiss_profile[:,1]
# integral = np.trapz(rho_n, Em_n)
# rho_n_scaled = 2*2.488923844830424 * rho_n / integral

# plt.plot(Em_n, rho_n_scaled, color = 'red', label = r'Neutron total $\rho(E_{m})$')
# plt.xlim(right=111.166899)
###################################