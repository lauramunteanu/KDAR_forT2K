import uproot as up
import csv
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

def Print(string):
    print("\033[92m[OUTPUT]\033[0m :: ", string)

try:
    from numpy.linalg import inv
except ImportError:
    Err("NumPy matrix inversion not available. Cannot run chi2 calculation")


def chi2_calc(evis, xsec, covar, pred):
    # Need to apply a mask to avoid low/zero count bins,
    # select Evis between 50 & 135 MeV
    mask = [x >= 45 and x < 125 for x in evis]
    nbins = sum(mask)
    diffs = [(x-p) for x, p, m in zip(xsec, pred, mask) if m]
    covar_temp = [[x for x, m2 in zip(row, mask) if m2] for row, m1 in zip(covar, mask) if m1]
    covar_inv = inv(covar_temp)
    chi2 = sum([d1*covar_inv[i, j]*d2 for i, d1 in enumerate(diffs) for j, d2 in enumerate(diffs)])
    return chi2, nbins

def read_in_covar_matrix(filename):
    covar = [[None]*26 for _ in range(26)]
    for i, j, v in list(csv.reader(open(filename)))[1:]:
        i, j, v = int(i), int(j), float(v)
        covar[i][j] = v
    return covar


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

    return mc_contents
    




# ---------------------------------- Set up axes ----------------------------------

fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)

ax_main = fig.add_subplot(gs[0])
ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)

format_axis(ax_main, ax_ratio, r"$E_{m}$", "MeV", r"Normalised $\frac{1}{\sigma}\frac{\text{d}\sigma}{\text{d}E_{m}}$", r"1 / MeV neutron", "")


# ---------------------------------- Main plot ----------------------------------

plt.sca(ax_main)

data_centers, data_contents, data_errors, bin_edges = plot_data_hist("kdar_analysis/kdar_emiss_hist.root")

mc_samples = ["neutout/nuisflat_neut_5.4.0_JSPS_C.root",
              "neutout/nuisflat_neut_5.4.0_JSPS_C_noFSI.root",
              "neutout/EDRMF/nuisflat_NEUT6_EDRMF_nominal_noFSI.root",
              "neutout/EDRMF/nuisflat_NEUT6_EDRMF_nominal.root",]

mc_names   = ["NEUT", "NEUT noFSI", "NEUT6.0.2 EDRMF noFSI", "NEUT6.0.2 EDRMF",]
mc_colors  = ["red", "blue", "orange", "purple"]#,"green", "lime", "pink", "green", "brown"]

mc_counts = []
for index, sample in enumerate(mc_samples):
    Log(f"Reading in: {sample}")
    mc_counts_temp = makeMCplot(ax_ratio, sample, mc_names[index], mc_colors[index], data_contents, bin_edges)
    mc_counts.append(mc_counts_temp)

mc_counts = np.array(mc_counts)
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


# Chi2 calculation
xsec_vals = list(csv.reader(open("kdar_analysis/xsec.csv")))[1:]
xsec_values = [float(x[2]) for x in xsec_vals]
evis_vals = [float(x[1]) for x in xsec_vals]
keys = ["stat", "energy", "birks", "generator"]
covar_mats = {k: read_in_covar_matrix("kdar_analysis/covar_%s.csv" % k) for k in keys}
total_covar = [[sum([covar_mats[k][i][j] for k in keys]) for j in range(26)] for i in range(26)]

for index, sample_count in enumerate(mc_counts):
    chi2, ndof = chi2_calc(evis_vals, xsec_values, total_covar, list([float(x) for x in sample_count]))
    chi2_string = r"$\chi^{2} / N_{\text{dof}} "
    Print(chi2_string + f"for {mc_names[index]} = {chi2:.2f}/{ndof}")


plt.show()
