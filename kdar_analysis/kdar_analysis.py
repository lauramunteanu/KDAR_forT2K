# -*- coding: utf-8 -*-
import csv
from array import array

import sys
import math
from ROOT import TFile
try:
    from numpy.linalg import inv
except ImportError:
    print("NumPy matrix inversion not available. Cannot run chi2 calculation")
try:
    import ROOT
except ImportError:
    print("Could not import ROOT. No plots will be produced")

BIN_WIDTH = 5.0  # Histogram bin widths [MeV]


#convert graph to th1 from chatgpt
def tgrapherrors_to_th1(graph, hist_name="hist_from_graph"):
    """
    Converts a TGraphErrors object into a TH1F histogram, preserving Y errors.

    Parameters:
    - graph: ROOT.TGraphErrors
    - hist_name: name for the resulting TH1F histogram

    Returns:
    - ROOT.TH1F object
    """
    n_points = graph.GetN()
    x_vals = graph.GetX()
    y_vals = graph.GetY()
    ex_vals = graph.GetEX()
    ey_vals = graph.GetEY()

    # Determine bin edges using x ± ex (for safety, if ex=0, use midpoint between points)
    bin_edges = []
    for i in range(n_points):
        x = x_vals[i]
        ex = ex_vals[i] if ex_vals[i] > 0 else 0.0
        left_edge = x - ex
        right_edge = x + ex

        if i == 0:
            bin_edges.append(left_edge)
        bin_edges.append(right_edge)

    # Remove duplicates and sort
    bin_edges = sorted(set(bin_edges))
    bin_edges_array = array('d', bin_edges)

    # Create histogram
    hist = ROOT.TH1F(hist_name, graph.GetTitle(), len(bin_edges_array) - 1, bin_edges_array)

    # Fill histogram with y values and set errors
    for i in range(n_points):
        x = x_vals[i]
        y = y_vals[i]
        ey = ey_vals[i]

        bin = hist.FindBin(x)
        hist.SetBinContent(bin, y)
        hist.SetBinError(bin, ey)

    return hist

# Example function for converting between visible energy & missing energy variable.
# Provided here just for clarity
def Evis2Emiss(Evis):
     """ Simple function to convert input visible energy value to equivalent missing energy value """
     E_nu = 235.532        # KDAR neutrino energy [MeV]
     M_mu = 105.6584       # Muon mass [MeV/c^2]
     M_proton = 938.2721   # Proton mass [MeV/c^2]
     M_neutron = 939.5654  # Neutron mass [MeV/c^2]
     return (E_nu + M_neutron) - (M_mu + M_proton) - Evis


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


def produce_unfolded_result_plots(evis, xsec, covar, pred):
    yaxis_title = "#frac{1}{#sigma} #frac{d#sigma}{dE_{vis}} / 5 MeV"
    h_pred = ROOT.TH1D("h_prediction", "Prediction;E_{vis} [MeV];%s" % yaxis_title,
                       26, 20, 150)

    g_xsec = ROOT.TGraphErrors()
    g_xsec.GetXaxis().SetTitle("E_{vis} [MeV]")
    g_xsec.GetYaxis().SetTitle(yaxis_title)

    # For simplicity just take the uncorrelated/diagonal variances.
    # A more sophisticated procedure, such as MC sampling a multivariate
    # Gaussian, should be employed to get error bars that include correlations.
    errs = [math.sqrt(covar[i][i]) for i in range(len(covar))]
    for i, (x, v, e, p) in enumerate(zip(evis, xsec, errs, pred)):
        g_xsec.SetPoint(i, x, v)
        g_xsec.SetPointError(i, BIN_WIDTH/2, e)
        h_pred.SetBinContent(i+1, p)
    return g_xsec, h_pred


def produce_signal_counts_plot(ereco, rates):
    g_signal = ROOT.TGraphAsymmErrors()
    g_signal.GetXaxis().SetTitle("Reconstructed Energy [MeV]")
    g_signal.GetYaxis().SetTitle("Counts")
    for i, (x, (y, yl, yh)) in enumerate(zip(ereco, rates)):
        g_signal.SetPoint(i, x, y)
        g_signal.SetPointEXhigh(i, 2.5)
        g_signal.SetPointEXlow(i, 2.5)
        g_signal.SetPointEYhigh(i, yh - y)
        g_signal.SetPointEYlow(i, y - yl)
    return g_signal


def read_in_covar_matrix(filename):
    covar = [[None]*26 for _ in range(26)]
    for i, j, v in list(csv.reader(open(filename)))[1:]:
        i, j, v = int(i), int(j), float(v)
        covar[i][j] = v
    return covar


if __name__ == "__main__":
    xsec_vals = list(csv.reader(open("xsec.csv")))[1:]
    evis_vals = [float(x[1]) for x in xsec_vals]
    emiss_vals = [Evis2Emiss(x) for x in evis_vals]
    xsec_vals = [float(x[2]) for x in xsec_vals]
    nuwro_pred_vals = [float(x[2]) for x in list(csv.reader(open("nuwro_prediction.csv")))[1:]]
    
    keys = ["stat", "energy", "birks", "generator"]
    covar_mats = {k: read_in_covar_matrix("covar_%s.csv" % k) for k in keys}

    # Estimate the total covariance matrix by summing the individual error matrices.
    # Here we combine error matrices by simply adding the individual matrices together.
    # A better way to combine matrices would be to MC sample from each matrix
    # and combine the samples, then use the summed samples to produce a new
    # covariance matrix.
    # We use the simpler method for this example script.
    total_covar = [[sum([covar_mats[k][i][j] for k in keys]) for j in range(26)] for i in range(26)]

#
    if 'numpy' not in sys.modules and 'ROOT' not in sys.modules:
        print("Cannot do anything without either NumPy or ROOT, sorry")
        sys.exit()
    if 'numpy' in sys.modules:
        chi2, ndof = chi2_calc(evis_vals, xsec_vals, total_covar, nuwro_pred_vals)
        print("NuWro chi2 = %f (%i)" % (chi2, ndof))

    if 'ROOT' in sys.modules:
        g_xsec, h_pred = produce_unfolded_result_plots(emiss_vals, xsec_vals, total_covar, nuwro_pred_vals)
        #g_xsec, h_pred = produce_unfolded_result_plots(evis_vals, xsec_vals, total_covar, nuwro_pred_vals)

        # Get the background subtracted signal rate data from csv file
        signal_rates = list(csv.reader(open("observed_counts.csv")))[1:]
        signal_rates = [[float(x[1]), float(x[2]), float(x[3]), float(x[4])] for x in signal_rates]
        ereco = [x[0] for x in signal_rates]
        signal_rates = [x[1:] for x in signal_rates]

        g_signal = produce_signal_counts_plot(ereco, signal_rates)

        c1 = ROOT.TCanvas()


        #convert graph to histogram:
        g_xsec_hist = tgrapherrors_to_th1(g_xsec)
        g_xsec_hist.SetLineColor(2)

        g_xsec.Draw("APE")
        h_pred.Draw("sameHIST")
        g_xsec_hist.Draw("E1 SAME")
        c1.Print("unfolded_result.png")

        c2 = ROOT.TCanvas()
        g_signal.Draw("APE")
        c2.Print("signal_rates.png")

        g_xsec_hist.Print("all")
        print(g_xsec_hist.GetBinLowEdge(1))
        print(g_xsec_hist.GetBinLowEdge(26))

        

        hist_file = TFile("kdar_emiss_hist.root", "RECREATE")
        g_xsec_hist.Write("g_xsec")
        hist_file.Close()
