/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <Math/ProbFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace {
double x_pos = 0.15f;
double title_x = x_pos;
double title_y = 0.8f;
double fit_title_x = x_pos;
double fit_title_y = 0.707f;
double gaus_fit_par_x = x_pos;
double gaus_fit_par_y = 0.64f;
double const_fit_par_x = x_pos;
double const_fit_par_y = 0.658f;
double tolerance_x = 0.7f;
double tolerance_y = 0.67f;
}  // namespace

auto get_tree(std::string name) {

    const std::string csv_name = name + ".csv";
    const std::string root_name = name + ".root";

    auto rdf = ROOT::RDF::MakeCsvDataFrame(csv_name);

    // Create root file
    rdf.Snapshot(name, root_name);

    auto f = TFile::Open(root_name.c_str(), "read");
    auto t = (TTree*)f->Get(name.c_str());

    return t;
}

std::pair<std::array<double, 3u>, std::array<double, 3u>> fit_pull(
    TH1D* h_pull) {

    // Function used for the fit.
    TF1 gaus{"gaus", "gaus", -5.f, 5.f};
    double fit_par[3];
    double fit_par_error[3];

    // Set the mean seed to 0
    gaus.SetParameters(1, 0);
    gaus.SetParLimits(1, -1., 1.);
    // Set the standard deviation seed to 1
    gaus.SetParameters(2, 1.0);
    gaus.SetParLimits(2, 0.5, 2.);

    auto res = h_pull->Fit("gaus", "Q0S");
    gaus.GetParameters(&fit_par[0]);

    std::array<double, 3u> par{fit_par[0], fit_par[1], fit_par[2]};
    std::array<double, 3u> error;
    error[0] = gaus.GetParError(0);
    error[1] = gaus.GetParError(1);
    error[2] = gaus.GetParError(2);

    return {par, error};
}

std::pair<double, double> fit_pval(TH1D* h_pval) {

    // Function used for the fit.
    TF1 unif{"uniform", "[0]", 0.f, 1.f};
    double fit_par[1];

    auto res = h_pval->Fit("uniform", "Q0S");
    unif.GetParameters(&fit_par[0]);
    double error = unif.GetParError(0);

    return {fit_par[0], error};
}

void set_yaxis_title(TH1D* h) {
    double bin_width = h->GetBinWidth(0u);
    std::string str = std::to_string(bin_width);
    str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    str.erase(str.find_last_not_of('.') + 1, std::string::npos);
    std::string y_axis_title = "Counts / (" + str + ")";
    h->GetYaxis()->SetTitle(y_axis_title.c_str());
    h->GetYaxis()->SetTitleSize(0.04);
}

void set_xaxis_title(TH1D* h) {

    std::string x_axis_title;

    const TString h_name = h->GetName();

    if (h_name.Contains("l0")) {
        x_axis_title = "l_{0} pull";
    } else if (h_name.Contains("l1")) {
        x_axis_title = "l_{1} pull";
    } else if (h_name.Contains("phi")) {
        x_axis_title = "#phi pull";
    } else if (h_name.Contains("theta")) {
        x_axis_title = "#theta pull";
    } else if (h_name.Contains("qop")) {
        x_axis_title = "#lambda pull";
    } else if (h_name.Contains("pval")) {
        x_axis_title = "p value";
    }

    h->GetXaxis()->SetTitle(x_axis_title.c_str());
    h->GetXaxis()->SetTitleSize(0.04);
}

void draw_title(const std::string& text, const double x, const double y) {

    TLatex* ttext = new TLatex(x, y, text.c_str());
    ttext->SetTextFont(22);
    ttext->SetTextSize(0.04);
    ttext->Draw();
}

void draw_fit_title(const std::string title, const double x, const double y) {

    TLatex* ttext = new TLatex(x, y, title.c_str());
    ttext->SetTextFont(22);
    ttext->SetTextSize(0.04);
    ttext->Draw();
    gPad->cd();
}

void draw_gaus_fit_par(const std::array<double, 3u>& fit_par,
                       const std::array<double, 3u>& fit_par_error,
                       const double x, const double y) {
    TLatex* ttext = new TLatex(x, y, "#splitline{Mean}{Sigma}");
    ttext->SetTextFont(132);
    ttext->SetTextSize(0.04);
    ttext->Draw();

    std::stringstream mean_stream;
    mean_stream << std::fixed << std::setprecision(3) << fit_par[1] << " #pm "
                << fit_par_error[1];
    std::stringstream sigma_stream;
    sigma_stream << std::fixed << std::setprecision(3) << fit_par[2] << " #pm "
                 << fit_par_error[2];

    TLatex* ttext2 = new TLatex(x + 0.1, y,
                                "#splitline{" + TString(mean_stream.str()) +
                                    "}{" + TString(sigma_stream.str()) + "}");
    ttext2->SetTextFont(132);
    ttext2->SetTextSize(0.04);
    ttext2->Draw();
}

void draw_const_fit_par(const double fit_par, const double fit_par_error,
                        const double x, const double y) {
    std::stringstream val_stream;
    val_stream << std::fixed << std::setprecision(3) << fit_par << " #pm "
               << fit_par_error;

    TLatex* ttext = new TLatex(x, y, "Value  " + TString(val_stream.str()));
    ttext->SetTextFont(132);
    ttext->SetTextSize(0.04);
    ttext->Draw();
}

void draw_tolerance(const double log10_rk_tolerance, const double x,
                    const double y) {
    std::stringstream val_stream;
    val_stream << "#tau = 10^{" << int(log10_rk_tolerance) << "}";

    TLatex* ttext = new TLatex(x, y, TString(val_stream.str()));
    ttext->SetTextFont(132);
    ttext->SetTextSize(0.052);
    ttext->Draw();
}

void draw_pull(TH1D* h_pull, const std::string& title_text,
               const double log10_rk_tol) {

    auto fit_res = fit_pull(h_pull);
    auto fit_pars = fit_res.first;
    auto fit_errors = fit_res.second;

    set_xaxis_title(h_pull);
    set_yaxis_title(h_pull);
    const double y_axis_max = h_pull->GetEntries() * 10.f;
    h_pull->GetYaxis()->SetRangeUser(1.f, y_axis_max);
    h_pull->GetYaxis()->SetMaxDigits(1);
    h_pull->Draw();
    TF1* gaus = new TF1(h_pull->GetName(), "gaus", -5, 5);
    gaus->SetParameters(fit_pars[0], fit_pars[1], fit_pars[2]);
    gaus->Draw("same");

    TPad* text_pad = new TPad("gaus_text_pad", "gaus_text_pad", 0, 0, 1, 1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    draw_title(title_text.c_str(), title_x, title_y);
    draw_fit_title("Gaussian fit", fit_title_x, fit_title_y);
    draw_gaus_fit_par(fit_pars, fit_errors, gaus_fit_par_x, gaus_fit_par_y);
    // draw_tolerance(log10_rk_tol, tolerance_x, tolerance_y);
}

void draw_pval(TH1D* h_pval, const std::string& title_text,
               const double log10_rk_tol) {

    auto fit_res = fit_pval(h_pval);
    auto fit_par = fit_res.first;
    auto fit_error = fit_res.second;
    set_xaxis_title(h_pval);
    set_yaxis_title(h_pval);
    const double y_axis_max = 2.f * fit_par;
    h_pval->GetYaxis()->SetRangeUser(0.f, y_axis_max);
    h_pval->Draw();

    TF1* unif = new TF1(h_pval->GetName(), "[0]", -5, 5);
    unif->SetParameters(&fit_par);
    unif->Draw("same");

    TPad* text_pad = new TPad("const_text_pad", "const_text_pad", 0, 0, 1, 1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    draw_title(title_text.c_str(), title_x, title_y);
    draw_fit_title("Constant fit", fit_title_x, fit_title_y);
    draw_const_fit_par(fit_par, fit_error, const_fit_par_x, const_fit_par_y);
    // draw_tolerance(log10_rk_tol, tolerance_x, tolerance_y);
}

std::string to_pdf(const std::string& name) {
    return name + ".pdf";
}

void read_tree(TTree* t, const std::string& tag, const std::string& title) {
    const std::array<float, 2> cdim{600, 500};

    float pull_min = -5.f;
    float pull_max = 5.f;
    int n_bins = 100;

    double pull_l0;
    double pull_l1;
    double pull_phi;
    double pull_theta;
    double pull_qop;
    double chi2;
    double log10_rk_tolerance;

    t->SetBranchAddress("pull_l0", &pull_l0);
    t->SetBranchAddress("pull_l1", &pull_l1);
    t->SetBranchAddress("pull_phi", &pull_phi);
    t->SetBranchAddress("pull_theta", &pull_theta);
    t->SetBranchAddress("pull_qop", &pull_qop);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("log10_rk_tolerance", &log10_rk_tolerance);

    std::string l0_name = tag + "_pull_l0";
    std::string l1_name = tag + "_pull_l1";
    std::string phi_name = tag + "_pull_phi";
    std::string theta_name = tag + "_pull_theta";
    std::string qop_name = tag + "_pull_qop";
    std::string chi2_name = tag + "_chi2";
    std::string pval_name = tag + "_pval";

    TH1D* h_l0 =
        new TH1D(l0_name.c_str(), l0_name.c_str(), n_bins, pull_min, pull_max);
    TH1D* h_l1 =
        new TH1D(l1_name.c_str(), l1_name.c_str(), n_bins, pull_min, pull_max);
    TH1D* h_phi = new TH1D(phi_name.c_str(), phi_name.c_str(), n_bins, pull_min,
                           pull_max);
    TH1D* h_theta = new TH1D(theta_name.c_str(), theta_name.c_str(), n_bins,
                             pull_min, pull_max);
    TH1D* h_qop = new TH1D(qop_name.c_str(), qop_name.c_str(), n_bins, pull_min,
                           pull_max);
    TH1D* h_chi2 =
        new TH1D(chi2_name.c_str(), chi2_name.c_str(), n_bins, 0.f, 50.f);
    TH1D* h_pval = new TH1D(pval_name.c_str(), pval_name.c_str(), 50, 0.f, 1.f);

    // Fill the histograms
    for (int i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        h_l0->Fill(pull_l0);
        h_l1->Fill(pull_l1);
        h_phi->Fill(pull_phi);
        h_theta->Fill(pull_theta);
        h_qop->Fill(pull_qop);
        h_chi2->Fill(chi2);
        h_pval->Fill(ROOT::Math::chisquared_cdf_c(chi2, 5.f));
    }

    auto c_l0 = new TCanvas(h_l0->GetName(), h_l0->GetName(), cdim[0], cdim[1]);
    c_l0->SetLogy();
    draw_pull(h_l0, title, log10_rk_tolerance);

    c_l0->SaveAs(to_pdf(l0_name).c_str());

    auto c_l1 = new TCanvas(h_l1->GetName(), h_l1->GetName(), cdim[0], cdim[1]);
    c_l1->SetLogy();
    draw_pull(h_l1, title, log10_rk_tolerance);
    c_l1->SaveAs(to_pdf(l1_name).c_str());

    auto c_phi =
        new TCanvas(h_phi->GetName(), h_phi->GetName(), cdim[0], cdim[1]);
    c_phi->SetLogy();
    draw_pull(h_phi, title, log10_rk_tolerance);
    c_phi->SaveAs(to_pdf(phi_name).c_str());

    auto c_theta =
        new TCanvas(h_theta->GetName(), h_theta->GetName(), cdim[0], cdim[1]);
    c_theta->SetLogy();
    draw_pull(h_theta, title, log10_rk_tolerance);
    c_theta->SaveAs(to_pdf(theta_name).c_str());

    auto c_qop =
        new TCanvas(h_qop->GetName(), h_qop->GetName(), cdim[0], cdim[1]);
    c_qop->SetLogy();
    draw_pull(h_qop, title, log10_rk_tolerance);
    c_qop->SaveAs(to_pdf(qop_name).c_str());

    auto c_chi2 =
        new TCanvas(h_chi2->GetName(), h_chi2->GetName(), cdim[0], cdim[1]);
    h_chi2->Draw();
    c_chi2->SaveAs(to_pdf(chi2_name).c_str());

    auto c_pval =
        new TCanvas(h_pval->GetName(), h_pval->GetName(), cdim[0], cdim[1]);
    draw_pval(h_pval, title, log10_rk_tolerance);
    c_pval->SaveAs(to_pdf(pval_name).c_str());
}

// ROOT Script for covariance file reading
void covariance_validation() {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    /************************
     *  Rectangular
     * **********************/

    std::string rect_name = "rect_cov_transport";
    auto rect_tree = get_tree(rect_name);
    const std::string rect_title =
        "#splitline{Bound-to-bound transport,}{RKN with an inhomogeneous "
        "field and a material}";
    read_tree(rect_tree, "bound_to_bound", rect_title);

    /************************
     *  Wire
     * **********************/

    std::string wire_name = "wire_cov_transport";
    auto wire_tree = get_tree(wire_name);
    const std::string wire_title =
        "#splitline{Perigee-to-perigee transport,}{RKN with an inhomogeneous "
        "field and a material}";
    read_tree(wire_tree, "perigee_to_perigee", wire_title);
}