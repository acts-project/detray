/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveLabel.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <iostream>
#include <sstream>
#include <vector>

bool is_tagged(std::string col_name, const std::string& tag, const char delim) {
    std::istringstream iss(col_name);
    std::string str;
    std::vector<std::string> strings;

    while (getline(iss, str, delim)) {
        strings.push_back(str);
    }

    if (strings.back() == tag) {
        return true;
    }
    return false;
}

std::vector<std::string> parse_columns(
    const std::vector<std::string>& column_names, const std::string& tag,
    const char delim = '_') {

    std::vector<std::string> ret;

    for (const auto& col : column_names) {
        if (is_tagged(col, tag, delim)) {
            ret.push_back(col);
        }
    }
    return ret;
}

std::vector<double> get_means(ROOT::RDataFrame& rdf) {

    auto col_names = rdf.GetColumnNames();

    // For absolute reltaive residuals
    auto cols_residual = parse_columns(col_names, "R");

    /// For evaulated jacobians
    auto cols_evaluate = parse_columns(col_names, "E");

    /// For numerical differntitaions
    auto cols_numerics = parse_columns(col_names, "D");

    std::vector<std::tuple<std::string, std::string, std::string>> cols;

    for (std::size_t i = 0u; i < cols_residual.size(); i++) {
        cols.push_back(std::make_tuple(cols_residual[i], cols_evaluate[i],
                                       cols_numerics[i]));
    }
    std::vector<double> ret;

    for (auto& cols_per_variable : cols) {

        const auto col_residual = std::get<0>(cols_per_variable);
        const auto col_evaluate = std::get<1>(cols_per_variable);
        const auto col_numerics = std::get<2>(cols_per_variable);

        const double residual_mean = *rdf.Mean<double>(col_residual);
        const double evaluate_mean = *rdf.Mean<double>(col_evaluate);
        const double numerics_mean = *rdf.Mean<double>(col_numerics);

        std::cout << col_residual << "  " << residual_mean << "  "
                  << col_evaluate << "  " << evaluate_mean << "  "
                  << col_numerics << "  " << numerics_mean << std::endl;

        // If evalulated jacobian is too small, set the residual mean to
        // infinitesiaml value.
        if (residual_mean == 0.f || std::abs(evaluate_mean) < 1e-16 ||
            std::abs(numerics_mean) < 1e-16) {
            ret.push_back(1e-99);
        } else {
            ret.push_back(residual_mean);
        }
    }

    return ret;
}

std::vector<std::string> create_labels() {

    std::vector<std::string> varI = {
        "{#partiall_{0i}}", "{#partiall_{1i}}", "{#partial#phi_{i}}",
        "{#partial#theta_{i}}", "{#partial#lambda_{i}}"};

    std::vector<std::string> varF = {
        "{#partiall_{0f}}", "{#partiall_{1f}}", "{#partial#phi_{f}}",
        "{#partial#theta_{f}}", "{#partial#lambda_{f}}"};

    std::vector<std::string> labels;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            std::string head = "#frac";
            labels.push_back(head + varF[i] + varI[j]);
        }
    }

    return labels;
}

void fill_histo(TH1D* hist, const std::vector<double>& means,
                const std::vector<string>& labels, const int n_labels,
                const int marker_style) {
    hist->SetStats(0);
    hist->SetFillColor(38);
    hist->SetCanExtend(TH1::kAllAxes);

    for (int i = 0; i < n_labels; i++) {
        hist->Fill(labels[i].c_str(), TMath::Log10(means[i]));
    }

    hist->LabelsDeflate();
    hist->SetMarkerSize(1.35);
    hist->SetMarkerStyle(marker_style);
}

TH1D* get_histogram(std::string name, const int n_labels,
                    const int marker_style) {

    std::cout << "Generating histogram for " << name << std::endl;

    auto labels = create_labels();

    const std::string csv_name = name + ".csv";
    const std::string root_name = name + ".root";
    const std::string histo_name = name + "_histo";

    auto rdf = ROOT::RDF::MakeCsvDataFrame(csv_name);
    auto rdf_means = get_means(rdf);
    TH1D* histo = new TH1D(histo_name.c_str(), histo_name.c_str(), 3, 0, 3);
    histo->GetYaxis()->SetRangeUser(-14, 0);
    histo->GetYaxis()->SetTitle("Mean log_{10}(Absolute, relative residual)");

    fill_histo(histo, rdf_means, labels, n_labels, marker_style);

    // No color
    histo->SetFillColor(kWhite);

    // Create root file
    rdf.Snapshot(name, root_name);

    return histo;
}

void draw_text(const std::string& text) {

    const float x1 = 1.2f;
    const float y1 = -1.8f;

    TText* ttext = new TText(0.f, 0.f, text.c_str());
    ttext->SetTextFont(132);
    ttext->SetTextSize(1.5);

    UInt_t w;
    UInt_t h;
    ttext->Modify();
    ttext->GetBoundingBox(w, h);

    TPaveLabel* plabel =
        new TPaveLabel(x1, y1, x1 + float(w) / gPad->GetWw(),
                       y1 + float(h) / gPad->GetWh(), text.c_str());
    plabel->SetTextFont(132);
    plabel->SetFillColor(kWhite);
    plabel->Draw();
}

// ROOT Script for jacboain file reading
void jacobian_comparison() {

    gStyle->SetOptTitle(0);
    const std::array<float, 2> cdim{1200, 600};
    const std::array<float, 4> ldim{0.591, 0.14, 0.889, 0.32};
    const std::string rect_text = "Bound-to-bound transport";
    const std::string wire_text = "Perigee-to-perigee transport";
    const std::string rect_pdf = "bound_to_bound_jacobian_comparison.pdf";
    const std::string wire_pdf = "perigee_to_perigee_jacobian_comparison.pdf";

    /************************
     *  Rectangular
     * **********************/

    auto rect_canvas =
        new TCanvas("rect_canvas", "rect_canvas", cdim[0], cdim[1]);
    rect_canvas->SetGridx();

    auto rect_legend = new TLegend(ldim[0], ldim[1], ldim[2], ldim[3]);
    rect_legend->SetMargin(0.2);

    auto inhom_rect_material_histo =
        get_histogram("inhom_rect_material", 25, 28);
    inhom_rect_material_histo->Draw("hist P ");
    rect_legend->AddEntry(inhom_rect_material_histo,
                          "RKN with an inhomogeneous field and a material",
                          "p");

    auto inhom_rect_histo = get_histogram("inhom_rect", 20, 26);
    inhom_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(inhom_rect_histo, "RKN with an inhomogeneous field",
                          "p");

    auto const_rect_histo = get_histogram("const_rect", 15, 25);
    const_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(const_rect_histo, "RKN with a homogeneous field",
                          "p");

    auto helix_rect_histo = get_histogram("helix_rect", 15, 24);
    helix_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(helix_rect_histo, "Helix with a homogeneous field",
                          "p");

    rect_legend->Draw();
    draw_text(rect_text);
    rect_canvas->Draw();

    rect_canvas->SaveAs(rect_pdf.c_str());

    /************************
     *  Wire
     * **********************/

    auto wire_canvas =
        new TCanvas("wire_canvas", "wire_canvas", cdim[0], cdim[1]);
    wire_canvas->SetGridx();

    auto wire_legend = new TLegend(ldim[0], ldim[1], ldim[2], ldim[3]);
    wire_legend->SetMargin(0.2);

    auto inhom_wire_material_histo =
        get_histogram("inhom_wire_material", 25, 28);
    inhom_wire_material_histo->Draw("hist P ");
    wire_legend->AddEntry(inhom_wire_material_histo,
                          "RKN with an inhomogeneous field and a material",
                          "p");

    auto inhom_wire_histo = get_histogram("inhom_wire", 20, 26);
    inhom_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(inhom_wire_histo, "RKN with an inhomogeneous field",
                          "p");

    auto const_wire_histo = get_histogram("const_wire", 15, 25);
    const_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(const_wire_histo, "RKN with a homogeneous field",
                          "p");

    auto helix_wire_histo = get_histogram("helix_wire", 15, 24);
    helix_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(helix_wire_histo, "Helix with a homogeneous field",
                          "p");

    wire_legend->Draw();
    draw_text(wire_text);
    wire_canvas->Draw();

    wire_canvas->SaveAs(wire_pdf.c_str());
}