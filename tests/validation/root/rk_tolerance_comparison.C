/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <array>
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

std::array<double, 25u> get_means(ROOT::RDataFrame& rdf) {

    std::array<double, 25u> ret;

    auto col_names = rdf.GetColumnNames();

    // For absolute reltaive residuals
    auto cols_residual = parse_columns(col_names, "R");

    std::size_t i = 0u;

    for (const auto& col : cols_residual) {
        const double residual_mean = *rdf.Mean<double>(col);
        ret[i++] = residual_mean;
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

std::map<std::string, std::vector<double>> get_means(
    const std::vector<std::string> labels, const std::string tag, const int min,
    const int max) {

    std::map<std::string, std::vector<double>> ret;

    int num = min;
    while (num <= max + 1e-3) {

        const std::string name = tag + "_" + std::to_string(num);
        const std::string csv_name = name + ".csv";

        const std::string sign = num >= 0 ? "p" : "m";
        const std::string root_name = name + ".root";

        std::cout << "Processing file: " << csv_name << std::endl;

        auto rdf = ROOT::RDF::MakeCsvDataFrame(csv_name);

        const std::array<double, 25u> means = get_means(rdf);

        for (unsigned int i = 0; i < 25u; i++) {
            ret[labels[i]].push_back(TMath::Log10(means[i]));
        }

        // Create root file
        rdf.Snapshot(tag, root_name);

        num = num + 2;
    }

    return ret;
}

std::vector<double> get_x_vector(const int min, const int max) {
    std::vector<double> ret;

    int num = min;
    while (num <= max + 1e-3) {
        ret.push_back(num);
        num = num + 2;
    }

    return ret;
}

void draw_graphs(const std::string header_title,
                 const std::vector<std::string> labels,
                 const std::vector<double> x_vec,
                 std::map<std::string, std::vector<double>> means) {

    TGraph* gr[25];
    TMultiGraph* mg = new TMultiGraph();

    const std::array<int, 5u> marker_styles = {7, 2, 5, 27, 32};
    const std::array<int, 5u> line_styles = {1, 3, 2, 7, 4};
    const std::array<int, 5u> hues = {kOrange + 1, kPink + 7, kBlue + 2,
                                      kCyan + 1, kGreen + 1};

    auto legend = new TLegend(0.15, 0.595, 0.87, 0.88);
    legend->SetHeader(header_title.c_str());
    legend->SetNColumns(5);
    legend->SetColumnSeparation(-0.2);
    legend->SetFillStyle(0);
    legend->SetMargin(0.6);

    for (int i = 0; i < 25; i++) {
        gr[i] = new TGraph(x_vec.size(), &x_vec[0], &means[labels[i]][0]);

        const int n = i / 5;
        const int m = i % 5;

        gr[i]->SetMarkerStyle(marker_styles[n]);
        gr[i]->SetMarkerSize(1.4);
        gr[i]->SetLineStyle(line_styles[m]);
        gr[i]->SetMarkerColor(hues[m]);

        mg->Add(gr[i]);
        legend->AddEntry(gr[i], labels[i].c_str(), "lp");
    }

    mg->GetXaxis()->SetTitle("log_{10}(#tau [mm])");
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitle("Mean log_{10}(Absolute, relative residual)");
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetRangeUser(-9, 9);
    mg->Draw("APL");

    TLegendEntry* header =
        (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextFont(22);
    header->SetTextSize(.03);

    legend->Draw();
}

// ROOT Script for jacboain file reading
void rk_tolerance_comparison(int min, int max) {
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(0.023);

    const std::array<float, 2> cdim{800, 1400};

    auto labels = create_labels();

    const auto x_vec = get_x_vector(min, max);

    /************************
     *  Rectangular
     * **********************/

    const std::string rect_header =
        "#splitline{Bound-to-bound transport,}{RKN with an inhomogeneous field "
        "and a material}";
    const std::string rect_pdf = "bound_to_bound_rk_tolerance.pdf";

    auto rect_canvas =
        new TCanvas("rect_canvas", "rect_canvas", cdim[0], cdim[1]);
    const auto rect_y_means =
        get_means(labels, "inhom_rect_material", min, max);
    draw_graphs(rect_header, labels, x_vec, rect_y_means);

    rect_canvas->SaveAs(rect_pdf.c_str());

    /************************
     *  Wire
     * **********************/

    const std::string wire_header =
        "#splitline{Perigee-to-perigee transport,}{RKN with an inhomogeneous "
        "field "
        "and a material}";
    const std::string wire_pdf = "perigee_to_perigee_rk_tolerance.pdf";

    auto wire_canvas =
        new TCanvas("wire_canvas", "wire_canvas", cdim[0], cdim[1]);
    const auto wire_y_means =
        get_means(labels, "inhom_wire_material", min, max);
    draw_graphs(wire_header, labels, x_vec, wire_y_means);

    wire_canvas->SaveAs(wire_pdf.c_str());
}