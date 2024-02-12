/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
double title_x = 0.16;
double title_y = 0.8;
double text_size = 0.045;
double label_font_size = 0.045;
double title_font_size = 0.045;
double marker_size = 1.75;
double pad_x0 = 0.f;
double pad_x1 = 1.f;
double pad_y0 = 0.f;
double pad_y1 = 1.f;

}  // namespace

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

std::vector<std::string> create_columns(const std::string& tag) {
    std::vector<std::string> varI = {"dl0", "dl1", "dphi", "dtheta", "dqop"};

    std::vector<std::string> varF = {"dl0", "dl1", "dphi", "dtheta", "dqop"};

    std::vector<std::string> columns;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            columns.push_back(varF[i] + varI[j] + "_" + tag);
        }
    }

    return columns;
}

std::array<double, 25u> get_means(ROOT::RDataFrame& rdf) {

    std::array<double, 25u> ret;

    auto col_names_R = create_columns("R");

    for (std::size_t i = 0u; i < 25u; i++) {
        const double residual_mean = *rdf.Mean<double>(col_names_R[i]);
        ret[i] = residual_mean;
    }

    return ret;
}

std::map<std::string, std::vector<double>> get_means(
    const std::vector<std::string> labels, const std::string tag, const int min,
    const int max, std::vector<double>& mean_step_sizes) {

    std::map<std::string, std::vector<double>> ret;

    int num = min;
    while (num <= max + 1e-3) {

        const std::string name = tag + "_" + std::to_string(num);
        const std::string csv_name = name + ".csv";

        const std::string sign = num >= 0 ? "p" : "m";
        const std::string root_name = name + ".root";

        std::cout << "Processing file: " << csv_name << std::endl;

        auto rdf = ROOT::RDF::FromCSV(csv_name);

        const std::array<double, 25u> means = get_means(rdf);

        for (unsigned int i = 0; i < 25u; i++) {
            ret[labels[i]].push_back(TMath::Log10(means[i]));
        }

        mean_step_sizes.push_back(
            TMath::Log10(*rdf.Mean<double>("average_step_size")));

        // Create root file
        rdf.Snapshot(tag.c_str(), root_name.c_str());

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

    TPad* gr_pad = new TPad("gr_pad", "gr_pad", 0, 0, 1, 1);
    gr_pad->Draw();
    gr_pad->cd();
    gr_pad->SetLeftMargin(80. / gr_pad->GetWw());
    gr_pad->SetTopMargin(50. / gr_pad->GetWh());
    gr_pad->SetBottomMargin(90. / gr_pad->GetWh());

    TGraph* gr[25];
    TMultiGraph* mg = new TMultiGraph();

    const std::array<int, 5u> marker_styles = {7, 2, 5, 27, 32};
    const std::array<int, 5u> line_styles = {1, 3, 2, 7, 4};
    const std::array<int, 5u> hues = {kOrange + 2, kPink + 5, kBlue + 2,
                                      kCyan + 2, kGreen + 2};

    auto legend = new TLegend(0.15, 0.63, 0.87, 0.945);
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
        gr[i]->SetMarkerSize(marker_size);
        gr[i]->SetLineStyle(line_styles[m]);
        gr[i]->SetMarkerColor(hues[m]);
        gr[i]->SetLineColor(hues[m]);

        mg->Add(gr[i]);
        legend->AddEntry(gr[i], labels[i].c_str(), "lp");
    }

    mg->GetXaxis()->SetTitle("log_{10}(#tau [mm])");
    mg->GetXaxis()->SetLabelOffset(-0.005);
    mg->GetXaxis()->SetTitleOffset(0.8);
    mg->GetYaxis()->SetTitle(
        "log_{10}(Mean of absolute and relative residuals)");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetYaxis()->SetRangeUser(-10, 8);
    mg->Draw("APL");

    TLegendEntry* header =
        (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextFont(22);
    header->SetTextSize(.03);

    legend->Draw();
}

void draw_mean_step_size(const std::string header_title,
                         const std::vector<double>& x_vec,
                         const std::vector<double>& means) {
    TPad* step_pad = new TPad("step_pad", "step_pad", 0, 0, 1, 1);
    step_pad->Draw();
    step_pad->cd();

    step_pad->SetLeftMargin(80. / step_pad->GetWw());
    step_pad->SetBottomMargin(60. / step_pad->GetWh());

    TGraph* gr = new TGraph(x_vec.size(), &x_vec[0], &means[0]);

    gr->SetMarkerStyle(4);
    gr->GetXaxis()->SetTitle("log_{10}(#tau [mm])");
    gr->GetYaxis()->SetTitle("log_{10}(Mean of average step sizes [mm])");
    gr->GetXaxis()->SetLimits(x_vec.front() - 0.5, x_vec.back() + 0.5);
    gr->GetYaxis()->SetRangeUser(means.front() - 0.2, means.back() + 0.95);
    gr->GetXaxis()->SetLabelSize(label_font_size);
    gr->GetYaxis()->SetLabelSize(label_font_size);
    gr->GetXaxis()->SetTitleSize(title_font_size);
    gr->GetYaxis()->SetTitleSize(title_font_size);
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetTitleOffset(1.1);

    gr->Draw();

    TPad* text_pad =
        new TPad("text_pad", "text_pad", pad_x0, pad_y0, pad_x1, pad_y1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    TLatex* ttext = new TLatex(title_x, title_y, header_title.c_str());
    ttext->SetTextFont(22);
    ttext->SetTextSize(text_size);
    ttext->Draw();
}

// ROOT Script for jacboain file reading
void rk_tolerance_comparison(int min, int max) {
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(0.02625);

    const std::array<float, 2> cdim1{800, 1350};
    const std::array<float, 2> cdim2{700, 500};

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
        new TCanvas("rect_canvas", "rect_canvas", cdim1[0], cdim1[1]);

    std::vector<double> rect_mean_step_sizes;
    const auto rect_y_means = get_means(labels, "inhom_rect_material", min, max,
                                        rect_mean_step_sizes);
    draw_graphs(rect_header, labels, x_vec, rect_y_means);

    rect_canvas->SaveAs(rect_pdf.c_str());

    auto rect_canvas2 =
        new TCanvas("rect_canvas2", "rect_canvas2", cdim2[0], cdim2[1]);
    const std::string rect_mean_step_pdf = "bound_to_bound_mean_step_size.pdf";

    draw_mean_step_size(rect_header, x_vec, rect_mean_step_sizes);
    rect_canvas2->SaveAs(rect_mean_step_pdf.c_str());

    /************************
     *  Wire
     * **********************/

    const std::string wire_header =
        "#splitline{Perigee-to-perigee transport,}{RKN with an inhomogeneous "
        "field "
        "and a material}";
    const std::string wire_pdf = "perigee_to_perigee_rk_tolerance.pdf";

    auto wire_canvas =
        new TCanvas("wire_canvas", "wire_canvas", cdim1[0], cdim1[1]);
    std::vector<double> wire_mean_step_sizes;
    const auto wire_y_means = get_means(labels, "inhom_wire_material", min, max,
                                        wire_mean_step_sizes);
    draw_graphs(wire_header, labels, x_vec, wire_y_means);

    wire_canvas->SaveAs(wire_pdf.c_str());

    auto wire_canvas2 =
        new TCanvas("wire_canvas2", "wire_canvas2", cdim2[0], cdim2[1]);
    const std::string wire_mean_step_pdf =
        "perigee_to_perigee_mean_step_size.pdf";

    draw_mean_step_size(wire_header, x_vec, wire_mean_step_sizes);
    wire_canvas2->SaveAs(wire_mean_step_pdf.c_str());
}