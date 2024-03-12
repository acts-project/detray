/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveLabel.h>
#include <TROOT.h>
#include <TStyle.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace {
double labelx_font_size = 0.06;
double labelx_offset = 0.006;
double labely_font_size = 0.06;
double title_font_size = 0.06;
double title_offset = 0.65;
double marker_size = 1.3875;
double legend_margin = 0.105;
int title_font = 132;
int label_font = 132;
int legend_font = 132;
double legend_font_size = 0.050;
double y_min = -14;
double y_max = 10;
double header_size = 0.05;
std::array<float, 4> ldim{0.508347, 0.605, 0.942404, 0.882353};

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

    // Residual columns
    auto col_names_R = create_columns("R");

    // Evaluate columns
    auto col_names_E = create_columns("E");

    std::array<double, 25u> ret;

    for (std::size_t i = 0u; i < 25u; i++) {

        const double residual_mean = *rdf.Mean<double>(col_names_R[i]);
        const double evaluate_mean = *rdf.Mean<double>(col_names_E[i]);

        // dtheta/d(l0,l1,phi,qop) is analytically zero for homogeneous field
        // dqop(dl0,dl1,dphi,dtheta) is analytically zero for free space without
        // material so those should not appear on the plots
        if (std::abs(evaluate_mean) < 1e-20f) {
            ret[i] = std::numeric_limits<double>::min();
        } else {
            ret[i] = residual_mean;
        }
    }

    return ret;
}

TString get_legend_header(const double log10_rk_tol) {
    std::stringstream val_stream;
    val_stream << "#tau = 10^{" << int(log10_rk_tol)
               << "} mm as an RKN error tolerance";

    return TString(val_stream.str());
}

void fill_histo(TH1D* hist, const std::array<double, 25u>& means,
                const std::vector<string>& labels, const int n_labels,
                const int marker_style) {
    hist->SetStats(0);
    hist->SetFillColor(38);
    hist->SetCanExtend(TH1::kAllAxes);

    for (int i = 0; i < n_labels; i++) {
        hist->Fill(labels[i].c_str(), TMath::Log10(means[i]));
    }

    hist->LabelsDeflate();
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerStyle(marker_style);
}

void draw_lines() {
    for (int i = 1; i < 25; i++) {
        if (i % 5 == 0) {
            TLine* line = new TLine(i, y_min, i, y_max);
            line->SetLineColor(kBlack);
            line->SetLineWidth(1);
            line->Draw();
        } else {
            TLine* line = new TLine(i, y_min, i, y_max);
            line->SetLineColor(kBlack);
            line->SetLineWidth(1);
            line->SetLineStyle(3);
            line->Draw();
        }
    }
}

TH1D* get_histogram(std::string name, const int n_labels,
                    const int marker_style, double& log10_rk_tol) {

    std::cout << "Generating histogram for " << name << std::endl;

    auto labels = create_labels();

    const std::string csv_name = name + ".csv";
    const std::string root_name = name + ".root";
    const std::string histo_name = name + "_histo";

    auto rdf = ROOT::RDF::FromCSV(csv_name);
    auto rdf_means = get_means(rdf);
    TH1D* histo = new TH1D(histo_name.c_str(), histo_name.c_str(), 3, 0, 3);
    histo->GetYaxis()->SetRangeUser(y_min, y_max);
    histo->GetYaxis()->SetTitle("log_{10}(Mean of #font[12]{#Omega_{R}})");
    histo->GetYaxis()->SetTitleOffset(0.68);
    histo->GetYaxis()->SetNdivisions(406, kFALSE);
    histo->GetYaxis()->SetTitleSize(title_font_size);
    histo->GetYaxis()->SetTitleFont(title_font);
    histo->GetYaxis()->SetLabelSize(labely_font_size);
    histo->GetXaxis()->SetLabelSize(labelx_font_size);
    histo->GetXaxis()->SetLabelOffset(labelx_offset);
    histo->GetXaxis()->SetLabelFont(label_font);
    histo->GetYaxis()->CenterTitle(true);

    fill_histo(histo, rdf_means, labels, n_labels, marker_style);

    // No color
    histo->SetFillColor(kWhite);

    // Create root file
    rdf.Snapshot(name, root_name);

    // Get log10(rk_tolerance)
    if (!TString(name).Contains("helix")) {
        log10_rk_tol = *rdf.Mean<double>("log10_rk_tolerance");
    }

    return histo;
}

void draw_pad(const std::string& pad_name) {
    TPad* apad = new TPad(pad_name.c_str(), pad_name.c_str(), 0, 0, 1, 1);
    apad->Draw();
    apad->cd();

    apad->SetLeftMargin(100. / apad->GetWw());
    apad->SetRightMargin(60. / apad->GetWw());
    apad->SetBottomMargin(60. / apad->GetWh());
}

void draw_text(const std::string& text) {

    const float x1 = 1.23;
    const float y1 = 6.05;

    TLatex* ttext = new TLatex(0.f, 0.f, text.c_str());
    ttext->SetTextFont(132);
    ttext->SetTextSize(3.5);

    UInt_t w;
    UInt_t h;
    ttext->GetBoundingBox(w, h);

    TPaveLabel* plabel =
        new TPaveLabel(x1, y1, x1 + float(w) / gPad->GetWw() * 0.62,
                       y1 + float(h) / gPad->GetWh() * 1.15, text.c_str());

    plabel->SetTextFont(132);
    plabel->SetFillColor(kWhite);
    plabel->Draw();
}

void draw_text2() {

    const double y = 7.5;

    TLatex* ttext1 = new TLatex(
        15.2, y, "#splitline{Measurable with}{an inhomgeneous field}");
    TLatex* ttext2 =
        new TLatex(20.2, y, "#splitline{Measurable with}{a material}");
    ttext1->SetTextFont(132);
    ttext1->SetTextSize(0.045);
    ttext2->SetTextFont(132);
    ttext2->SetTextSize(0.045);

    ttext1->Draw();
    ttext2->Draw();
}

// ROOT Script for jacboain file reading
void jacobian_comparison() {

    gStyle->SetOptTitle(0);
    const std::array<float, 2> cdim{1200, 500};
    const std::array<int, 4> markers{kOpenCross, kOpenTriangleUp, kOpenSquare,
                                     kOpenCircle};
    const std::array<int, 4> hues{kOrange + 8, kMagenta + 1, kAzure,
                                  kGreen + 2};
    const std::string rect_pdf = "bound_to_bound_jacobian_comparison.pdf";
    const std::string wire_pdf = "perigee_to_perigee_jacobian_comparison.pdf";

    double log10_rk_tolerance_rect;
    double log10_rk_tolerance_wire;
    double dummy;

    /************************
     *  Rectangular
     * **********************/

    auto rect_canvas =
        new TCanvas("rect_canvas", "rect_canvas", cdim[0], cdim[1]);
    draw_pad("rect_pad");

    auto rect_legend = new TLegend(ldim[0], ldim[1], ldim[2], ldim[3]);
    rect_legend->SetMargin(legend_margin);
    rect_legend->SetTextFont(legend_font);
    rect_legend->SetTextSize(legend_font_size);
    rect_legend->SetBorderSize(4);

    std::string rect_text = "Bound-to-bound transport";
    auto inhom_rect_material_histo = get_histogram(
        "inhom_rect_material", 25, markers[0u], log10_rk_tolerance_rect);
    // rect_legend->SetHeader("Bound-to-bound transport");
    inhom_rect_material_histo->SetMarkerColor(hues[0u]);
    inhom_rect_material_histo->Draw("hist P ");
    rect_legend->AddEntry(inhom_rect_material_histo,
                          "RKN with an inhomogeneous field and a material",
                          "p");

    auto inhom_rect_histo = get_histogram("inhom_rect", 20, markers[1u], dummy);
    inhom_rect_histo->SetMarkerColor(hues[1u]);
    inhom_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(inhom_rect_histo, "RKN with an inhomogeneous field",
                          "p");

    auto const_rect_histo = get_histogram("const_rect", 15, markers[2u], dummy);
    const_rect_histo->SetMarkerColor(hues[2u]);
    const_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(const_rect_histo, "RKN with a homogeneous field",
                          "p");

    auto helix_rect_histo = get_histogram("helix_rect", 15, markers[3u], dummy);
    helix_rect_histo->SetMarkerColor(hues[3u]);
    helix_rect_histo->Draw("hist P same");
    rect_legend->AddEntry(helix_rect_histo, "Helix with a homogeneous field",
                          "p");

    /*
    TLegendEntry* rect_header =
        (TLegendEntry*)rect_legend->GetListOfPrimitives()->First();
    rect_header->SetTextFont(22);
    rect_header->SetTextSize(header_size);
    */

    draw_lines();
    // draw_text2();
    rect_legend->Draw();
    draw_text(rect_text);
    rect_canvas->Draw();

    rect_canvas->SaveAs(rect_pdf.c_str());

    /************************
     *  Wire
     * **********************/

    auto wire_canvas =
        new TCanvas("wire_canvas", "wire_canvas", cdim[0], cdim[1]);
    draw_pad("wire_pad");

    auto wire_legend = new TLegend(ldim[0], ldim[1], ldim[2], ldim[3]);
    wire_legend->SetMargin(legend_margin);
    wire_legend->SetTextFont(legend_font);
    wire_legend->SetTextSize(legend_font_size);
    wire_legend->SetBorderSize(4);

    std::string wire_text = "Perigee-to-perigee transport";
    auto inhom_wire_material_histo = get_histogram(
        "inhom_wire_material", 25, markers[0u], log10_rk_tolerance_wire);
    // wire_legend->SetHeader("Perigee-to-perigee transport");
    inhom_wire_material_histo->SetMarkerColor(hues[0u]);
    inhom_wire_material_histo->Draw("hist P ");
    wire_legend->AddEntry(inhom_wire_material_histo,
                          "RKN with an inhomogeneous field and a material",
                          "p");

    auto inhom_wire_histo = get_histogram("inhom_wire", 20, markers[1u], dummy);
    inhom_wire_histo->SetMarkerColor(hues[1u]);
    inhom_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(inhom_wire_histo, "RKN with an inhomogeneous field",
                          "p");

    auto const_wire_histo = get_histogram("const_wire", 15, markers[2u], dummy);
    const_wire_histo->SetMarkerColor(hues[2u]);
    const_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(const_wire_histo, "RKN with a homogeneous field",
                          "p");

    auto helix_wire_histo = get_histogram("helix_wire", 15, markers[3u], dummy);
    helix_wire_histo->SetMarkerColor(hues[3u]);
    helix_wire_histo->Draw("hist P same");
    wire_legend->AddEntry(helix_wire_histo, "Helix with a homogeneous field",
                          "p");

    /*
    TLegendEntry* wire_header =
        (TLegendEntry*)wire_legend->GetListOfPrimitives()->First();
    wire_header->SetTextFont(22);
    wire_header->SetTextSize(header_size);
    */

    draw_lines();
    // draw_text2();
    wire_legend->Draw();
    draw_text(wire_text);
    wire_canvas->Draw();

    wire_canvas->SaveAs(wire_pdf.c_str());
}