import ROOT

def set_hist_style(hist, color):
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(20)
    hist.SetStats(0)

def set_canvas_style(canvas):
    canvas.SetGrid()
    canvas.SetTickx()
    canvas.SetTicky()
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)

def set_legend_style(legend):
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)

def draw_ratio_plot(_hist1, _hist2, canvas, legend, title):
    hist1 = _hist1.Clone("hist1")
    hist1.SetDirectory(0)
    hist2 = _hist2.Clone("hist2")
    hist2.SetDirectory(0)

    # add overflow bin to the last bin
    hist1.SetBinContent(hist1.GetNbinsX(), hist1.GetBinContent(hist1.GetNbinsX()) + hist1.GetBinContent(hist1.GetNbinsX() + 1))
    hist2.SetBinContent(hist2.GetNbinsX(), hist2.GetBinContent(hist2.GetNbinsX()) + hist2.GetBinContent(hist2.GetNbinsX() + 1))
    upper_pad = ROOT.TPad("upper_pad", "upper_pad", 0, 0.3, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "lower_pad", 0, 0, 1, 0.3)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
    upper_pad.Draw()
    lower_pad.Draw()

    upper_pad.cd()
    upper_pad.SetLogy()
    hist1.Draw("E")
    hist2.Draw("E SAME")
    legend.Draw()

    ratio = hist1.Clone("ratio")
    ratio.Divide(hist2)
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(20)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.55)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(15)
    ratio.GetXaxis().SetTitleSize(20)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(4.)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(15)

    # Draw line at 1.0

    lower_pad.cd()
    ratio.Draw("E")
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0, ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("same")

    canvas.SaveAs(title)

def compare_and_reweight():
    # Open the ROOT file
    # file_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_3/output_zmumu/zmumu_v18_reweight_sig/hists/Cut_014_has_more_two_jets_jet_lead_score_dxysig.root"
    file_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_3/output_zmumu/zmumu_v19_maxdxy_v2/hists/Cut_015_after_jet_def_jet_lead_score_maxdxy.root"
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print("Error opening file!")
        return

    # Retrieve the histograms
    h_SingleMuon_OS = root_file.Get("SingleMuon/iso/OS/hist")
    h_Inclusive_DYLO_M50_OS = root_file.Get("Inclusive_DYLO_M-50/iso/OS/hist")
    h_Inclusive_DYLO_M50_OS.Scale(6424000 * 59.7 / 76736130)

    # rebin the histograms
    h_SingleMuon_OS.Rebin2D(2, 2)
    h_Inclusive_DYLO_M50_OS.Rebin2D(2, 2)
    h_SingleMuon_OS_score_ratio = h_SingleMuon_OS.Clone("h_SingleMuon_OS_score_ratio")
    h_SingleMuon_OS_score_ratio.Divide(h_Inclusive_DYLO_M50_OS)
    c0 = ROOT.TCanvas("c1", "Tagger Score Comparison", 800, 800)
    c0.SetLogz()
    h_SingleMuon_OS_score_ratio.Draw("colz")
    c0.SaveAs("output/2D.png")
    # h_SingleMuon_OS.GetYaxis().SetRangeUser(0., 50)
    # h_Inclusive_DYLO_M50_OS.GetYaxis().SetRangeUser(0., 50)


    if not h_SingleMuon_OS or not h_Inclusive_DYLO_M50_OS:
        print("Error retrieving histograms!")
        root_file.Close()
        return

    # Project histograms onto X-axis (tagger score) and Y-axis (dxy)
    h_SingleMuon_OS_score = h_SingleMuon_OS.ProjectionX("h_SingleMuon_OS_score")
    h_SingleMuon_OS_dxy = h_SingleMuon_OS.ProjectionY("h_SingleMuon_OS_dxy")
    h_Inclusive_DYLO_M50_OS_score = h_Inclusive_DYLO_M50_OS.ProjectionX("h_Inclusive_DYLO_M50_OS_score")
    h_Inclusive_DYLO_M50_OS_dxy = h_Inclusive_DYLO_M50_OS.ProjectionY("h_Inclusive_DYLO_M50_OS_dxy")

    # Add 0.2*10^4 to every bin in h_Inclusive_DYLO_M50_OS_score
    # for i in range(0, h_Inclusive_DYLO_M50_OS_score.GetNbinsX() + 2):
    #     h_Inclusive_DYLO_M50_OS_score.SetBinContent(i, h_Inclusive_DYLO_M50_OS_score.GetBinContent(i) + 0.2 * 10 ** 4)

    # Reweight the Inclusive_DYLO_M50_OS to match the dxy distribution of SingleMuon_OS
    reweight_factor = h_SingleMuon_OS_dxy.Clone("reweight_factor")
    # reweight_factor.Print("all")
    reweight_factor.Divide(h_Inclusive_DYLO_M50_OS_dxy)
    # reweight_factor.Print("all")

    print(h_Inclusive_DYLO_M50_OS.GetNbinsX())
    print(h_Inclusive_DYLO_M50_OS.GetNbinsY())
    print(reweight_factor.GetNbinsX())

    h_SingleMuon_OS_dxy.Print("all")
    h_Inclusive_DYLO_M50_OS_dxy.Print("all")
    reweight_factor.Print("all")

    # print h_Inclusive_DYLO_M50_OS before reweighting
    # h_Inclusive_DYLO_M50_OS.Print("all")
    # Apply the reweighting to the tagger score distribution for Inclusive_DYLO_M50_OS including overflow bins
    h_Inclusive_DYLO_M50_OS_reweighted = h_Inclusive_DYLO_M50_OS.Clone("h_Inclusive_DYLO_M50_OS_reweighted")
    for i in range(0, h_Inclusive_DYLO_M50_OS_reweighted.GetNbinsX() + 2):
        for j in range(0, h_Inclusive_DYLO_M50_OS_reweighted.GetNbinsY() + 2):
            h_Inclusive_DYLO_M50_OS_reweighted.SetBinContent(i, j, h_Inclusive_DYLO_M50_OS_reweighted.GetBinContent(i, j) * reweight_factor.GetBinContent(j))
    # h_Inclusive_DYLO_M50_OS_reweighted.Print("all")

    # Set styles for histograms
    set_hist_style(h_SingleMuon_OS_score, ROOT.kRed)
    set_hist_style(h_SingleMuon_OS_dxy, ROOT.kRed)
    set_hist_style(h_Inclusive_DYLO_M50_OS_score, ROOT.kBlue)
    set_hist_style(h_Inclusive_DYLO_M50_OS_dxy, ROOT.kBlue)

    # Create legends
    legend_score = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    set_legend_style(legend_score)
    legend_score.AddEntry(h_SingleMuon_OS_score, "SingleMuon OS", "l")
    legend_score.AddEntry(h_Inclusive_DYLO_M50_OS_score, "Inclusive DYLO M-50 OS", "l")

    legend_dxy = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    set_legend_style(legend_dxy)
    legend_dxy.AddEntry(h_SingleMuon_OS_dxy, "SingleMuon OS", "l")
    legend_dxy.AddEntry(h_Inclusive_DYLO_M50_OS_dxy, "Inclusive DYLO M-50 OS", "l")

    # create folder output if not exists
    import os
    if not os.path.exists("output"):
        os.makedirs("output")

    # Create canvases to plot and compare distributions
    c1 = ROOT.TCanvas("c1", "Tagger Score Comparison", 800, 800)
    set_canvas_style(c1)
    draw_ratio_plot(h_SingleMuon_OS_score, h_Inclusive_DYLO_M50_OS_score, c1, legend_score, "output/score_comparison_with_ratio.png")

    c2 = ROOT.TCanvas("c2", "dxy Comparison", 800, 800)
    set_canvas_style(c2)
    draw_ratio_plot(h_SingleMuon_OS_dxy, h_Inclusive_DYLO_M50_OS_dxy, c2, legend_dxy, "output/dxy_comparison_with_ratio.png")

    # do the projection of score
    h_Inclusive_DYLO_M50_OS_score_reweighted = h_Inclusive_DYLO_M50_OS_reweighted.ProjectionX("h_Inclusive_DYLO_M50_OS_score_reweighted")

    # Set style for reweighted histogram
    set_hist_style(h_Inclusive_DYLO_M50_OS_score_reweighted, ROOT.kGreen)

    # Create legend for reweighted plot
    legend_reweighted = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    set_legend_style(legend_reweighted)
    legend_reweighted.AddEntry(h_SingleMuon_OS_score, "SingleMuon OS", "l")
    legend_reweighted.AddEntry(h_Inclusive_DYLO_M50_OS_score_reweighted, "Inclusive DYLO M-50 OS Reweighted", "l")

    # Plot and compare the reweighted score distributions
    c3 = ROOT.TCanvas("c3", "Reweighted Tagger Score Comparison", 800, 800)
    set_canvas_style(c3)
    draw_ratio_plot(h_SingleMuon_OS_score, h_Inclusive_DYLO_M50_OS_score_reweighted, c3, legend_reweighted, "output/score_comparison_with_ratio_reweighted.png")


    # do projection over dxy to make sure reweighting is correct
    h_Inclusive_DYLO_M50_OS_dxy_reweighted = h_Inclusive_DYLO_M50_OS_reweighted.ProjectionY("h_Inclusive_DYLO_M50_OS_dxy_reweighted")
    h_Inclusive_DYLO_M50_OS_dxy_reweighted.Print("all")
    # Set style for reweighted histogram
    set_hist_style(h_Inclusive_DYLO_M50_OS_dxy_reweighted, ROOT.kGreen)

    # Create legend for reweighted plot
    legend_reweighted_dxy = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    set_legend_style(legend_reweighted_dxy)
    legend_reweighted_dxy.AddEntry(h_SingleMuon_OS_dxy, "SingleMuon OS", "l")
    legend_reweighted_dxy.AddEntry(h_Inclusive_DYLO_M50_OS_dxy_reweighted, "Inclusive DYLO M-50 OS Reweighted", "l")

    # Plot and compare the reweighted dxy distributions
    c4 = ROOT.TCanvas("c4", "Reweighted dxy Comparison", 800, 800)
    set_canvas_style(c4)
    draw_ratio_plot(h_SingleMuon_OS_dxy, h_Inclusive_DYLO_M50_OS_dxy_reweighted, c4, legend_reweighted_dxy, "output/dxy_comparison_with_ratio_reweighted.png")


    # Close the ROOT file
    root_file.Close()

# Run the function
compare_and_reweight()