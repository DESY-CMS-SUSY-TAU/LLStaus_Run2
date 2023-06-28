import dataclasses
import numpy
import os
import sys
import ROOT

import utils.cms_lumi as CMS_lumi

def ColorIterator(index : int, scale : int) -> int:
    # kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    # kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    # kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
    cs = scale
    colow_wheel = [ cs + ROOT.kBlue,
                    cs + ROOT.kRed,
                    cs + ROOT.kGreen,
                    cs + ROOT.kMagenta,
                    cs + ROOT.kYellow,
                    cs + ROOT.kCyan ]
    color_idx =  colow_wheel +\
                [x - 5 for x in colow_wheel] +\
                [x - 10 for x in colow_wheel]
    return color_idx[index]

def get_canvas(ratio = False) :

    ROOT.gROOT.LoadMacro(os.path.split(os.path.realpath(__file__))[0]+"/tdrstyle.C")
    ROOT.gROOT.ProcessLine("setTDRStyle()")
    
    ROOT.gROOT.SetStyle("tdrStyle")
    ROOT.gROOT.ForceStyle(True)
    
    canvas = ROOT.TCanvas("canvas", "canvas", 1600, 1300)
    canvas.UseCurrentStyle()
    
    #canvas.SetLeftMargin(0.16)
    #canvas.SetRightMargin(0.05)
    #canvas.SetTopMargin(0.1)
    #canvas.SetBottomMargin(0.135)
    
    if (ratio) :
        
        canvas.Divide(1, 2)
        
        canvas.cd(1).SetPad(0, 0.32, 1, 1)
        canvas.cd(1).SetTopMargin(0.075)
        canvas.cd(1).SetBottomMargin(0)
        
        canvas.cd(2).SetPad(0, 0.0, 1, 0.3)
        canvas.cd(2).SetTopMargin(0.05)
        canvas.cd(2).SetBottomMargin(0.285)
    
    canvas.cd(1).SetLeftMargin(0.125)
    canvas.cd(1).SetRightMargin(0.05)
    
    if (ratio) :
        
        canvas.cd(2).SetLeftMargin(0.125)
        canvas.cd(2).SetRightMargin(0.05)
    
    return canvas

def get_min_max(hist):
    min_ =  10E+36
    max_ = -10E-36
    for i in range(1, hist.GetNbinsX(), 1):
        for j in range(1, hist.GetNbinsY(), 1):
            bin_z = hist.GetBinContent(i, j)
            if bin_z !=0:
                if bin_z > max_:
                    max_ = bin_z
                if bin_z < min_:
                    min_ = bin_z
    return min_, max_

def root_plot1D(
    l_hist,
    outfile,
    xrange, yrange,
    l_hist_overlay = [],
    logx = False, logy = False,
    include_overflow = False,
    title = "",
    xtitle = "", ytitle = "",
    xtitle_ratio = "", ytitle_ratio = "",
    yrange_ratio = (0.001, 1.1),
    logx_ratio = False, logy_ratio = True,
    centertitlex = True, centertitley = True,
    centerlabelx = False, centerlabely = False,
    gridx = False, gridy = False,
    ndivisionsx = None, ndivisionsy = None,
    ndivisionsy_ratio = (5, 5, 0),
    stackdrawopt = "nostack",
    normilize = False,
    normilize_overlay = True,
    legendpos = "UR",
    legendncol = 1,
    legendtextsize = 0.065,
    legendtitle = "",
    legendheightscale = 1.0,
    legendwidthscale = 1.0,
    ratio_num_den_pairs = [],
    ratio_mode = "B",
    signal_to_background_ratio = False,
    CMSextraText = "Private work (CMS simulation)",
    lumiText = "(13 TeV)",
    draw_errors = False
) :
    """
    l_hist: list of TH1 to be stacked according to `stackdrawopt`.
    l_hist_overlay: list of TH1 to be overlaid on the stack.
    stackdrawopt: pass empty string to stack
    ratio_num_den_pairs: list of (numerator TH1, denominator TH1) pairs of to be plotted as ratios: [(num1, den1), (num2, den2), ...].
    Note that the desired plotting styles and colors (like FillStyle/Color, LineSize/Style/Color, MarkerSize/Style/Color, SetOption) need to be set for the stack and overlay histograms before calling this function.
    """
    
    # canvas = get_canvas(ratio = len(ratio_num_den_pairs))
    canvas = get_canvas(ratio = signal_to_background_ratio)
    
    canvas.cd(1)
    
    
    legendHeight = legendheightscale * 0.065 * (len(l_hist) + 1.5*(len(legendtitle)>0))
    legendWidth = legendwidthscale * 0.4
    
    padTop = 1 - canvas.GetTopMargin() - 1*ROOT.gStyle.GetTickLength("y")
    padRight = 1 - canvas.GetRightMargin() - 0.6*ROOT.gStyle.GetTickLength("x")
    padBottom = canvas.GetBottomMargin() + 0.6*ROOT.gStyle.GetTickLength("y")
    padLeft = canvas.GetLeftMargin() + 0.6*ROOT.gStyle.GetTickLength("x")
    
    if(legendpos == "UR") :
        
        legend = ROOT.TLegend(padRight-legendWidth, padTop-legendHeight, padRight, padTop)
    
    elif(legendpos == "LR") :
        
        legend = ROOT.TLegend(padRight-legendWidth, padBottom, padRight, padBottom+legendHeight)
    
    elif(legendpos == "LL") :
        
        legend = ROOT.TLegend(padLeft, padBottom, padLeft+legendWidth, padBottom+legendHeight)
    
    elif(legendpos == "UL") :
        
        legend = ROOT.TLegend(padLeft, padTop-legendHeight, padLeft+legendWidth, padTop)
    
    else :
        
        print("Wrong legend position option:", legendpos)
        exit(1)
    
    
    legend.SetNColumns(legendncol)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetHeader(legendtitle)
    legend.SetTextSize(legendtextsize)
    
    stack = ROOT.THStack()
    stack_integral = 0.0
    
    accume_hist = None

    for hist in l_hist :
        hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if include_overflow:
            hist.SetBinContent(1, hist.GetBinContent(1) +  hist.GetBinContent(0)) # lower bin
            hist.SetBinContent(hist.GetNbinsX() , hist.GetBinContent(hist.GetNbinsX()) +  hist.GetBinContent(hist.GetNbinsX() + 1)) # upper bin
            single_bkgr_int = hist.Integral()+hist.GetBinContent(hist.GetNbinsX()+1)+hist.GetBinContent(0)
            stack_integral += single_bkgr_int
        else:
            single_bkgr_int = hist.Integral()
            stack_integral += single_bkgr_int
                    
        if stackdrawopt == 'nostack' and normilize:
            if single_bkgr_int != 0:
                hist.Scale(1.0/single_bkgr_int)

    for hist in l_hist :
        if stackdrawopt != 'nostack' and normilize:
            hist.Scale(1.0/stack_integral)
        stack.Add(hist, "hist")
        legend.AddEntry(hist, hist.GetTitle(), "LPFE")
        if accume_hist==None:
            hist.SetDirectory(0)
            accume_hist = hist.Clone()
            accume_hist.SetDirectory(0)
        else:
            hist.SetDirectory(0)
            accume_hist.Add(hist)

    stack.Draw(stackdrawopt)
    
    accume_hist.SetFillStyle(3004)
    accume_hist.SetFillColor(1)
    # accume_hist.SetLineColor(15)
    # accume_hist.SetLineWidth(0)
    accume_hist.SetMarkerStyle(21)
    accume_hist.SetMarkerSize(0)
 
    if draw_errors:
        accume_hist.Draw("e2same")
    
    stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    stack.SetMinimum(yrange[0])
    stack.SetMaximum(yrange[1])

    for hist in l_hist_overlay :
        if include_overflow:
            hist.SetBinContent(1, hist.GetBinContent(1) +  hist.GetBinContent(0))
            hist.SetBinContent(hist.GetNbinsX() , hist.GetBinContent(hist.GetNbinsX() ) +  hist.GetBinContent(hist.GetNbinsX() + 1))
            
        if stackdrawopt == 'nostack':
            if hist.Integral()!=0: hist.Scale(1.0/(hist.Integral()+hist.GetBinContent(hist.GetNbinsX()+1)+hist.GetBinContent(0)))
        elif normilize_overlay and include_overflow:
            if hist.Integral()!=0: hist.Scale(1.0/(hist.Integral()+hist.GetBinContent(hist.GetNbinsX()+1)+hist.GetBinContent(0)))
        elif normilize_overlay:
            if hist.Integral()!=0: hist.Scale(1.0/(hist.Integral()))
        
        if ratio_mode=="DATA":
            hist.SetOption("HISTP")
        else:
            hist.SetOption("HIST")
        hist.Draw(f"same {hist.GetOption()}")
        legend.AddEntry(hist, hist.GetTitle(), "LPFE")
    
    legend.Draw()
    
    if (ndivisionsx is not None) :
        
        stack.GetXaxis().SetNdivisions(ndivisionsx[0], ndivisionsx[1], ndivisionsx[2], False)
    
    if (ndivisionsy is not None) :
        
        stack.GetYaxis().SetNdivisions(ndivisionsy[0], ndivisionsy[1], ndivisionsy[2], False)
    
    stack.GetXaxis().SetTitle(xtitle)
    #stack.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize("X") * xTitleSizeScale)
    #stack.GetXaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset("X") * 1.1)
    
    stack.GetYaxis().SetTitle(ytitle)
    #stack.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize("Y") * yTitleSizeScale)
    stack.GetYaxis().SetTitleOffset(1)
    
    #stack.SetTitle(title)

    stack.GetXaxis().CenterTitle(centertitlex)
    stack.GetYaxis().CenterTitle(centertitley)
    
    stack.GetXaxis().CenterLabels(centerlabelx)
    stack.GetYaxis().CenterLabels(centerlabely)
    
    canvas.cd(1).SetLogx(logx)
    canvas.cd(1).SetLogy(logy)
    
    canvas.cd(1).SetGridx(gridx)
    canvas.cd(1).SetGridy(gridy)
    
    CMS_lumi.lumiTextSize = 0.9
    CMS_lumi.cmsTextSize = 0.9
    CMS_lumi.relPosX = 0.045
    CMS_lumi.CMS_lumi(pad = canvas.cd(1), iPeriod = 0, iPosX = 0, CMSextraText = CMSextraText, lumiText = lumiText)
    
    accume_hist.SetDirectory(0)
    if signal_to_background_ratio:
        
        canvas.cd(2)
        stack_ratio = ROOT.THStack()
        #h1_xRange_ratio = h1_xRange.Clone()
        #stack_ratio.Add(h1_xRange_ratio)
        for hist in l_hist_overlay :
            hist.SetDirectory(0)
            h1_ratio = hist.Clone()
            h1_ratio.SetDirectory(0)
            h1_ratio.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            if ratio_mode=="B":
                h1_ratio = h1_ratio.Divide(accume_hist)
                stack_ratio.Add(h1_ratio, "HIST")
            elif ratio_mode=="SB":
                SandB = accume_hist.Clone()
                SandB.SetDirectory(0)
                SandB.Add(hist)
                # h1_ratio.Divide(SandB)
                for bin_i in range(hist.GetNcells()):
                    if SandB.GetBinContent(bin_i) != 0:
                        h1_ratio.SetBinContent(bin_i, h1_ratio.GetBinContent(bin_i) / ROOT.TMath.Sqrt(SandB.GetBinContent(bin_i)))
                stack_ratio.Add(h1_ratio, "HIST")
            elif ratio_mode=="DATA":
                SandB = accume_hist.Clone()
                SandB.SetDirectory(0)
                h1_ratio.Divide(SandB)
                stack_ratio.Add(h1_ratio, "HISTPE1")
    
        stack_ratio.Draw("nostack")
        stack_ratio.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        stack_ratio.SetMinimum(yrange_ratio[0])
        stack_ratio.SetMaximum(yrange_ratio[1])
        
        if (ndivisionsx is not None) :
            
            stack_ratio.GetXaxis().SetNdivisions(ndivisionsx[0], ndivisionsx[1], ndivisionsx[2], False)
        
        if (ndivisionsy_ratio is not None) :
            
            stack_ratio.GetYaxis().SetNdivisions(ndivisionsy_ratio[0], ndivisionsy_ratio[1], ndivisionsy_ratio[2], False)
        
        stack_ratio.GetXaxis().CenterLabels(centerlabelx)
        stack_ratio.GetYaxis().CenterLabels(centerlabely)
        
        stack_ratio.GetXaxis().SetLabelSize(0.1)
        stack_ratio.GetYaxis().SetLabelSize(0.1)
        
        stack_ratio.GetXaxis().SetTitle(xtitle_ratio)
        stack_ratio.GetXaxis().SetTitleSize(0.13)
        stack_ratio.GetXaxis().SetTitleOffset(0.91)
        stack_ratio.GetXaxis().CenterTitle(centertitlex)
        
        stack_ratio.GetYaxis().SetTitle(ytitle_ratio)
        stack_ratio.GetYaxis().SetTitleSize(0.126)
        stack_ratio.GetYaxis().SetTitleOffset(0.45)
        stack_ratio.GetYaxis().CenterTitle(centertitley)
        
        canvas.cd(2).SetGridx(gridx)
        canvas.cd(2).SetGridy(gridy)
        
        canvas.cd(2).SetLogx(logx_ratio)
        canvas.cd(2).SetLogy(logy_ratio)
    
    
    if ("/" in outfile) :
        
        outdir = outfile
        outdir = outdir[0: outdir.rfind("/")]
        
        os.system("mkdir -p %s" %(outdir))
    
    # ROOT.gStyle.SetImageScaling(2.)
    canvas.SaveAs(outfile)
    
    return 0


def root_plot2D(
    l_hist,
    outfile,
    xrange,
    yrange,
    l_hist_overlay = [],
    logx = False, logy = False, logz = False,
    title = "",
    xtitle = "", ytitle = "",
    xtitle_ratio = "", ytitle_ratio = "",
    centertitlex = True, centertitley = True,
    centerlabelx = False, centerlabely = False,
    gridx = False, gridy = False,
    ndivisionsx = None, ndivisionsy = None,
    ndivisionsy_ratio = (5, 5, 0), 
    stackdrawopt = "nostack",
    ratio_mode="SB",
    normilize = False,
    normilize_overlay = True,
    legendpos = "UR",
    legendncol = 1,
    legendtextsize = 0.045,
    legendtitle = "",
    legendheightscale = 1.0, legendwidthscale = 1.0,
    ratio_num_den_pairs = [],
    fill_empty_bins_ratio = False,
    signal_to_background_ratio = False,
    CMSextraText = "Private work (CMS simulation)",
    lumiText = "(13 TeV)",
    text_colz = False
) :

    
    ROOT.gROOT.LoadMacro(os.path.split(os.path.realpath(__file__))[0]+"/tdrstyle.C")
    ROOT.gROOT.ProcessLine("setTDRStyle()")
    
    ROOT.gROOT.SetStyle("tdrStyle")
    ROOT.gROOT.ForceStyle(True)
    
    canvas = ROOT.TCanvas("canvas", "canvas", 1300, 1000)
    canvas.UseCurrentStyle()
    
    ROOT.gStyle.SetPaintTextFormat("1.3e")
     
    canvas.Divide(2, 2)
    
    canvas.cd(1).SetPad(0, 0.5, 0.5, 1.0)
    canvas.cd(2).SetPad(0.5, 0.5, 1.0, 1.0)
    canvas.cd(3).SetPad(0, 0, 0.5, 0.5)
    canvas.cd(4).SetPad(0.5, 0, 1.0, 0.5)

    
    canvas.cd(1)
       
    accume_hist = None
    for hist in l_hist :
        if accume_hist==None:
            hist.SetDirectory(0)
            accume_hist = hist.Clone()
            accume_hist.SetDirectory(0)
        else:
            accume_hist.Add(hist)
    
    accume_hist.GetZaxis().SetLabelFont(10)
    if text_colz:
        accume_hist.Draw("colz text")
    else:
        accume_hist.Draw("colz")

    accume_hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    accume_hist.GetYaxis().SetRangeUser(yrange[0], yrange[1])
    
    accume_hist.GetXaxis().SetTitle(xtitle)
    accume_hist.GetYaxis().SetTitle(ytitle)
    accume_hist.GetYaxis().SetTitleOffset(1)
    accume_hist.GetXaxis().CenterTitle(centertitlex)
    accume_hist.GetYaxis().CenterTitle(centertitley)
    accume_hist.GetXaxis().CenterLabels(centerlabelx)
    accume_hist.GetYaxis().CenterLabels(centerlabely)
    
    canvas.cd(1).SetLogx(logx)
    canvas.cd(1).SetLogy(logy)
    canvas.cd(1).SetLogz(logz)
    
    canvas.cd(1).SetGridx(gridx)
    canvas.cd(1).SetGridy(gridy)
    
    CMS_lumi.lumiTextSize = 0.9
    CMS_lumi.cmsTextSize = 0.9
    CMS_lumi.relPosX = 0.045
    CMS_lumi.CMS_lumi(pad = canvas.cd(1), iPeriod = 0, iPosX = 0, CMSextraText = CMSextraText + "  (Bkgr.)", lumiText = lumiText)
    
    
    canvas.cd(4)
    
    accume_hist_sig = None
    for hist in l_hist_overlay :
        
        hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if accume_hist_sig==None:
            hist.SetDirectory(0)
            accume_hist_sig = hist.Clone()
            accume_hist_sig.SetDirectory(0)
        else:
            hist.SetDirectory(0)
            accume_hist_sig.Add(hist)
  
    accume_hist_sig.GetZaxis().SetLabelFont(10)
    if text_colz:
        accume_hist_sig.Draw("colz text")
    else:
        accume_hist_sig.Draw("colz")
    
    accume_hist_sig.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    accume_hist_sig.GetYaxis().SetRangeUser(yrange[0], yrange[1])

    accume_hist_sig.GetXaxis().SetTitle(xtitle)
    accume_hist_sig.GetYaxis().SetTitle(ytitle)
    accume_hist_sig.GetYaxis().SetTitleOffset(1)
    accume_hist_sig.GetXaxis().CenterTitle(centertitlex)
    accume_hist_sig.GetYaxis().CenterTitle(centertitley)
    accume_hist_sig.GetXaxis().CenterLabels(centerlabelx)
    accume_hist_sig.GetYaxis().CenterLabels(centerlabely)
    
    canvas.cd(4).SetLogx(logx)
    canvas.cd(4).SetLogy(logy)
    canvas.cd(4).SetLogz(logz)
    
    canvas.cd(4).SetGridx(gridx)
    canvas.cd(4).SetGridy(gridy)
    
    CMS_lumi.lumiTextSize = 0.9
    CMS_lumi.cmsTextSize = 0.9
    CMS_lumi.relPosX = 0.045
    CMS_lumi.CMS_lumi(pad = canvas.cd(4), iPeriod = 0, iPosX = 0, CMSextraText = CMSextraText + "  (Signal)", lumiText = lumiText)
    
    accume_hist.SetDirectory(0)
    accume_hist_sig.SetDirectory(0)
    if signal_to_background_ratio:
        
        canvas.cd(2)
        
        accume_hist_ratio = accume_hist_sig.Clone()
        accume_hist_ratio.SetDirectory(0)

        if ratio_mode=="B":
            accume_hist_ratio.Divide(accume_hist)
        elif ratio_mode=="SB":
            SandB = accume_hist.Clone()
            SandB.Add(accume_hist_ratio)
            accume_hist_ratio.Divide(SandB) 
         
        _min, _max = get_min_max(accume_hist_ratio)
        nLevels = 10000
        levels_int = numpy.linspace(0, nLevels, num=nLevels)
        levels = _min + (_max - _min) / (nLevels - 1) * levels_int
        accume_hist_ratio.SetContour(nLevels, numpy.array(levels))
        accume_hist_ratio.GetZaxis().SetRangeUser(_min, _max)
        accume_hist_ratio.GetZaxis().SetLabelFont(10)
        
        if text_colz:
            accume_hist_ratio.Draw("colz text")
        else:
            accume_hist_ratio.Draw("colz")
        
        accume_hist_ratio.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        accume_hist_ratio.GetYaxis().SetRangeUser(yrange[0], yrange[1])
        accume_hist_ratio.GetXaxis().SetTitle(xtitle)
        accume_hist_ratio.GetYaxis().SetTitle(ytitle)
        accume_hist_ratio.GetYaxis().SetTitleOffset(1)
        accume_hist_ratio.GetXaxis().CenterTitle(centertitlex)
        accume_hist_ratio.GetYaxis().CenterTitle(centertitley)
        accume_hist_ratio.GetXaxis().CenterLabels(centerlabelx)
        accume_hist_ratio.GetYaxis().CenterLabels(centerlabely)
    
        canvas.cd(2).SetGridx(gridx)
        canvas.cd(2).SetGridy(gridy)
        canvas.cd(2).SetLogx(logx)
        canvas.cd(2).SetLogy(logy)
        canvas.cd(2).SetLogz(logz)
        CMS_lumi.lumiTextSize = 0.9
        CMS_lumi.cmsTextSize = 0.9
        CMS_lumi.relPosX = 0.045
        CMS_lumi.CMS_lumi(pad = canvas.cd(2), iPeriod = 0, iPosX = 0, CMSextraText = CMSextraText + "  (S/(S+B))", lumiText = lumiText)
    
    
    if ("/" in outfile) :
        
        outdir = outfile
        outdir = outdir[0: outdir.rfind("/")]
        
        os.system("mkdir -p %s" %(outdir))
    
    canvas.SaveAs(outfile)
    
    return 0
