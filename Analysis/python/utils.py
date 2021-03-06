import dataclasses
import numpy
import os

import ROOT

import CMS_lumi


def get_canvas(ratio = False) :
    
    ROOT.gROOT.LoadMacro("utils/tdrstyle.C")
    ROOT.gROOT.ProcessLine("setTDRStyle()")
    
    ROOT.gROOT.SetStyle("tdrStyle")
    ROOT.gROOT.ForceStyle(True)
    
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    canvas.UseCurrentStyle()
    
    #canvas.SetLeftMargin(0.16)
    #canvas.SetRightMargin(0.05)
    #canvas.SetTopMargin(0.1)
    #canvas.SetBottomMargin(0.135)
    
    if (ratio) :
        
        canvas.Divide(1, 2);
        
        canvas.cd(1).SetPad(0, 0.32, 1, 1);
        canvas.cd(1).SetTopMargin(0.075)
        canvas.cd(1).SetBottomMargin(0)
        
        canvas.cd(2).SetPad(0, 0.0, 1, 0.3);
        canvas.cd(2).SetTopMargin(0.05)
        canvas.cd(2).SetBottomMargin(0.285)
    
    canvas.cd(1).SetLeftMargin(0.125)
    canvas.cd(1).SetRightMargin(0.05)
    
    if (ratio) :
        
        canvas.cd(2).SetLeftMargin(0.125)
        canvas.cd(2).SetRightMargin(0.05)
    
    return canvas


def root_plot1D(
    l_hist,
    outfile,
    xrange,
    yrange,
    logx = False, logy = False,
    title = "",
    xtitle = "", ytitle = "",
    xtitle_ratio = "", ytitle_ratio = "",
    yrange_ratio = (0, 1),
    centertitlex = True, centertitley = True,
    centerlabelx = False, centerlabely = False,
    gridx = False, gridy = False,
    ndivisionsx = None, ndivisionsy = None,
    ndivisionsy_ratio = (5, 5, 0), 
    stackdrawopt = "nostack",
    legendpos = "UR",
    legendncol = 1,
    legendtextsize = 0.045,
    legendtitle = "",
    legendheightscale = 1.0, legendwidthscale = 1.0,
    ratio_num_den_pairs = [],
    CMSextraText = "Simulation Preliminary",
    lumiText = "(13 TeV)"
) :
    
    canvas = get_canvas(ratio = len(ratio_num_den_pairs))
    
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
    
    for hist in l_hist :
        
        hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        #hist.SetFillStyle(0)
        
        stack.Add(hist, "hist")
        legend.AddEntry(hist, hist.GetTitle(), "LP")
    
    # Add a dummy histogram so that the X-axis range can be beyond the histogram range
    h1_xRange = ROOT.TH1F("h1_xRange", "h1_xRange", 1, xrange[0], xrange[1])
    stack.Add(h1_xRange)
    
    stack.Draw(stackdrawopt)
    legend.Draw()
    
    stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    stack.SetMinimum(yrange[0])
    stack.SetMaximum(yrange[1])
    
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
    
    canvas.SetLogx(logx)
    canvas.SetLogy(logy)
    
    canvas.SetGridx(gridx)
    canvas.SetGridy(gridy)
    
    CMS_lumi.lumiTextSize = 0.9
    CMS_lumi.cmsTextSize = 0.9
    CMS_lumi.relPosX = 0.045
    CMS_lumi.CMS_lumi(pad = canvas.cd(1), iPeriod = 0, iPosX = 0, CMSextraText = CMSextraText, lumiText = lumiText)
    
    
    if (len(ratio_num_den_pairs)) :
        
        canvas.cd(2)
        
        stack_ratio = ROOT.THStack()
        
        #h1_xRange_ratio = h1_xRange.Clone()
        #stack_ratio.Add(h1_xRange_ratio)
        
        for h1_num, h1_den in ratio_num_den_pairs :
            
            h1_ratio = h1_num.Clone()
            h1_ratio.Divide(h1_den)
            
            h1_ratio.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            
            #print(h1_num.Integral())
            #print(h1_den.Integral())
            
            stack_ratio.Add(h1_ratio, "hist")
        
        
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
        
        canvas.SetGridx(gridx)
        canvas.SetGridy(gridy)
    
    
    if ("/" in outfile) :
        
        outdir = outfile
        outdir = outdir[0: outdir.rfind("/")]
        
        os.system("mkdir -p %s" %(outdir))
    
    canvas.SaveAs(outfile)
    
    return 0
