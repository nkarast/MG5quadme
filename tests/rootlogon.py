import ROOT

style = ROOT.TStyle('Default','Default style')

style.SetPalette(1)

# use plain black on white colors
iColor = iMode = 0
style.SetFrameBorderMode(iMode)
style.SetCanvasBorderMode(iMode)
style.SetPadBorderMode(iMode)

style.SetLineWidth(2)

style.SetPadColor(iColor)
style.SetCanvasColor(iColor)
style.SetStatColor(iColor)

# set the paper & margin sizes
style.SetPaperSize(20,26)
style.SetPadTopMargin(0.075)
style.SetPadRightMargin(0.15)
style.SetPadBottomMargin(0.15)
style.SetPadLeftMargin(0.15)
style.SetCanvasDefH(600)
style.SetCanvasDefW(700)

# use large fonts
iFont = 42
TextSize = 0.05
style.SetTextFont(iFont)
style.SetTextSize(TextSize)

iFont = 42
TextSize = 0.045
style.SetLabelFont(iFont,'x')
style.SetTitleFont(iFont,'x')
style.SetLabelFont(iFont,'y')
style.SetTitleFont(iFont,'y')
style.SetLabelFont(iFont,'z')
style.SetTitleFont(iFont,'z')
style.SetLabelSize(TextSize,'x')
style.SetTitleSize(TextSize,'x')
style.SetLabelSize(TextSize,'y')
style.SetTitleSize(TextSize,'y')
style.SetLabelSize(TextSize,'z')
style.SetTitleSize(TextSize,'z')

# use bold lines and markers
style.SetMarkerStyle(20)
style.SetMarkerSize(1.2)
style.SetHistLineWidth(2)
style.SetHatchesLineWidth(1)
style.SetLineStyleString(2,'[12 12]') # postscript dashes

# do not display any of the standard histogram decorations
style.SetOptTitle(0)
style.SetOptStat(0)
style.SetOptFit(1111)

# put tick marks on top and RHS of plots
style.SetPadTickX(1)
style.SetPadTickY(1)

# stats box configuration
style.SetStatW(0.20)
style.SetStatX(0.975)
style.SetStatY(0.975)
style.SetStatFormat('4.2g')
style.SetStatBorderSize(1)
style.SetStatFont(iFont)

style.SetEndErrorSize(0)

style.SetTickLength(0.02, 'XYZ')
style.SetNdivisions(406, 'XYZ')

ROOT.gROOT.SetStyle('Default')
ROOT.gROOT.ForceStyle()
