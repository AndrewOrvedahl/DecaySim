from __future__ import division, print_function
from ROOT import TH1F, TH2F, TTree, TNtuple, TCanvas, TFile
from ROOT import gStyle, gROOT
import sys
import os


class Plotter(object):
    """Histogram plotter object"""


    def __init__(self, title='', filename=None, filepath='plots/'):

        if( filepath.endswith('/') ):
            self.filepath = filepath
        else:
            self.filepath = filepath + '/'
        try:
            os.mkdir(filepath)
        except OSError:
            #The directory already exists
            pass
        if( filename is not None ):
            self.tfile = TFile(self.filepath+filename, 'RECREATE', title)
        else:
            self.tfile = None
        self.plots1D = {}
        self.plots2D = {}
        return


    def SetPalette(self, optstat=1111111):
        """Sets default plot settings."""
        gROOT.Reset()
        gROOT.SetStyle('Plain')
        gStyle.SetOptStat(optstat)
        gStyle.SetPalette(1)
        gStyle.SetLineColor(4)
        return


    def BookNTuple(self, title, vals):
        "Book a TNtuple."
        self.ntuple = TNtuple(title, title, ':'.join(vals))
        return


    def FillNTuple(self, values):
        """Fill the ntuple."""
        exec('self.ntuple.Fill%s' % str(tuple(values)))
        return
    

    def AddPlot1D(self, name, nbins, minimum, maximum, title='', xtitle='',
                  ytitle=''):
        """Book a 1D histogram."""
        hist = TH1F(name, title, nbins, minimum, maximum)
        hist.Sumw2()
        hist.SetXTitle(xtitle)
        hist.SetYTitle(ytitle)
        self.plots1D[name] = hist

    
    def AddPlot2D(self, name, nxbins, xmin, xmax, nybins, ymin, ymax, title='',
                  xtitle='', ytitle=''):
        """Book a 2D histogram."""
        hist = TH2F(name, title, nxbins, xmin, xmax, nybins, ymin, ymax)
        hist.Sumw2()
        hist.SetXTitle(xtitle)
        hist.SetYTitle(ytitle)
        self.plots2D[name] = hist
        return
    

    def FillHist1D(self, name, val):
        """Fill a 1D histogram."""
        self.plots1D[name].Fill(val)
        return


    def FillHist2D(self, name, xval, yval):
        """Fill a 2D histogram."""
        self.plots2D[name].Fill(xval, yval)
        return


    def NormalizePlots(self):
        """Normalize plots to 1."""
        if( sys.version_info[0] == 2 ):
            for name,hist in self.plots1D.iteritems():
                try:
                    hist.Scale(1/hist.Integral())
                except ZeroDivisionError:
                    pass
            for name,hist in self.plots2D.iteritems():
                try:
                    hist.Scale(1/hist.Integral())
                except ZeroDivisionError:
                    pass

        if ( sys.version_info[0] == 3 ):
            for name,hist in self.plots1D.items():
                try:
                    hist.Scale(1/hist.Integral())
                except ZeroDivisionError:
                    pass
            for name,hist in self.plots2D.items():
                try:
                    hist.Scale(1/hist.Integral())
                except ZeroDivisionError:
                    pass

        return
    

    def PlotsToPDF(self, filetype='pdf'):
        if( sys.version_info[0] == 2 ):
            for name,hist in self.plots1D.iteritems():
                c = TCanvas()
                hist.Draw('E')
                c.SetLeftMargin(0.15)
                c.SetBottomMargin(0.15)
                hist.GetYaxis().SetTitleOffset(1.6)
                hist.GetXaxis().SetTitleOffset(1.6)
                exec('c.SaveAs("%s."+filetype, filetype)' %
                     (self.filepath+name))
            for name,hist in self.plots2D.iteritems():
                c = TCanvas()
                hist.Draw('LEGO2 E')
                c.SetLeftMargin(0.15)
                c.SetBottomMargin(0.15)
                hist.GetYaxis().SetTitleOffset(1.6)
                hist.GetXaxis().SetTitleOffset(1.6)
                exec('c.SaveAs("%s."+filetype, filetype)' %
                     (self.filepath+name))

        if( sys.version_info[0] == 3 ):
            for name,hist in self.plots1D.items():
                c = TCanvas()
                hist.Draw('E')
                c.SetLeftMargin(0.15)
                c.SetBottomMargin(0.15)
                hist.GetYaxis().SetTitleOffset(1.6)
                hist.GetXaxis().SetTitleOffset(1.6)
                exec('c.SaveAs("%s."+filetype, filetype)' %
                     (self.filepath+name))
            for name,hist in self.plots2D.items():
                c = TCanvas()
                hist.Draw('LEGO2 E')
                c.SetLeftMargin(0.15)
                c.SetBottomMargin(0.15)
                hist.GetYaxis().SetTitleOffset(1.6)
                hist.GetXaxis().SetTitleOffset(1.6)
                exec('c.SaveAs("%s."+filetype, filetype)' %
                     (self.filepath+name))
        return
    

    def SaveFile(self):
        """Save plots to file."""
        if( self.tfile ):
            self.tfile.Write()
        return

    def CloseFile(self):
        if( self.tfile ):
            self.tfile.Close()
            self.tfile=None
        return
