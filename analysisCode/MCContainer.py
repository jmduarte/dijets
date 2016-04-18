import ROOT
import sys

from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D


from optparse import OptionParser
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 4.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.45,help="mass of LSP", metavar="MLSP")

(options, args) = parser.parse_args()



class MCContainer:

	def __init__( self , filename, lumi, name, scaleFactor = 1):

		self._fn = filename;
		self._tf = ROOT.TFile( self._fn );
		self._tt = self._tf.Get("otree");
		self._scaleFactor = scaleFactor;

		# Define histogram
		self.h_jetpt = ROOT.TH1F("h_jetpt"+name,"; jet pT (GeV);", 50, 200, 2000);
		self.h_jeteta = ROOT.TH1F("h_jeteta"+name,"; jet #eta;", 20, -3, 3);
		self.h_jetmsd = ROOT.TH1F("h_jetmsd"+name,"; soft drop mass (GeV);", 30, 0, 300);
		self.h_rhoDDT = ROOT.TH1F("h_rhoDDT"+name,"; #rho^{DDT};", 20, options.rholo, options.rhohi);

		# Loop
		self.loop();

	def loop( self ):

		# looping
		nent = self._tt.GetEntries();
		for i in range(self._tt.GetEntries()):

			self._tt.GetEntry(i);

			# preamble
			self._tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			
			if i % self._scaleFactor != 0: continue;		

			jpt = getattr(self._tt,"bst8_PUPPIjet0_pt");
			jeta = getattr(self._tt,"bst8_PUPPIjet0_eta");
			jmsd = getattr(self._tt,"bst8_PUPPIjet0_msd");		
			weight = getattr(self._tt,"scale1fb");

			self.h_jetpt.Fill( jpt, weight );
			self.h_jeteta.Fill( jeta, weight );
			self.h_jetmsd.Fill( jmsd, weight );

		print "\n";




