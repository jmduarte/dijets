import ROOT
import sys
import math

from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 4.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.45,help="mass of LSP", metavar="MLSP")

(options, args) = parser.parse_args()

class ralphabet:

	def __init__( self , filename, lumi, name, scaleFactor = 1, extractTFs = True):

		self._filename = filename;
		self._tf = ROOT.TFile( self._filename );
		self._tt = self._tf.Get("otree");
		self._sf = scaleFactor;

		tfopt = "READ";
		if extractTFs: tfopt = "RECREATE"
		self._storageFile = ROOT.TFile("storedHistos.root",tfopt);

		if extractTFs: 
			self.FillTFMap();
			self._storageFile.cd();
			hys = makeCanvas2D( self.TF_pafa,"map","plots" );
			for h in hys: h.Write();
			self.TF_pafa.Write();

		# else:
		# 	# get stored information

		self._storageFile.Close();

	def FillTFMap( self ):

		DDTCUT = float(options.DDTcut);

		self.TF_fail = ROOT.TH2F("TF_fail",";rho;pT",20,float(options.rholo),float(options.rhohi),10,500,1500)
		self.TF_pafa = ROOT.TH2F("TF_pafa",";rho;pT",20,float(options.rholo),float(options.rhohi),10,500,1500)	
		self.TF_fail.Sumw2();
		self.TF_pafa.Sumw2();

		t = self._tt;

		# looping
		nent = t.GetEntries();
		for i in range(t.GetEntries()):

			# preamble
			t.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			if i % self._sf != 0: continue;

			# cutting
			jpt = getattr(t,"bst8_PUPPIjet0_pt");
			jmsd = getattr(t,"bst8_PUPPIjet0_msd");		
			weight = getattr(t,"scale1fb");

			if jmsd > 0:

				jt21 = getattr(t,"bst8_PUPPIjet0_tau21");
				rhP = math.log(jmsd*jmsd/jpt);			
				jt21P = jt21 + 0.063*rhP;

				if jt21P < DDTCUT: self.TF_pafa.Fill( rhP, jpt, weight );
				else: self.TF_fail.Fill( rhP, jpt, weight );

		print "\n";

		#####!!!! this uncertinaty may be wrong, binomial instead of poisson
		self.TF_pafa.Divide( self.TF_fail );

