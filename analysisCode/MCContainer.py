import ROOT
import sys
import math

from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D

sys.path.append('../fitting')
from hist import hist


from optparse import OptionParser
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--doMCLooping', action='store_true', dest='doMCLooping', default=False, help='go!')
parser.add_option('--doRhalphabet', action='store_true', dest='doRhalphabet', default=False, help='go!')
parser.add_option('--doData', action='store_true', dest='doData', default=False, help='go!')
parser.add_option('--doPlots', action='store_true', dest='doPlots', default=False, help='go!')
parser.add_option('--doCards', action='store_true', dest='doCards', default=False, help='go!')

parser.add_option("--lumi", dest="lumi", default = 0.44,help="mass of LSP", metavar="MLSP")
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 6.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.38,help="mass of LSP", metavar="MLSP")
parser.add_option('--qcdClosure', action='store_true', dest='qcdClosure', default=False, help='go!')

(options, args) = parser.parse_args()



class MCContainer:

	def __init__( self , filename, lumi, name, tag, scaleFactor = 1,isData=False):

		self._fn = filename;
		self._tf = ROOT.TFile( self._fn );
		self._tt = self._tf.Get("otree");
		self._scaleFactor = scaleFactor;
		self._lumi = lumi;
		self._isData = isData;
		self._name = name;
		self._tag = tag;


		# Define histogram
		self.h_jetpt = ROOT.TH1F("h_jetpt"+name,"; jet pT (GeV);", 50, 200, 2000);
		self.h_jeteta = ROOT.TH1F("h_jeteta"+name,"; jet #eta;", 20, -3, 3);
		self.h_jetmsd = ROOT.TH1F("h_jetmsd"+name,"; soft drop mass (GeV);", 60, 0, 300);
		self.h_jett21 = ROOT.TH1F("h_jett21"+name,"; #tau_{21};", 20, 0, 1);
		self.h_jett21DDT = ROOT.TH1F("h_jett21DDT"+name,"; #tau_{21}^{DDT};", 20, 0, 1.2);
		self.h_rhoDDT = ROOT.TH1F("h_rhoDDT"+name,"; #rho^{DDT};", 20,float(options.rholo),float(options.rhohi));

		self.h_jetmsd_passcut = ROOT.TH1F("h_jetmsd_passcut"+name,"; soft drop mass (GeV);", 60, 30, 330);
		self.h_rhoDDT_passcut = ROOT.TH1F("h_rhoDDT_passcut"+name,"; #rho^{DDT};", 20,float(options.rholo),float(options.rhohi));

		# Loop
		self.loop();

	def loop( self ):

		DDTCUT = float(options.DDTcut);

		# looping
		nent = self._tt.GetEntries();
		for i in range(self._tt.GetEntries()):

			# preamble
			self._tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			
			if i % self._scaleFactor != 0: continue;		

			jpt = getattr(self._tt,"bst8_PUPPIjet0_pt");
			jeta = getattr(self._tt,"bst8_PUPPIjet0_eta");
			jmsd = getattr(self._tt,"bst8_PUPPIjet0_msd");		
			if jmsd == 0.: jmsd = 0.01;
			weight = self._scaleFactor*self._lumi*getattr(self._tt,"scale1fb")*getattr(self._tt,"kfactor");
			if self._isData: weight = 1;

			if jpt < 500: continue;

			jt21 = getattr(self._tt,"bst8_PUPPIjet0_tau21");
			rhP = math.log(jmsd*jmsd/jpt);                  
			jt21P = jt21 + 0.063*rhP;

			self.h_jetpt.Fill( jpt, weight );
			self.h_jeteta.Fill( jeta, weight );
			self.h_jett21.Fill( jt21, weight );
			self.h_jett21DDT.Fill( jt21P, weight );
			self.h_jetmsd.Fill( jmsd, weight );

			if rhP > 0 and jt21P < DDTCUT: 
				self.h_jetmsd_passcut.Fill( jmsd, weight );
				self.h_rhoDDT_passcut.Fill( rhP, weight );				

		print "\n";

		print "MC yields = ", self.h_jetmsd.Integral(), self.h_jetmsd_passcut.Integral();

	def morphSignal(self,newname,mass,mass_shift,mass_shift_unc,mass_res,mass_res_unc):

		setattr(self,newname, getattr(self,"h_jetmsd_passcut").Clone());
		hist_container = hist( [mass],[getattr(self,newname)] );

		# get new central value
		shift_val = mass - mass*mass_shift;
		tmp_shifted_h = hist_container.shift( getattr(self,newname), shift_val);
		# get new central value and new smeared value
		smear_val = mass_res - 1;
		tmp_smeared_h = hist_container.smear( tmp_shifted_h[0], smear_val)
		if smear_val <= 0: setattr(self,newname+"_central",tmp_smeared_h[1])
		else: setattr(self,newname+"_central",tmp_smeared_h[0])
		
		# get shift up/down
		shift_unc = mass*mass_shift*mass_shift_unc;
		hsys_shift = hist_container.shift( getattr(self,newname+"_central"), shift_unc);
		# get res up/down
		hsys_smear = hist_container.smear( getattr(self,newname+"_central"), mass_res_unc);	


		print shift_val, smear_val, shift_unc, mass_res_unc
		setattr(self,newname+"_shiftUp",hsys_shift[0])
		setattr(self,newname+"_shiftDn",hsys_shift[1])
		setattr(self,newname+"_smearUp",hsys_smear[0])
		setattr(self,newname+"_smearDn",hsys_smear[1])




