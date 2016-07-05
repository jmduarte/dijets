import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
from array import array
import os
import glob

import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.08);
ROOT.gStyle.SetPadLeftMargin(0.11);
ROOT.gStyle.SetPadTopMargin(0.10);

ROOT.gStyle.SetPalette(1);

########################################################################################
########################################################################################
def main():

	sXSgb1 = [1.394e+05,8.419e+04,4.481e+04,2.641e+04,1.939e+04,1.462e+04,9976,7870,5707,4254,3233,2320,1131, 620]
	masses = [       50,       60,       75,       90,      100,      110, 125, 135, 150, 165, 180, 200, 250, 300];


	for i,m in enumerate(masses):

		print i,m
		# if i > 2: break;

		nui_cmmd = "python diffNuisances.py datacards/mlfitZprime%i.root -g datacards/output_%03i.root" % (m,m);
		os.system(nui_cmmd);
		os.system("cp nuisances.pdf limplots/nuisances_%i.pdf" % (m));
		postfitPlot_cmmd = "python plotPostfit.py --data datacards/combine_Zprime%i_0.root --input datacards/mlfitZprime%i.root --name limplots/postfitresults_%i" % (m,m,m)
		os.system(postfitPlot_cmmd);


if __name__ == '__main__':

	main();