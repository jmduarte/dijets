import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
from array import array
import os
import glob

def main():
	
	f = ROOT.TFile("mlfit.root")
	t = f.Get("Events");

	

if __name__ == '__main__':

	main();