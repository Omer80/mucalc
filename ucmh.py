"""
Omer Tzuk, February 2014
Version 0.1

"""
import numpy as np
from PPPC4DMID_Reader import *
from mssm_data_Reader import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.table import Table
import argparse
from progressbar import Bar, ETA, Percentage, ProgressBar
import mucalc
import h5py
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def main():
	pass
	
	
# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
