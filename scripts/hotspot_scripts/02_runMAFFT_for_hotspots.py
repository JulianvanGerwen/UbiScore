#This script runs MAFFT on all sequences prior to hotspot detection analysis
# This must be run from my mafft environment

import os
import glob
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

import MAFFT_fns

#Run MAFFT
MAFFT_fns.runMAFFT_and_unravel("../../output/data/hotspots/domains/domainsPreAlign", "../../output/data/hotspots/domains/domainsPostAlign", "../../output/data/hotspots/domains/domainsPostAlignUnwrap")