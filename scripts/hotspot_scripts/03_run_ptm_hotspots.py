# Script to run PTM hotspot analysis
# Must be run from my ptm_hotspots environment

import subprocess

#Define commans as lists
#Run both hotspots and resiude-level results
commandHotspots = [
    "python", "../ptm_hotspots/ptm_hotspots.py",
    "-o", "../../output/data/hotspots/domains/hotspotOutput/hotspotOutput.csv",
    "--dir", "../../output/data/hotspots/domains/domainsPostAlignUnwrap",
    "--ptmfile", "../../output/data/hotspots/ubnSites_forHotspots",
    "--aa_letters", "K"
]
commandResidues = [
    "python", "../ptm_hotspots/ptm_hotspots.py",
    "-o", "../../output/data/hotspots/domains/hotspotOutput/hotspotResidueOutput.csv",
    "--dir", "../../output/data/hotspots/domains/domainsPostAlignUnwrap",
    "--ptmfile", "../../output/data/hotspots/ubnSites_forHotspots",
    "--aa_letters", "K",
    "--printSites"
]


#Run commands
subprocess.run(commandHotspots, check=True)
subprocess.run(commandResidues, check=True)