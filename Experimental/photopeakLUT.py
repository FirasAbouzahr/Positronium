# Firas Abouzhr
# Hello, Yale, if you are reading this... please accept me into your graduate schools. Thanks (:
'''

The 3 photons from an o-Ps self-annihilation event must have an energy (within some uncertainity) summed energy of 1022 keV.
The data we are using is not energy calibrated, meaning that charge is read out in DAQ units. This means that each SiPM, due
to their inherent gain differences, would all read, say 511 keV, at different DAQ values. Hence when imposing a 1022 keV
energy cut on our triple coincidences, we must consider each SiPM involved in detection and calculate a relative 1022 keV
value that is specfic to the three SiPMs.

Here we create a look-up-table of all the photopeak (so 511 keV) location for each SiPM, so we may reference it for the energy
cuts later down the pipeline.
'''


import sys
sys.path.append('/Users/feef/Documents/GitHub/Positronium')

from positroniumHeader import *
from statHeader import *
save_as_csv = True



# read-in triple coincidence data
df = getTripleCoincidenceDataFrame("/Users/feef/Positronium/TCoinc-Tests/LightTightness7-21-23_quadRegion_Tripcoinc.dat",toGeo = False)
print("Reading in triples...")

# read in double coincidence data so we may get the photopeak locations for all of our channels
file = "/Users/feef/Positronium/TCoinc-Tests/LightTightness7-21-23_coinc.dat"
doubleCoincDf  = pd.read_csv(file,sep='\t',usecols=[2,3,4,7,8,9],chunksize = 10000000)
print("Reading in doubles...")
print("This may take awhile... but it's okay we only need to do this once for an dataset, so be patient and consider the following:")
print("where the hell are all the right-handed neutrinos yet!?")


framelist = []
for frames in doubleCoincDf:
    framelist.append(frames)
doubleCoincDf = pd.concat(framelist)
doubleCoincDf.columns = ['TimeL', 'ChargeL', 'ChannelIDL', 'TimeR', 'ChargeR', 'ChannelIDR']

# this takes awhile *sigh* but is necessary.. Positronium analysis should always be done in recoID
doubleCoincDf["ChannelIDL"] = doubleCoincDf["ChannelIDL"].apply(toRecoChannelID)
doubleCoincDf["ChannelIDR"] = doubleCoincDf["ChannelIDR"].apply(toRecoChannelID)

# get photopeak locations for every channel across all triple LORs
photopeakDict = {"Channel":[],"Photopeak Location":[],'std':[]}
triplecoincidencechannels = np.unique(np.concatenate((np.unique(df.ChannelID1),np.unique(df.ChannelID2),np.unique(df.ChannelID3))))
stackeddf = getStackedDoubleDataFrame(doubleCoincDf)

print("Locating photopeaks and generating the LUT")
for chans in tqdm(triplecoincidencechannels):
    p = photoPeakFinder(stackeddf,chans,[100,0,40])
    photopeakDict["Channel"]+=[chans]
    photopeakDict["Photopeak Location"]+=[p[1]]
    photopeakDict["std"]+=[p[2]]
    
photopeakDf = pd.DataFrame(photopeakDict)
if save_as_csv == True:
    photopeakDf.to_csv("photopeakLUT.csv")
    
