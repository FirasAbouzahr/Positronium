# Firas Abouzahr

import sys
sys.path.append('/Users/feef/Documents/GitHub/Positronium')
from positroniumHeader import *
from statHeader import *

df = getTripleCoincidenceDataFrame("/Users/feef/Positronium/TCoinc-Tests/LightTightness7-21-23_quadRegion_Tripcoinc.dat",toGeo = False)

try:
    photopeakDf = pd.read_csv("photopeakLUT.csv")
except:
    print("A photopeak LUT does not yet exist... \n Please run photopeakLUT.py")

# grab triple coincidences that meet annhilation energy threshold
error = []
ene = []
sigma = 2 # we filter out triple coincidences with energy that does not meet this threshold

energycutdf = pd.DataFrame()
for trips in tqdm(df.index):
    tempdf = df[df.index == trips]
    ID1,ID2,ID3 = tempdf.ChannelID1.iloc[0],tempdf.ChannelID2.iloc[0],tempdf.ChannelID3.iloc[0]
    annhilationEnergy = 2/3*(photopeakDf[photopeakDf.Channel == ID1]["Photopeak Location"].iloc[0] + photopeakDf[photopeakDf.Channel == ID2]["Photopeak Location"].iloc[0] + photopeakDf[photopeakDf.Channel == ID3]["Photopeak Location"].iloc[0])
    uncertainity = 2/3*additionErrorPropagation([photopeakDf[photopeakDf.Channel == ID1]["std"].iloc[0]*sigma,photopeakDf[photopeakDf.Channel == ID2]["std"].iloc[0]*sigma,photopeakDf[photopeakDf.Channel == ID3]["std"].iloc[0]*sigma])
    
    ene.append(annhilationEnergy)
    error.append(uncertainity)
    
    totalCharge = tempdf.Charge1.iloc[0] + tempdf.Charge2.iloc[0] + tempdf.Charge3.iloc[0]
    
    # only accept triple coincidences whose summed energy meets the 2*photopeak threshold within some uncertainity
    if totalCharge >= annhilationEnergy - uncertainity and totalCharge <= annhilationEnergy + uncertainity:
        energycutdf = pd.concat([energycutdf,tempdf])


error_limits = np.mean(np.array(error)/np.array(ene) * 1022)
print('The required uncertainity is within ± {} keV of 1022 keV'.format(error_limits))

energycutdf.to_csv("triple_coincidences_energycut_uncertainity={}.csv".format(int(error_limits)),index = False)
