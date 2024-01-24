from positroniumHeader import *


# math functions for fitting, momenta reconstruction, uncertainity

def gaussian(x,A,mu,sig):
    return A * np.exp(-((x-mu)/sig)**2)

def getAngle(v1,v2):
    mag_v1 = np.linalg.norm(v1)
    mag_v2 = np.linalg.norm(v2)
    numer = np.dot(v1,v2)
    denom = mag_v1 * mag_v2
    theta = np.arccos(numer/denom)
    return theta

def additionErrorPropagation(errlist):
    err = np.sqrt(np.sum([x**2 for x in errlist]))
    return err


# analysis functions

# bins should be put in as a list as [number of bins, min_value,max_value]
def getDoubleCoincidenceEnergySpectrum(df,channelID,bins,guess = [1,1,1],fitcut = True,display = False):
    fig,ax = plt.subplots()
       
    if channelID in np.unique(df.ChannelIDL):
        df_by_chan = df[df.ChannelIDL == channelID]
        energy = df_by_chan.ChargeL.to_numpy()
    
    else:
        df_by_chan = df[df.ChannelIDR == channelID]
        energy = df_by_chan.ChargeR.to_numpy()
    
    y,x = np.histogram(energy,bins[0],(bins[1],bins[2]))
    centers = (x[:-1] + x[1:]) / 2
    
    if fitcut == True:
        fitcut = centers[np.where(y == max(y))[0][0]] - 5
        energy_temp = energy[energy >= fitcut]
        y,x = np.histogram(energy_temp,bins[0],(bins[1],bins[2]))
        centers = (x[:-1] + x[1:]) / 2
    
    if guess == [1,1,1]:
        guess = [max(y),centers[np.where(y == max(y))[0][0]],np.std(energy_temp)]
        
    try:
        p,c = curve_fit(gaussian,centers,y,p0=guess)
        xspace = np.linspace(p[1]-2.5*p[2],p[1]+2.5*p[2],500)
        ax.plot(xspace,gaussian(xspace,*p),color = 'red')
        ax.hist(energy,bins = np.linspace(bins[1],bins[2],bins[0]),color = 'C0')
    except:
        p = [1,1,1]
        ('Fit Failed')
        
    if display == False:
        plt.close()
        
    return p

# using this w/ photoPeakFinder saves an insane amount of time!
# im talking a 5 hour run time to a 5 minute run time... cheers, Firas
def getStackedDoubleDataFrame(olddf):
    leftdf = olddf[["TimeL","ChargeL","ChannelIDL"]]
    leftdf.columns = ['Time','Charge','ChannelID']
    rightdf = olddf[["TimeR","ChargeR","ChannelIDR"]]
    rightdf.columns = ['Time','Charge','ChannelID']
    df = pd.concat([leftdf,rightdf])
    df.index = np.arange(0,np.shape(df)[0],1)
    return df

# depreciated for use w/ normal double coincidence frames
'''MUST be used with a stacked dataframe!'''
def photoPeakFinder(stacked_df,channelID,bins,guess = [1,1,1]):

    df_by_chan = stacked_df[stacked_df.ChannelID == channelID]
    energy = df_by_chan.Charge.to_numpy()
    
    y,x = np.histogram(energy,bins[0],(bins[1],bins[2]))
    centers = (x[:-1] + x[1:]) / 2
    
    fitcut = centers[np.where(y == max(y))[0][0]] - 5
    energy_temp = energy[energy >= fitcut]
    y,x = np.histogram(energy_temp,bins[0],(bins[1],bins[2]))
    
    centers = (x[:-1] + x[1:]) / 2
    
    if guess == [1,1,1]:
        guess = [max(y),centers[np.where(y == max(y))[0][0]],np.std(energy_temp)]
        
    try:
        p,c = curve_fit(gaussian,centers,y,p0=guess)
    except:
        p = [1,1,1]
        
    return p
    
def histProjection(data1,data2,bins,labels = ['','']):
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.02, hspace=0.02)
    ax = fig.add_subplot(gs[1, 0])

    # plot theta_nk vs theta_nl
    ax.hist2d(data1,data2,bins = bins)
    ax.set_xlabel(labels[0],fontsize = 19)
    ax.tick_params('x',labelsize = 14)
    ax.set_ylabel(labels[1],fontsize = 19)
    ax.tick_params('y',labelsize = 14)

    # create x & y projection histograms
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx.hist(data1,bins = bins)
    ax_histy.hist(data2, orientation='horizontal',bins = bins)
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    ax_histx.tick_params('y',labelsize = 14)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    ax_histy.tick_params('x',labelsize = 14)
    plt.show()

# same thing as histProjection but scatter plot weighted with gaussian kernal density
# sometimes nicer for visualization but very computational expensive for large datasets, so use w/ caution or small subsets
# seriously **DO NOT USE** for any large datasets...
def scatterProjection(data1,data2,bins,labels = ['','']):
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.02, hspace=0.02)
    ax = fig.add_subplot(gs[1, 0])

    # Calculate the point density
    xy = np.vstack([data1,data2])
    density = gaussian_kde(xy)(xy)
    
    ax.scatter(data1,data2,c=density)
    ax.set_xlabel(labels[0],fontsize = 19)
    ax.tick_params('x',labelsize = 14)
    ax.set_ylabel(labels[1],fontsize = 19)
    ax.tick_params('y',labelsize = 14)

    # create x & y projection histograms
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx.hist(data1,bins = bins)
    ax_histy.hist(data2, orientation='horizontal',bins = bins)
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    ax_histx.tick_params('y',labelsize = 14)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    ax_histy.tick_params('x',labelsize = 14)
    plt.show()
