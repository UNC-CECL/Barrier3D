# Script for dune height histogram of batch simulations

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 1
# Updated: 13 September 2020


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, signal


#%% Initialize

Plot = True
Save = False


paths = ['C:/Barrier3D/Output/BatchSims_Group30/BatchSims_2020_1003_2346',
          'C:/Barrier3D/Output/BatchSims_Group30/BatchSims_2020_1003_2351',
          'C:/Barrier3D/Output/BatchSims_Group30/BatchSims_2020_1009_0431',
          'C:/Barrier3D/Output/BatchSims_Group30/BatchSims_2020_1011_0945']



PSNums = [25,
          25,
          25,
          25]


SimNum = np.sum(PSNums) * 25
Stats = np.zeros([SimNum, 4]) 


# Initialize
DuneHeight_5 = []
SCR_5 = []
FastHd_5 = []
SlowHd_5 = []

DuneHeight_13 = []
SCR_13 = []
FastHd_13 = []
SlowHd_13 = []

DuneHeight_21 = []
SCR_21 = []
FastHd_21 = []
SlowHd_21 = []


#%% LOAD DATA 

for q in range(len(paths)):
    
    print(paths[q])
    
    SimN = PSNums[q] * 25 + 1
    print(SimN)
    
    for Sim in range(1,SimN):
           
        # Create File Path
        filename = paths[q] + '/SimData_' + str(Sim) + '.npz'
    
        # Load File
        SimData = np.load(filename, allow_pickle = True)
        
        # Load  Data
        x_s_TS = SimData['x_s_TS']
        SimParams = SimData['SimParams']
        RSLR = SimParams[1]
        TMAX = len(x_s_TS)
        Hd_AverageTS = SimData['Hd_AverageTS']
        BermEl = SimParams[3] 
        DShoreface = SimParams[6] 
        LShoreface = SimParams[7] 
        
        
        
        #%%
        # Calculate shoreline change due exclusively to RSLR
        R_SLR = RSLR[0] * (LShoreface/DShoreface) * 10
        
        # Shoreline Change & Change Rate Over Time    
        scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
        
        # Average Dune Height
        aHd = [a * 10 for a in Hd_AverageTS]
        
        # Filter
        win = 31 #or 35?
        poly = 3
        F_scts = signal.savgol_filter(scts, win, poly) # Savitsky-Golay Filter
        der1 = (signal.savgol_filter(scts, win, poly, deriv=1))
        
        # Optional: Subtract Effect of RSLR on shoreline change rate
        # der1 = der1 - R_SLR
        
        HitDown = [] # Slow-downs
        HitUp = [] # Speed-ups
        
        SlowHd = []
        FastHd = []
        
        window1 = 3 # Maximum allowed length for gaps in slow periods
        window2 = 30 # Minimum length required for slow periods, including gaps
        buffer = 3
        thresh1 = 0.5 # Max slope for slow periods
        thresh2 = 1
        
        # Find slow periods
        der_under = np.where(der1 < thresh1)[0] 
         
        if len(der_under) > 0:
                                            
            gaps = np.diff(der_under) > window1
            peak_start = np.insert(der_under[1:][gaps], 0, der_under[0])
            peak_stop = np.append(der_under[:-1][gaps], der_under[-1])
        
            # for n in range(len(peak_stop)):
            for n in range(1,len(peak_stop)):
                if peak_stop[n] - peak_start[n] > window2:
                    if len(HitDown) == 0:
                        if peak_start[n] > buffer:
                            HitDown.append(peak_start[n])
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp.append(peak_stop[n])
                    else:
                        a = HitUp[-1]
                        b = peak_start[n]   
                        gap_length = b - a
                        gap_slope = np.average(der1[a:b])
                        # gap_slope = np.max(der1[a:b])
                        if gap_length < window2 and gap_slope < thresh2: 
                            if peak_stop[n] < len(scts) - buffer:
                                HitUp[-1] = peak_stop[n]
                            else:
                                del HitUp[-1]
                        else:
                            if peak_start[n] > buffer:
                                HitDown.append(peak_start[n])
                            if peak_stop[n] < len(scts) - buffer:
                                HitUp.append(peak_stop[n])    
        
        # # Find average dune heights for slow and fast periods
        # if len(HitDown) > 0 and len(HitUp) > 0:
        #     if min(HitDown) < min(HitUp):
        #         if len(HitDown) == len(HitUp):
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
                    
        #             HitDown_end = HitDown[1:] + [TMAX]
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown_end[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 FastHd.append(AvgHd)
        #         else:
        #             HitUp_end = HitUp[:] + [TMAX]
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp_end[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
                    
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown[n+1]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 FastHd.append(AvgHd)
        #         FastHd = [np.mean(aHd[0:HitDown[0]])] + FastHd
        #     else:
        #         if len(HitUp) == len(HitDown):
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 FastHd.append(AvgHd)
                    
        #             HitUp_end = HitUp[1:] + [TMAX]
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp_end[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
        #         else:
        #             HitDown_end = HitDown[:] + [TMAX]
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown_end[n]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 FastHd.append(AvgHd)
                    
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp[n+1]
        #                 AvgHd = np.mean(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
        #         SlowHd = [np.mean(aHd[0:HitUp[0]])] + SlowHd
        
        
        
        # # Find median dune heights for slow and fast periods
        # if len(HitDown) > 0 and len(HitUp) > 0:
        #     if min(HitDown) < min(HitUp):
        #         if len(HitDown) == len(HitUp):
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
                    
        #             HitDown_end = HitDown[1:] + [TMAX]
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown_end[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 FastHd.append(AvgHd)
        #         else:
        #             HitUp_end = HitUp[:] + [TMAX]
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp_end[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
                    
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown[n+1]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 FastHd.append(AvgHd)
        #         FastHd = [np.mean(aHd[0:HitDown[0]])] + FastHd
        #     else:
        #         if len(HitUp) == len(HitDown):
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 FastHd.append(AvgHd)
                    
        #             HitUp_end = HitUp[1:] + [TMAX]
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp_end[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
        #         else:
        #             HitDown_end = HitDown[:] + [TMAX]
        #             for n in range(len(HitUp)):
        #                 start = HitUp[n]
        #                 stop = HitDown_end[n]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 FastHd.append(AvgHd)
                    
        #             for n in range(len(HitDown)):
        #                 start = HitDown[n]
        #                 stop = HitUp[n+1]
        #                 AvgHd = np.median(aHd[start:stop])
        #                 SlowHd.append(AvgHd)
        #         SlowHd = [np.mean(aHd[0:HitUp[0]])] + SlowHd
                
                
                
        # Find all dune heights for slow and fast periods
        if len(HitDown) > 0 and len(HitUp) > 0:
            if min(HitDown) < min(HitUp):
                if len(HitDown) == len(HitUp):
                    for n in range(len(HitDown)):
                        start = HitDown[n]
                        stop = HitUp[n]
                        Hd = aHd[start:stop]
                        SlowHd += Hd
                    
                    HitDown_end = HitDown[1:] + [TMAX]
                    for n in range(len(HitUp)):
                        start = HitUp[n]
                        stop = HitDown_end[n]
                        Hd = aHd[start:stop]
                        FastHd += Hd
                else:
                    HitUp_end = HitUp[:] + [TMAX]
                    for n in range(len(HitDown)):
                        start = HitDown[n]
                        stop = HitUp_end[n]
                        Hd = aHd[start:stop]
                        SlowHd += Hd
                    
                    for n in range(len(HitUp)):
                        start = HitUp[n]
                        stop = HitDown[n+1]
                        Hd = aHd[start:stop]
                        FastHd += Hd
                FastHd = [np.mean(aHd[0:HitDown[0]])] + FastHd
            else:
                if len(HitUp) == len(HitDown):
                    for n in range(len(HitUp)):
                        start = HitUp[n]
                        stop = HitDown[n]
                        Hd = aHd[start:stop]
                        FastHd += Hd
                    
                    HitUp_end = HitUp[1:] + [TMAX]
                    for n in range(len(HitDown)):
                        start = HitDown[n]
                        stop = HitUp_end[n]
                        Hd = aHd[start:stop]
                        SlowHd += Hd
                else:
                    HitDown_end = HitDown[:] + [TMAX]
                    for n in range(len(HitUp)):
                        start = HitUp[n]
                        stop = HitDown_end[n]
                        Hd = aHd[start:stop]
                        FastHd += Hd
                    
                    for n in range(len(HitDown)):
                        start = HitDown[n]
                        stop = HitUp[n+1]
                        Hd = aHd[start:stop]
                        SlowHd += Hd
                SlowHd = [np.mean(aHd[0:HitUp[0]])] + SlowHd
                
                
                
        
    #%% STORE HISTOGRAM DATA    
            
        N = Sim % 25 - 1
        if N < 0: N = 24
        
        # Average Dune Height
        aHd = [a * 10 for a in Hd_AverageTS] 
    
        # Shorline Change Rate
        scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
        SCrate = [0]
        for k in range(1,len(scts)):
            SCrate.append(scts[k]- scts[k-1])
        
        if N == 5 -1:
            DuneHeight_5 = DuneHeight_5 + aHd[1:]
            SCR_5 = SCR_5 + SCrate[1:]
            SlowHd_5 = SlowHd_5 + SlowHd
            FastHd_5 = FastHd_5 + FastHd
            
        if N == 13 -1:
            DuneHeight_13 = DuneHeight_13 + aHd[1:]
            SCR_13 = SCR_13 + SCrate[1:]
            SlowHd_13 = SlowHd_13 + SlowHd
            FastHd_13 = FastHd_13 + FastHd
            
        if N == 21 -1:
            DuneHeight_21 = DuneHeight_21 + aHd[1:]
            SCR_21 = SCR_21 + SCrate[1:]
            SlowHd_21 = SlowHd_21 + SlowHd
            FastHd_21 = FastHd_21 + FastHd



#%% CALCULATE RUNNING MEANS

N = 40 # Number of years for running mean

# Shoreline change rm
SCR_5_rm = np.convolve(SCR_5, np.ones((N,))/N, mode='same')
SCR_13_rm = np.convolve(SCR_13, np.ones((N,))/N, mode='same')
SCR_21_rm = np.convolve(SCR_21, np.ones((N,))/N, mode='same')

# Dune rm
DuneHeight_5_rm = np.convolve(DuneHeight_5, np.ones((N,))/N, mode='same')
DuneHeight_13_rm = np.convolve(DuneHeight_13, np.ones((N,))/N, mode='same')
DuneHeight_21_rm = np.convolve(DuneHeight_21, np.ones((N,))/N, mode='same')


    
        
#%% PLOT HISTOGRAMS      
    
if Plot:
    
    ### Dunes
    # Bin Specifications
    BinStart = 0
    BinStop = 1.6
    BinWidth = 0.1
    BinNum = int(((BinStop - BinStart) / BinWidth) + 2)
    Bin = np.linspace(BinStart, BinStop, BinNum)
    
    # 5
    Fig = plt.figure(figsize=(20,5))
    plt.rcParams.update({'font.size':16})
    ax = Fig.add_subplot(1,3,1)
    ax.hist(DuneHeight_5, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')  
    plt.ylabel('PDF', labelpad=15)
    plt.title('5')

    # 13
    ax = Fig.add_subplot(1,3,2)
    ax.hist(DuneHeight_13, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel('Dune Height (m)', labelpad=15)
    plt.title('13')

    # 21
    ax = Fig.add_subplot(1,3,3)
    ax.hist(DuneHeight_21, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')
    plt.title('21')
    
    ### Shoreline Change
    # Bin Specifications
    BinStart = -2.5
    BinStop = 15
    BinWidth = 0.5
    BinNum = int(((BinStop - BinStart) / BinWidth) + 2)
    Bin = np.linspace(BinStart, BinStop, BinNum)
    
    # 5
    Fig = plt.figure(figsize=(20,5))
    plt.rcParams.update({'font.size':16})
    ax = Fig.add_subplot(1,3,1)
    ax.hist(SCR_5, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')
    plt.ylabel('PDF', labelpad=15)
    plt.title('5')

    # 13
    ax = Fig.add_subplot(1,3,2)
    ax.hist(SCR_13, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel('Shoreline Change (m/yr)', labelpad=15)
    plt.title('13')

    # 21
    ax = Fig.add_subplot(1,3,3)
    ax.hist(SCR_21, bins=Bin, density=True)
    ax.xaxis.set_ticks_position('bottom')
    plt.title('21')    
    
    
    ### SCR - Dune Scatter
    
    CI = 0.05 # 95% confidence interval for p_values 
    CM = 'inferno'
    VM = 75 # Warning: Vmax is hard-wired!
    
    Fig = plt.figure(figsize=(20,5))
    
    # 5
    x = SCR_5_rm
    y = DuneHeight_5_rm
    ax = Fig.add_subplot(1,3,1)
    xy = np.vstack([x, y])
    z = stats.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    plt.scatter(x, y, marker='o', s=10, c=z, edgecolors='none', cmap=CM)#, vmax=VM) 
    plt.ylabel('DuneHeight (m)', labelpad=15)
    axes = plt.gca()
    axes.set_ylim([0.05,1.55])
    axes.set_xlim([-0.5,2.5])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    Ktau, p_value_Ktau = stats.kendalltau(x, y)
    if p_value < 0.05 and p_value_Ktau < CI:
        r_sqr = r_value**2
        x1 = np.linspace(np.min(x), np.max(x), 500)
        y1 = slope * x1 + intercept
        plt.plot(x1, y1, 'r')
        plt.title('5: ' + '$R^2$ = ' + str(round(r_sqr,3)) + ', ' + r'$\tau$ = ' + str(round(Ktau,3)), fontsize='small', position = [0.5, 1], color = 'r')
    print(max(z))
    
    # 13
    x = SCR_13_rm
    y = DuneHeight_13_rm
    ax = Fig.add_subplot(1,3,2)
    xy = np.vstack([x, y])
    z = stats.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    plt.scatter(x, y, marker='o', s=10, c=z, edgecolors='none', cmap=CM)#, vmax=VM)
    plt.xlabel('Shoreline Change (m/yr)', labelpad=15)
    axes = plt.gca()
    axes.set_ylim([0.05,1.55])
    axes.set_xlim([-0.5,2.5])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    Ktau, p_value_Ktau = stats.kendalltau(x, y)
    if p_value < 0.05 and p_value_Ktau < CI:
        r_sqr = r_value**2
        x1 = np.linspace(np.min(x), np.max(x), 500)
        y1 = slope * x1 + intercept
        plt.plot(x1, y1, 'r')
        plt.title('13: ' + '$R^2$ = ' + str(round(r_sqr,3)) + ', ' + r'$\tau$ = ' + str(round(Ktau,3)), fontsize='small', position = [0.5, 1], color = 'r')
    print(max(z))
    
    # 21
    x = SCR_21_rm
    y = DuneHeight_21_rm
    ax = Fig.add_subplot(1,3,3)
    xy = np.vstack([x, y])
    z = stats.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    plt.scatter(x, y, marker='o', s=10, c=z, edgecolors='none', cmap=CM)#, vmax=VM)
    axes = plt.gca()
    axes.set_ylim([0.05,1.55])
    axes.set_xlim([-0.5,2.5])
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    Ktau, p_value_Ktau = stats.kendalltau(x, y)
    if p_value < 0.05 and p_value_Ktau < CI:
        r_sqr = r_value**2
        x1 = np.linspace(np.min(x), np.max(x), 500)
        y1 = slope * x1 + intercept
        plt.plot(x1, y1, 'r')
        plt.title('21: ' + '$R^2$ = ' + str(round(r_sqr,3)) + ', ' + r'$\tau$ = ' + str(round(Ktau,3)), fontsize='small', position = [0.5, 1], color = 'r')
    print(max(z))    
   



    ### Fast/Slow - Dune Histogram
    
    Bin = np.linspace(0,1.6,17)
    Fig = plt.figure(figsize=(20,5))
    
    ax = Fig.add_subplot(1,3,1)
    plt.hist(SlowHd_5, color='g', bins = Bin, alpha = 0.5, density=True)
    plt.hist(FastHd_5, color='r', bins = Bin, alpha = 0.5, density=True)
    plt.ylabel('PDF')
    
    ax = Fig.add_subplot(1,3,2)
    plt.hist(SlowHd_13, color='g', bins = Bin, alpha = 0.5, density=True)   
    plt.hist(FastHd_13, color='r', bins = Bin, alpha = 0.5, density=True)
    plt.xlabel('Average Dune Height (m)')
    
    ax = Fig.add_subplot(1,3,3)
    plt.hist(SlowHd_21, color='g', bins = Bin, alpha = 0.5, density=True)
    plt.hist(FastHd_21, color='r', bins = Bin, alpha = 0.5, density=True)
    plt.legend(['Immobile', 'Transgressive'])
    
    