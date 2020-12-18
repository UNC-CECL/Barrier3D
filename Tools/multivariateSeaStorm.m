% ------------------------multivariateSeaStorm---------------------------%
% Purpose: This function creates synthetic time series from tidal guage and
% WIS buoy data using the method of Wahl et al., 2016
%
% SEE ALSO: t_tide
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       07/09/20        KAnarde             Original code
%       07/29/20        IReeves             Tidal sampling & simulation TWL calculations
%
%% preamble

% NEED TO MAKE THIS HAVE MORE FUCNTIONATLITY BEFORE PUBLISHING: boolean for
% different copulas, sub-functions, routines to get data from the WIS 
% server etc.

%%% load WIS data 1980-2014
% KA: from Wahl et al., 2016 - need Hs, Tp, and direction theta 
% (they used offshore wave buoy in 28 m water depth, 1980-2013)
% for VCR - Wis ST63183 - the water depth is 22 m

cd /Users/KatherineAnardeWheels/Research/BARis/UNC/VCR/SyntheticStorms
rData = load('ST63183_v03.onlns');
dtH   = datetime(string(rData(:,1)),'InputFormat','yyyyMMddHHmmss');
rHs   = rData(:,10);
rTp   = rData(:,12);  % 11 = TPD, 12 = TP
rWavD = rData(:,16);

%%% load Water Level from CSV (none of this worked, just ended up copying/pasting
% KA: from Wahl et al., 2016 - need one hour time records of water level 
% (they used tidal gauge)
% Trim off first and last readings from tide gauge to match waves
%dtSL = readtable('DateTimes-8631044-Raw.csv');
%dtSL = dtSL(2:end-1, 1)
load('Tide-8631044-Raw.mat', 'rSL', 'dtSL') % dt is in datetime

%% process the tidal data (b/c of data gaps)

rSL_nan = nan(size(dtH));       % for VCR, dtH (2014) and dtSL (2015)
dtSL_nan = nan(size(dtH));
[~, loc] = ismember(dtSL,dtH);  % what is the id for mapping dtSL to dtH
idLoc = loc(loc>0);
rSL_nan(idLoc) = rSL(1:length(idLoc));
dtSL_nan(idLoc) = datenum(dtSL(1:length(idLoc))); % this is for t-tide only

% figure; plot(dtH, rSL_nan, dtSL, rSL)  % for debugging

%% remove non-stationarity

% KA: Wahl et al., 2016 - 30 (or 365) day running medians
%N = 24 * 30;  % For 30 day running mean
N = 24 * 365;  % For 1 yr running mean

%rSL_rm = movmean(rSL_nan, N, 'includenan');
rSL_rm = rSL_nan - movmedian(rSL_nan, N, 'omitnan');
rHs_rm = rHs - movmedian(rHs, N, 'omitnan');
rTp_rm = rTp - movmedian(rTp, N, 'omitnan');

% what is the median of the last 3 years? add this median back to the 
% corrected time series so it is representative of recent climate
idYrs = find(dtH == (dtH(end)-calmonths(12*3)));
rSL_corr = rSL_rm + median(rSL_nan(idYrs:end), 'omitnan');
rHs_corr = rHs_rm + median(rHs(idYrs:end));
rTp_corr = rTp_rm + median(rTp(idYrs:end));

% debugging
min(rTp_corr) % these get negative with 30-day medians, set to zero just in case
min(rHs_corr)
rHs_corr(rHs_corr<0) = 0;
rTp_corr(rTp_corr<0) = 0;

figure
% SL
subplot(4,2,1)
plot(dtH, rSL_nan, dtH, rSL_corr) %, dtH, movmedian(rSL_nan, N, 'omitnan'))
ylabel('Sea level [mNAVD88]')

% Hs
subplot(4,2,2)
plot(dtH, rHs, dtH, rHs_corr) %, dtH, movmedian(rHs, N, 'omitnan'))
ylabel('Hs [m]')
ylim([0 10])

% Tp
% KA: note, I think this Tp data is bad between 2013-2014 
subplot(4,2,3)
plot(dtH, rTp, dtH, rTp_corr) %, dtH, movmedian(rTp, N, 'omitnan'))
ylabel('Tp [s]')
ylim([0 20])

% WavD
subplot(4,2,4)
plot(dtH, rWavD)
ylabel('Wave Direction [degree]')

hold on

%% calculate non-tidal and tidal residuals

% KA: From the corrected SL time series, use t-tide python to perform a 
% year-by-year tidal analysis

% split data separated by NaNs into chunks
rChunk = [rSL_corr, dtSL_nan];  % C1 = corrected SL, C2 = dt
idx = all(isnan(rChunk),2);
idy = 1+cumsum(idx);
idz = 1:size(rChunk,1);
cChunk = accumarray(idy(~idx),idz(~idx),[],@(r){rChunk(r,:)});
cChunkTide = cell(size(cChunk));

% do t-tide for each chunk in one year intervals
for iChunk = 1:length(cChunk)
    
    if  ~isempty(cChunk{iChunk})
        
        % number of years in chunk (sadly we will lose data at the tail; 
        % will want to fix in future)
        tmpSLcorr = cChunk{iChunk}(:,1);
        tmpDT     = cChunk{iChunk}(:,2);
        nStart    = datetime(tmpDT(1),'ConvertFrom','datenum');
        nEnd      = datetime(tmpDT(end),'ConvertFrom','datenum');
        nYrs      = floor(days(nEnd - nStart) / 365); 
        iTideSt   = 1;
        rTideOut  = [];
        
        if nYrs > 0 % I would really like to not throw out data...come back to this
            
            for iYear = 1:nYrs

                iTideEnd = iYear * 24 * 365; % probably a better way to do this with leap years

                % do t-tide predictions
                datetime(tmpDT(iTideSt), 'ConvertFrom', 'datenum') % for debugging
                datetime(tmpDT(iTideEnd), 'ConvertFrom', 'datenum')
                [~, rTideOut(iTideSt:iTideEnd)] = ...
                                     t_tide(tmpSLcorr(iTideSt:iTideEnd), ...
                                     'interval', 1, ...               % hours
                                     'start', tmpDT(iTideSt),...      % datenum
                                     'latitude', 37.5);               % lat

                iTideSt = iTideEnd + 1; 

            end
            
            % save tidal output to new cell array  
            cChunkTide{iChunk} = [tmpSLcorr(1:iTideEnd), ... % subset corr SL
                                  tmpDT(1:iTideEnd), ...     % subset datenum
                                  rTideOut', ...             % tidal prediction
                                  tmpSLcorr(1:iTideEnd)-rTideOut']; % nontidal residual

        end   
        
    end
end

% condense data to remove empty cells
tmpChunk = cell2mat(cChunkTide);
rSLcorr_sub = tmpChunk(:,1);
dtSLcorr_sub = datetime(tmpChunk(:,2), 'ConvertFrom', 'datenum');
rAT = tmpChunk(:,3);
rNTR = tmpChunk(:,4);

% now do as before and map to dtH 
[rNTR_nan, rAT_nan] = deal(nan(size(dtH)));    
[~, loc] = ismember(dtSLcorr_sub,dtH);  % what is the id for mapping dtSLcorr_sub to dtH
idLoc = loc(loc>0);
rNTR_nan(idLoc) = rNTR(1:length(idLoc));
rAT_nan(idLoc) = rAT(1:length(idLoc));

% for debugging
% figure; plot(dtSLcorr_sub, rSLcorr_sub, ... % corr SL
%             dtSLcorr_sub, rAT, ... % pred tide
%             dtSLcorr_sub, rNTR) % residual

% update plots
% rNTR
subplot(4,2,5)
plot(dtH, rNTR_nan)
%plot(dtSLcorr_sub, rNTR)
ylabel('\eta_{NTR} [m]')

% rAT
subplot(4,2,6)
plot(dtH, rAT_nan)
%plot(dtSLcorr_sub, rAT)
ylabel('\eta_{A} [m NAVD88]')

hold on

%% calculate R2% and add to SL to get the TWL (currently only corrected data)

rBeta = 0.04;                      % beach slope, Hog Island    
rL0   = (9.8 * rTp_corr.^2) / (2 * pi); % wavelength       

% KA: is this Stockdon 2006 broken down into components?
rSetup = 0.35 * rBeta * sqrt(rHs_corr .* rL0); 
rSin   = 0.75 * rBeta * sqrt(rHs_corr .* rL0);  % incident band swash
rSig   = 0.06 * sqrt(rHs_corr .* rL0) ;         % infragravity band swash
rSwash = sqrt((rSin.^2) + (rSig.^2));      % total swash
rR2    = 1.1 * (rSetup + (rSwash/2));      % R2%

% KA: not clear from Wahl if the TWL is the corrected SL+R2...do both?
rTWL  = rSL_corr + rR2;       % corrected for nonstationarity
%rTWL     = rSL_nan + rR2;    % observed
rRlow = (rTWL - (rSwash/2));  % this is just for Ian...use observed

% update plots
% R2
subplot(4,2,7)
plot(dtH, rR2)
ylabel('R2/% [m]')

% TWL
subplot(4,2,8)
plot(dtH, rTWL)
ylabel('TWL [m NAVD88]')
hold off

%% Extract sea-storm events from the observational record

%%% First, define TWL threshold exceedence. Then, how much the different 
%%% variables contributed to those events and if there is a dominant driver 
%%% that can be used for the event selection

% for each year, find when the TWL exceeds an erosion threshold
% Wahl used the 5th % of dune toe heights (free parameter for dune erosion)
rBermEl = 2.0;     % m NAVD88

% find the annual average TWL from all threshold exceedances from a 
% given year calculate annual averages of MSL (here 35-day mean), 
% tidal amplitude, residual, and R2% during the TWL exceedances
nYrs = floor(days(dtH(end) - dtH(1))/ 365);
iStart = 1;
[rHs_over_yearly, rR2_over_yearly, rNTR_over_yearly, rTWL_over_yearly] = ...
    deal(NaN(nYrs,1));

% use only corrected data
for iYear = 1 : nYrs
    
    iStop = iYear * 24 * 365;
    
    rHH = rHs_corr(iStart:iStop);      
    rRR = rR2(iStart:iStop); 
    rNN = rNTR_nan(iStart:iStop); 
    rTT = rTWL(iStart:iStop);     
    rTWL_over_yearly(iYear) = mean(rTT(rTT > rBermEl));
    rHs_over_yearly(iYear)  = mean(rHH(rTT > rBermEl));
    rR2_over_yearly(iYear)  = mean(rRR(rTT > rBermEl));
    rNTR_over_yearly(iYear) = mean(rNTR(rTT > rBermEl));
    
    iStart = iStop + 1;
    
end

% identify Hs threshold to qualify as a storm event, round nearest 0.05 m
%rHs_over_yearly(28) = NaN; % Remove year 2007 (anonymously low)?
nHs_min = min(rHs_over_yearly);
nHs_threshold = floor(nHs_min / 0.05) * 0.05     

% Alternative method: Define Hs from lower 2 sigma of all TWL exceedances
% nHs_min = mean(rHH(rTT > rBermEl)) - 2 * std(rHH(rTT > rBermEl));
% nHs_threshold = floor(nHs_min / 0.05) * 0.05

% visual check of threshold and drivers (this is hard coded to 1980, will 
% need to be changed for NC)
figure
subplot(2,1,1)
plot(1980 : 1 : 1980+nYrs-1, [rTWL_over_yearly, rR2_over_yearly, rNTR_over_yearly], '-o')
ylabel('Contribution to TWL')
legend('TWL', 'R2/%', 'NTR')
subplot(2,1,2)
plot(1980 : 1 : 1980+nYrs-1, rHs_over_yearly, '-o')
hold on
refline(0, nHs_threshold)
ylabel('Hs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find storms (translated from Ian Reeves); use corrected data only
[iStormStart, iStormStop, dtStormStart, dtStormStop, cStormHs, ...
    cStormDur, cStormTWL, cStormRlow, cStormTp, cStormNTR, cStormAT, ...
    cStormWavD, cStormNegSurgeDT, cStormNegSurgeNTR, dtYear] = deal(cell(0));

% NOTE for Ian, our storms are all significantly under these max thresholds; 
% these thresholds are important for the synthetic storms (see new bit
% of code in the copula section - deleted here)
dur_min = 8;   % hr, Magliocca et al. (2011)
dtH_yr  = year(dtH);

t = 1 ;
while t <= length(dtH)
    % discard storms where simultaneous surge is negative (will want to
    % check if we omit any, note that this omits many storms b/c of nans)
    if rHs_corr(t) >= nHs_threshold && rNTR_nan(t) >=0
        stormStart = t;
        dur = 1;
        t = t + 1;
        % If Hs drops below Hs_threshold for only 24 hrs or less, 
        % exceedence is assumed part of same weather system 
        % (Wahl et al., 2016; Li et al., 2014)
        while sum(rHs_corr(t:t+24) > nHs_threshold) > 0 
            if rHs_corr(t) > nHs_threshold
                dur = dur + 1;
                t = t + 1;
            else
                t = t + 1;
            end
        end
        % Minimum of an 8 hr storm (Magliocca et al., 2011)
        if dur > dur_min
            stormStop = t;
            iStormStart{end+1}  = stormStart;
            iStormStop{end+1}   = stormStop;
            dtStormStart{end+1} = datenum(dtH(stormStart));
            dtStormStop{end+1}  = datenum(dtH(stormStop));
            cStormDur{end+1}    = dur;
            cStormTWL{end+1}    = max(rTWL(stormStart:stormStop));
            cStormRlow{end+1}   = max(rRlow(stormStart:stormStop));
            % Ian, note that you need to find the max Hs and simultaneous
            % (not max Tp and WavD)
            [cStormHs{end+1}, iHs] = max(rHs_corr(stormStart:stormStop));
            tmpTp   = rTp_corr(stormStart:stormStop);
            tmpWavD = rWavD(stormStart:stormStop);
            cStormTp{end+1}     = tmpTp(iHs);
            cStormWavD{end+1}   = tmpWavD(iHs);
            % Ian, same here: find the max NTR and simultaneous AT,
            % otherwise you will only have positive tidal values
            [cStormNTR{end+1}, iNTR] = max(rNTR_nan(stormStart:stormStop));
            tmpAT   = rAT_nan(stormStart:stormStop);
            cStormAT{end+1}     = tmpAT(iNTR);  
            dtYear{end+1}       = dtH_yr(stormStart);
        end
        
        t = t + 1;  % KA: need to check this 
    else
        
        % for debugging - see if there are any negative surge values
        % during large wave events
        if rHs_corr(t) >= nHs_threshold && rNTR_nan(t) <0 
            
            % save datetime, wave height, surge 
            cStormNegSurgeDT{end+1}  = datenum(dtH(t));
            cStormNegSurgeNTR{end+1} = rNTR_nan(t);
            
            % takeaway: indeed, there are many large wave events with
            % negative surge or nearly zero surge
            
        end
        
        t = t + 1;

    end
    
end         

% for debugging
% figure; scatter(datetime(cell2mat(cStormNegSurgeDT), 'ConvertFrom', 'datenum'), ...
%    cell2mat(cStormNegSurgeNTR))

% convert cells back to arrays
rStormTWL = cell2mat(cStormTWL)';
rStormRlow = cell2mat(cStormRlow);
rStormHs  = cell2mat(cStormHs)';
rStormDur = cell2mat(cStormDur)';
rStormTp  = cell2mat(cStormTp)';
rStormNTR = cell2mat(cStormNTR)';
rStormAT  = cell2mat(cStormAT)';
rStormStart = cell2mat(iStormStart)';
rStormStop  = cell2mat(iStormStop)';
rStormStart_dt = cell2mat(dtStormStart)';
rStormStop_dt  = cell2mat(dtStormStop)';
rYear = cell2mat(dtYear)';

% Create matrix of all storms and parameters
[len, ~] = size(rStormTWL);
Storms = zeros(len, 12);
Storms(:,1) = rStormStart;
Storms(:,2) = rStormStop;
Storms(:,3) = rStormStart_dt;
Storms(:,4) = rStormStop_dt;
Storms(:,5) = rStormHs;
Storms(:,6) = rStormDur;
Storms(:,7) = rStormTWL;
Storms(:,8) = rStormNTR;
Storms(:,9) = rStormTp;
Storms(:,10) = rStormAT;
Storms(:,11) = rStormRlow;
Storms(:,12) = rYear;

% print number of storms
nStorms = length(rStormTWL)

% Plot storm TWL, Hs, Dur, Tp, NTR, AT histogram
figure
subplot(2,3,1)
hist(rStormTWL, 50)
ylabel('Storm TWL [m NAVD88]')
title('Berm Elev = 2 m NAVD88')  % hard-coded: will need to make more modular for NC
%plotFancyAxis

subplot(2,3,2)
hist(rStormHs, 50)
ylabel('Storm Hs [m]')
title('1980 - 2014')
%plotFancyAxis

subplot(2,3,3)
hist(rStormDur, 50)
ylabel('Storm Dur [hrs]')
%plotFancyAxis

subplot(2,3,4)
hist(rStormTp, 50)
ylabel('Storm Tp [s]')
%plotFancyAxis

subplot(2,3,5)
hist(rStormNTR, 50)
ylabel('Storm \eta_{NTR} [m]')
%plotFancyAxis

subplot(2,3,6)
hist(rStormAT, 50)
ylabel('Storm \eta_{A} [m NAVD88]')
%plotFancyAxis

%% USE COPULAS TO MODEL INTERDEPENDENCY BETWEEN VARIABLES

% From Matlab: To generate data Xsim with a distribution "just like" 
% (in terms of marginal distributions and correlations) the distribution of 
% data in the matrix X, you need to: 
%       1) fit marginal distributions to the columns of X
%       2) use appropriate cdf functions to transform X to U ([0,1] space) 
%       3) use copulafit to fit a copula to U 
%       4) generate new data Usim from the copula
%       5) use appropriate inverse cdf functions to transform Usim to Xsim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: fit marginal distributions to each variable 

% From MvCAT: struct D of fitted distributions/parameters, and struct PD  
% representing the fitted distributions (ProbDist class). 
[stD_U1, stPD_U1] = allfitdist(rStormNTR, 'PDF'); % Fit a distribution to NTR
[stD_U2, stPD_U2] = allfitdist(rStormHs, 'PDF');  % Fit a distribution to Hs
[stD_U3, stPD_U3] = allfitdist(rStormTp, 'PDF');  % Fit a distribution to Tp
[stD_U4, stPD_U4] = allfitdist(rStormDur, 'PDF'); % Fit a distribution to D
[stD_U5, stPD_U5] = allfitdist(rStormAT, 'PDF');  % Fit a distribution to AT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: use appropriate cdf functions to transform X to U ([0,1] space) %

% KA: Matlab example transforms the data to the copula scale (unit square)
% using a kernal estimator, which provide a smooth estimate of the CDF,
% however my gut (and I think what Wahl did) says to use the empirical CDF
% (below) -- try both?
rEP1 = cdf(stPD_U1{1}, rStormNTR);  % Rayleigh
rEP2 = cdf(stPD_U2{1}, rStormHs);   % Generalized Pareto
rEP3 = cdf(stPD_U3{1}, rStormTp);   % Generalized Extreme Value
rEP4 = cdf(stPD_U4{1}, rStormDur);  % Generalized Pareto
rEP5 = cdf(stPD_U5{1}, rStormAT);   % Nakagami Distribution? basically normal

% the following code produces the same result as ecdf() function 

% Find data ranks
n    = length(rStormNTR); 
data = [rStormNTR, rStormHs, rStormTp, rStormDur];
[rR1, rR2, rR3, rR4] = deal(nan(n,1));

for i = 1:n
    rR1(i,1) = sum(data(:,1) >= data(i,1));
    rR2(i,1) = sum(data(:,2) >= data(i,2));
    rR3(i,1) = sum(data(:,3) >= data(i,3));
    rR4(i,1) = sum(data(:,4) >= data(i,4));
end

% Transform to uniform marginals (the empirical CDF)
rU1 = (n-rR1+0.5)./n;
rU2 = (n-rR2+0.5)./n;
rU3 = (n-rR3+0.5)./n;
rU4 = (n-rR4+0.5)./n;

% for debugging, prove that ecdf() is the same as above and that the 
figure
ecdf(rStormNTR)
hold on
scatter(rStormNTR,rU1)
hold on
scatter(rStormNTR, rEP1)  % the CDF estimate from the marginal distribution

% Compute Kendal's Corelation Coefficient for each pair using the ECDF
% NOTE: if pval(a,b) is small (less than 0.05), then the correlation 
% rho(a,b) is significantly different from zero
[rTau, pval] = corr([rU1, rU2, rU3, rU4], 'type', 'kendall');
%[rTau, pval] = corr(data, 'type', 'kendall'); % exactly the same as above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: use copulafit to fit a copula to U 

% To account for interdependencies, fit copulas to the transformed 
% four-dimensional data sets. Copulas are great b/c we can mix various 
% marginal distributions. Wahl used elliptical copulas (Gaussian and 
% t-student); these are not capable of modeling tail dependence (that's why 
% we should try Vine, Archimedian, or EV copulas)...but they are easily 
% transformed to 4 dimensions (i.e., NTR, Hs, Tp, and Dur)

% for fitting vine copula models in R
dlmwrite('U_mssmVCR_2m.txt',[rU1, rU2, rU3, rU4],'delimiter','\t','precision',12)

% returns an estimate, rhohat, of the matrix of linear correlation 
% parameters for a gaussian and t copula, and an estimate of the dof parameter, nuhat, 
% given the data in [0,1] space
%rRhoHat = copulafit('Gaussian',[rEP1 rEP2 rEP3 rEP4]); 
%rRhoHat = copulafit('Gaussian',[rU1, rU2, rU3, rU4]);  % the tau values look better for this (the ecdf, what Wahl did)
%[rRhoHat, rNuHat, rNuCI] = copulafit('t',[rU1, rU2, rU3, rU4]); % the tau values look best (of elliptical) for this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: generate new data Usim from the copula

% Generate random samples from the copula ([0,1] space, sampled 
% from a continuous uniform distribution)
nSimNum = 10000; % Number of simulated storms to create
%rU = copularnd('Gaussian', rRhoHat, nSimNum);  % gaussian copula, elliptical
%rU = copularnd('t', rRhoHat, rNuHat, nSimNum); % t-student copula, also elliptical
                  
% load outputs from vine model in R
rU = load('Usim_mssmVCR-Cvine.mat', 'uSim');
rU = rU.uSim;  % best tau values (closest to empirical)

% kendall's coefficient for the simulated data
rTauSim = corr([rU(:,1), rU(:,2), rU(:,3), rU(:,4)],'type','Kendall');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: use appropriate inverse cdf functions to transform Usim to Xsim

% use the inverse of the fitted marginal CDFs to transform the simulated 
% data from the unit hypercube space back to the original scale of the data
rSimNTR = icdf(stPD_U1{1}, rU(:,1));
rSimHs  = icdf(stPD_U2{1}, rU(:,2));
rSimTp  = icdf(stPD_U3{1}, rU(:,3));
rSimDur = icdf(stPD_U4{1}, rU(:,4));

% this workflow simulated 10000 quadruplets of NTR, Hs, Tp, and Dur in the 
% unit hypercube (preserves the interdependencies between variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, simulate the tide randomly from its ecdf

% KA: Ian, so really they just sample randomly from the empirical cdf 
% (observations)...we were making it complicated
[f,rU5] = ecdf(rStormAT);
randAT  = randi(length(rU5), nSimNum, 1);
rSimAT  = rU5(randAT);
%[len, ~] = size(rIEP5);
%randAT = randi(len,nSimNum,1);
%rSimAT = rIEP5(randAT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot Wahl Figure 6
figure
subplot(5,6,1)
hist(rStormNTR, 50)
ylabel('\eta_{NTR} [m]')
subplot(5,6,2)
hist(rSimNTR, 50)
ylabel('\etaSim_{NTR} [m]')

subplot(5,6,7)
scatter(rSimNTR, rSimAT)
hold on
scatter(rStormNTR, rStormAT)
ylabel('\eta_{AT} [m]')
%scatterhist(rStormNTR, rStormHs, 'Direction', 'out')

subplot(5,6,13)
scatter(rSimNTR, rSimHs)
hold on
scatter(rStormNTR, rStormHs)
text(0.5,1,sprintf('t = %0.2f', rTau(2,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,1)))
ylabel('Hs [m]')

subplot(5,6,19)
scatter(rSimNTR, rSimTp)
hold on
scatter(rStormNTR, rStormTp)
text(0.5,1,sprintf('t = %0.2f', rTau(3,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(3,1)))
ylabel('Tp [s]')

subplot(5,6,25)
scatter(rSimNTR, rSimDur)
hold on
scatter(rStormNTR, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(4,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(4,1)))
ylabel('D [h]')
xlabel('\eta_{NTR} [m]')

subplot(5,6,8)
hist(rStormAT, 50)
ylabel('\eta_{AT} [m]')
subplot(5,6,9)
hist(rSimAT, 50)
ylabel('\etaSim_{AT} [m]')

subplot(5,6,14)
scatter(rSimAT, rSimHs)
hold on
scatter(rStormAT, rStormHs)
ylabel('Hs [m]')

subplot(5,6,20)
scatter(rSimAT, rSimTp)
hold on
scatter(rStormAT, rStormTp)
ylabel('Tp [s]')

subplot(5,6,26)
scatter(rSimAT, rSimDur)
hold on
scatter(rStormAT, rStormDur)
ylabel('D [h]')
xlabel('\eta_{AT} [m]')

subplot(5,6,15)
hist(rStormHs, 50)
ylabel('Hs [m]')
subplot(5,6,16)
hist(rSimHs, 50)
ylabel('HsSim [m]')

subplot(5,6,21)
scatter(rSimHs, rSimTp)
hold on
scatter(rStormHs, rStormTp)
text(0.5,1,sprintf('t = %0.2f', rTau(2,3)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,3)))
ylabel('Tp [s]')

subplot(5,6,27)
scatter(rSimHs, rSimDur)
hold on
scatter(rStormHs, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(2,4)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,4)))
ylabel('D [h]')
xlabel('Hs [m]')

subplot(5,6,22)
hist(rStormTp, 50)
ylabel('Tp [s]')
subplot(5,6,23)
hist(rSimTp, 50)
ylabel('TpSim [s]')

subplot(5,6,28)
scatter(rSimTp, rSimDur)
hold on
scatter(rStormTp, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(3,4)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(3,4)))
ylabel('D [h]')
xlabel('Tp [s]')

subplot(5,6,29)
hist(rStormDur, 50)
ylabel('D [h]')
xlabel('D [h]')
subplot(5,6,30)
hist(rSimDur, 50)
ylabel('D [h]')
xlabel('D [h]')

%% calculate simulated R2% and add to SL to get the simulated TWL

rBeta    = 0.04;                         % beach slope, Hog Island    
rSimL0   = (9.8 * rSimTp.^2) / (2 * pi); % wavelength       
rSimSetup = 0.35 * rBeta * sqrt(rSimHs .* rSimL0); 
rSimSin   = 0.75 * rBeta * sqrt(rSimHs .* rSimL0);  % incident band swash
rSimSig   = 0.06 * sqrt(rSimHs .* rSimL0) ;         % infragravity band swash
rSimSwash = sqrt((rSimSin.^2) + (rSimSig.^2));      % total swash
rSimR2    = 1.1 * (rSimSetup + (rSimSwash/2));      % R2%

rSimTWL  = rSimNTR + rSimR2 + rSimAT;       
rSimRlow = (rSimTWL - (rSimSwash/2));  

%% Lastly, apply max thresholds for synthetics

nZ = 22; % m
nHs_max  = 0.5 * sqrt(2) * nZ; % m, Thorton and Guza [1982]
nTp_max  = 30;
nDur_max = 240; % hr, Wahl et al. (2016) (NOTE, this seems arbitrary to me)

% only save synthetic storms below thresholds
% Find storms (translated from Ian Reeves); use corrected data only
[cSimHs, cSimDur, cSimTp, cSimNTR, cSimAT, cSimRlow, cSimTWL] = deal(cell(0));

for iSim = 1 : length(rSimDur)
    if rSimHs(iSim)<nHs_max && rSimTp(iSim)<nTp_max && rSimDur(iSim)<nDur_max
        
        cSimHs{end+1} = rSimHs(iSim);
        cSimTp{end+1} = rSimTp(iSim);
        cSimDur{end+1} = rSimDur(iSim);
        cSimNTR{end+1} = rSimNTR(iSim);
        cSimAT{end+1}  = rSimAT(iSim);
        cSimRlow{end+1} = rSimRlow(iSim);
        cSimTWL{end+1}  = rSimTWL(iSim);
        
    end
end

% Create matrix of all simulated storms and parameters
SimStorms = zeros(length(cSimTWL), 7);

SimStorms(:,1) = cell2mat(cSimHs);
SimStorms(:,2) = cell2mat(cSimDur);
SimStorms(:,3) = cell2mat(cSimTWL);
SimStorms(:,4) = cell2mat(cSimNTR);
SimStorms(:,5) = cell2mat(cSimTp);
SimStorms(:,6) = cell2mat(cSimAT);
SimStorms(:,7) = cell2mat(cSimRlow);

% % Plot final storm TWL, Hs, Dur, Tp, NTR, AT histogram
% figure
% subplot(2,3,1)
% hist(SimStorms(:,3), 50)
% ylabel('Simulated TWL [m]')
% title('10,000 storms')
% 
% subplot(2,3,2)
% hist(SimStorms(:,1), 50)
% ylabel('Simulated Hs [m]')
% 
% subplot(2,3,3)
% hist(SimStorms(:,2), 50)
% ylabel('Simulated Dur [hrs]')
% 
% subplot(2,3,4)
% hist(SimStorms(:,5), 50)
% ylabel('Simulated Tp [s]')

% plot Wahl Figure 6 (again)
figure
subplot(5,6,1)
hist(rStormNTR, 50)
ylabel('\eta_{NTR} [m]')
subplot(5,6,2)
hist(SimStorms(:,4), 50)
ylabel('\etaSim_{NTR} [m]')

subplot(5,6,7)
scatter(SimStorms(:,4), SimStorms(:,6))
hold on
scatter(rStormNTR, rStormAT)
ylabel('\eta_{AT} [m]')
%scatterhist(rStormNTR, rStormHs, 'Direction', 'out')

subplot(5,6,13)
scatter(SimStorms(:,4), SimStorms(:,1))
hold on
scatter(rStormNTR, rStormHs)
text(0.5,1,sprintf('t = %0.2f', rTau(2,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,1)))
ylabel('Hs [m]')

subplot(5,6,19)
scatter(SimStorms(:,4), SimStorms(:,5))
hold on
scatter(rStormNTR, rStormTp)
text(0.5,1,sprintf('t = %0.2f', rTau(3,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(3,1)))
ylabel('Tp [s]')

subplot(5,6,25)
scatter(SimStorms(:,4), SimStorms(:,2))
hold on
scatter(rStormNTR, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(4,1)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(4,1)))
ylabel('D [h]')
xlabel('\eta_{NTR} [m]')

subplot(5,6,8)
hist(rStormAT, 50)
ylabel('\eta_{AT} [m]')
subplot(5,6,9)
hist(SimStorms(:,6), 50)
ylabel('\etaSim_{AT} [m]')

subplot(5,6,14)
scatter(SimStorms(:,6), SimStorms(:,1))
hold on
scatter(rStormAT, rStormHs)
ylabel('Hs [m]')

subplot(5,6,20)
scatter(SimStorms(:,6), SimStorms(:,5))
hold on
scatter(rStormAT, rStormTp)
ylabel('Tp [s]')

subplot(5,6,26)
scatter(SimStorms(:,6), SimStorms(:,2))
hold on
scatter(rStormAT, rStormDur)
ylabel('D [h]')
xlabel('\eta_{AT} [m]')

subplot(5,6,15)
hist(rStormHs, 50)
ylabel('Hs [m]')
subplot(5,6,16)
hist(SimStorms(:,1), 50)
ylabel('HsSim [m]')

subplot(5,6,21)
scatter(SimStorms(:,1), SimStorms(:,5))
hold on
scatter(rStormHs, rStormTp)
text(0.5,1,sprintf('t = %0.2f', rTau(2,3)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,3)))
ylabel('Tp [s]')

subplot(5,6,27)
scatter(SimStorms(:,1), SimStorms(:,2))
hold on
scatter(rStormHs, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(2,4)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(2,4)))
ylabel('D [h]')
xlabel('Hs [m]')

subplot(5,6,22)
hist(rStormTp, 50)
ylabel('Tp [s]')
subplot(5,6,23)
hist(SimStorms(:,5), 50)
ylabel('TpSim [s]')

subplot(5,6,28)
scatter(SimStorms(:,5), SimStorms(:,2))
hold on
scatter(rStormTp, rStormDur)
text(0.5,1,sprintf('t = %0.2f', rTau(3,4)))
text(1.5,1, sprintf('t_{sim} = %0.2f', rTauSim(3,4)))
ylabel('D [h]')
xlabel('Tp [s]')

subplot(5,6,29)
hist(rStormDur, 50)
ylabel('D [h]')
xlabel('D [h]')
subplot(5,6,30)
hist(SimStorms(:,2), 50)
ylabel('D [h]')
xlabel('D [h]')
