function [stStorms, stSimStorms] = multivariateSeaStorm(sCopula, ...
    sWIS_filename, sWaterLevel_filename, fBeta, fBermEl, nSimStorm, bPlot)
% [stStorms, stSimStorms] = multivariateSeaStorm("c-vine", ...
%     'ST63183_v03.onlns', 'Tide-8631044-Combined.txt', 0.04, 1.9, 10000, true)
%
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
%       04/30/21        KAnarde             Made more user friendly with
%                                           structures
%
% Variable naming rules:
%       r - real array
%       n - integer value
%       f - float value
%       b - boolean
%       s - string
%       dt - datetime
%
% Inputs:
%       sCopula              - options are "c-vine", "d-vine", "gaussian", 
%                              or "t"
%       sWIS_filename        - .onlns file downloaded from WIS; must  
%                              contain hourly records of wave height (m)
%       sWaterLevel_filename - must contain hourly records of total water 
%                              levels in m NAVD88 as second column, first
%                              column is datetime; downloaded from NOAA;
%                              must be either the same length or longer
%                              than WIS time record
%       fBeta                - beach slope 
%       nSimStorm            - number of simulated storms to create
%       fBermEl              - erosion threshold; Wahl used 5% of dune toe 
%                              heights; we use the average berm elevation 
%                              (m NAVD88)
%       bPlot                - boolean for plotting
%
%% -----------------------------------------------------------------------%

nDays = 365; % for running mean calculation Wahl et al., used 30 
nMedianYears = 3; % yrs, rep. of recent climate, for time series correction

% thresholds for identifying storms
nMinDuration = 8; % hr, minimum duration of a storm 
nZ = 22; % m
fHs_max  = 0.5 * sqrt(2) * nZ; % m, max wave height, Thorton & Guza [1982]
nTp_max  = 30; % sec, max peak period
nDur_max = 240; % hr, Wahl et al. (2016) (this seems arbitrary)

nStartYear = 1980; % just used for plotting
    
%% -----------------------------------------------------------------------%

% load WIS a tide data
stObs = load_data(sWIS_filename, sWaterLevel_filename);  
    
% add nan for data gap, shorten water level time series to match waves 
stObs = process_tides(stObs); 

% remove nonstationarity
stObs = remove_nonstationarity(nDays, nMedianYears, stObs);

% calculate the tidal and non-tidal residuals
stObs = tide_residuals(stObs);

% calculate the TWL and R2% from observational record
stObs = calculate_TWL(fBeta, stObs, fBermEl);

% extract sea-storm events from the observational record
stStorms = extract_sea_storms_from_obs(fBermEl, nMinDuration, stObs);

% use copulas to model interdependencies, create synthetic storms
[stStorms, stSimStorms] = mssm(sCopula, nSimStorm, stStorms);

% calculate R2% and add to SL to get the TWL for synthetic storms
stSimStorms = calculate_simulated_TWL(fBeta, stSimStorms);

% lastly, apply max thresholds for synthetics
stSimStorms = apply_max_thresholds(fHs_max, nTp_max, ...
        nDur_max, stSimStorms, stStorms);

%% -----------------------------------------------------------------------%
    
function stObs = load_data(sWIS_filename, sWaterLevel_filename)
    
    % save as structure
    stObs = struct();

    % WIS data 1980-2014 (one hour time records, time zone = ?)
    % - from Wahl et al., 2016 - used offshore wave buoy in 28 m water depth, 1980-2013
    % - for Barrier3D Reeves et al., 2021 - WIS ST63183, Virgina Coastal Reserve - 22 m water depth
    rWISdata = load(sWIS_filename);
    stObs.dtH = datetime(string(rWISdata(:,1)), 'InputFormat', 'yyyyMMddHHmmss');
    stObs.rHs = rWISdata(:,10);
    stObs.rTp = rWISdata(:,12);  % 11 = TPD, 12 = TP
    stObs.rWavD = rWISdata(:,16);

    % Water levels from tide gauge (one hour time records, doesn't have to be the same as WIS dt, 
    % time zone = ?)
    % - for Barrier3D Reeves et al., 2021 - NOAA 8631044, Watchapreague VA
    fid = fopen(sWaterLevel_filename);
    rSLdata = textscan(fid, '%s %f', 'delimiter', '\t');
    fclose(fid); 
    dtSL = rSLdata{1};
    stObs.dtSL = datetime(regexprep(dtSL, '''', ''), 'InputFormat', 'dd-MMM-y HH:mm:ss');
    stObs.rSL = rSLdata{2};
 
end

function stObs = process_tides(stObs)

    % add nan for data gap, shorten time series of SL if need to
    stObs.rSL_nan = nan(size(stObs.dtH));  % for VCR, dtH goes to 2014, dtSL 2015
    stObs.dtSL_nan = nan(size(stObs.dtH));
    [~, loc] = ismember(stObs.dtSL,stObs.dtH);  % id for mapping dtSL to dtH
    idLoc = loc(loc>0);
    stObs.rSL_nan(idLoc) = stObs.rSL(1:length(idLoc));
    stObs.dtSL_nan(idLoc) = datenum(stObs.dtSL(1:length(idLoc))); % this is for t-tide only

      % for debugging
%     if bPlot
%         figure; plot(stObs.dtSL, stObs.rSL, stObs.dtH, stObs.rSL_nan)  
%         ylabel('Sea level [mNAVD88]')
%         xlabel('time')
%         legend('full time series', 'with nan for data gap')
%     end

end

function stObs = remove_nonstationarity(nDays, nMedianYears, stObs)

    % remove non-stationarity
    N = 24 * nDays;  % for running mean

    rSL_rm = stObs.rSL_nan - movmedian(stObs.rSL_nan, N, 'omitnan');
    rHs_rm = stObs.rHs - movmedian(stObs.rHs, N, 'omitnan');
    rTp_rm = stObs.rTp - movmedian(stObs.rTp, N, 'omitnan');

    % what is the median of the last 3 years? add this median back to the 
    % corrected time series so it is representative of recent climate
    idYrs = find(stObs.dtH == (stObs.dtH(end)-calmonths(12*nMedianYears)));
    stObs.rSL_corr = rSL_rm + median(stObs.rSL_nan(idYrs:end), 'omitnan');
    stObs.rHs_corr = rHs_rm + median(stObs.rHs(idYrs:end));
    stObs.rTp_corr = rTp_rm + median(stObs.rTp(idYrs:end));

    % for debugging
%     min(stObs.rTp_corr) % these get negative with 30-day medians, set to 0 just in case
%     min(stObs.rHs_corr)
    stObs.rHs_corr(stObs.rHs_corr<0) = 0;
    stObs.rTp_corr(stObs.rTp_corr<0) = 0;

end

function stObs = tide_residuals(stObs)

    % from the corrected SL time series, use t-tide to perform a 
    % year-by-year tidal analysis

    % split data separated by NaNs into chunks
    rChunk = [stObs.rSL_corr, stObs.dtSL_nan];  % C1 = corrected SL, C2 = dt
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
    stObs.rAT = tmpChunk(:,3);
    stObs.rNTR = tmpChunk(:,4);

    % now do as before and map to dtH 
    [stObs.rNTR_nan, stObs.rAT_nan] = deal(nan(size(stObs.dtH)));    
    [~, loc] = ismember(dtSLcorr_sub,stObs.dtH);  % what is the id for mapping dtSLcorr_sub to dtH
    idLoc = loc(loc>0);
    stObs.rNTR_nan(idLoc) = stObs.rNTR(1:length(idLoc));
    stObs.rAT_nan(idLoc) = stObs.rAT(1:length(idLoc));

      % for debugging
%     if bPlot
%         figure; plot(dtSLcorr_sub, rSLcorr_sub, ... % corr SL
%                     dtSLcorr_sub, stObs.rAT, ... % pred tide
%                     dtSLcorr_sub, stObs.rNTR) % residual
%     end

end

function stObs = calculate_TWL(fBeta, stObs, fBermEl)

    % calculate R2% and add to SL to get the TWL (currently only corrected data)     

    % KA: is this Stockdon 2006 broken down into components?
    stObs.rL0   = (9.8 * stObs.rTp_corr.^2) / (2 * pi); % wavelength  
    rSetup = 0.35 * fBeta * sqrt(stObs.rHs_corr .* stObs.rL0); 
    rSin   = 0.75 * fBeta * sqrt(stObs.rHs_corr .* stObs.rL0);  % incident band swash
    rSig   = 0.06 * sqrt(stObs.rHs_corr .* stObs.rL0) ;         % infragravity band swash
    rSwash = sqrt((rSin.^2) + (rSig.^2));      % total swash
    stObs.rR2    = 1.1 * (rSetup + (rSwash/2));      % R2%

    % KA: not clear from Wahl if the TWL is the corrected SL+R2...do both?
    stObs.rTWL  = stObs.rSL_corr + stObs.rR2;       % corrected for nonstationarity
    %stObs.rTWL     = stObs.rSL_nan + stObs.rR2;    % observed
    stObs.rRlow = (stObs.rTWL - (rSwash/2));  % this is just for Ian...use observed

    if bPlot        
        figure
        % SL
        subplot(4,2,1)
        plot(stObs.dtH, stObs.rSL_nan, stObs.dtH, stObs.rSL_corr) %, dtH, movmedian(rSL_nan, N, 'omitnan'))
        ylabel('Sea level [mNAVD88]')
        legend('observed', 'corrected for non-stationarity')
        title(sprintf('Berm Elev = %d m NAVD88', fBermEl))

        % Hs
        subplot(4,2,2)
        plot(stObs.dtH, stObs.rHs, stObs.dtH, stObs.rHs_corr) %, dtH, movmedian(rHs, N, 'omitnan'))
        ylabel('Hs [m]')
        ylim([0 10])

        % Tp
        % KA: note, I think this Tp data is bad between 2013-2014 
        subplot(4,2,3)
        plot(stObs.dtH, stObs.rTp, stObs.dtH, stObs.rTp_corr) %, dtH, movmedian(rTp, N, 'omitnan'))
        ylabel('Tp [s]')
        ylim([0 20])

        % WavD
        subplot(4,2,4)
        plot(stObs.dtH, stObs.rWavD)
        ylabel('Wave Direction [degree]')

        % rNTR
        subplot(4,2,5)
        plot(stObs.dtH, stObs.rNTR_nan)
        %plot(stObs.dtSLcorr_sub, stObs.rNTR)
        ylabel('\eta_{NTR} [m]')

        % rAT
        subplot(4,2,6)
        plot(stObs.dtH, stObs.rAT_nan)
        %plot(dtSLcorr_sub, rAT)
        ylabel('\eta_{A} [m NAVD88]')
        
        % R2
        subplot(4,2,7)
        plot(stObs.dtH, stObs.rR2)
        ylabel('R2/% [m]')

        % TWL
        subplot(4,2,8)
        plot(stObs.dtH, stObs.rTWL)
        ylabel('TWL [m NAVD88]')
    end
end

function stStorms = extract_sea_storms_from_obs(fBermEl,nMinDuration,stObs)

    %%% First, define TWL threshold exceedence. Then, how much the different 
    %%% variables contributed to those events and if there is a dominant driver 
    %%% that can be used for the event selection

    % for each year, find when the TWL exceeds an erosion threshold

    % find the annual average TWL from all threshold exceedances from a 
    % given year calculate annual averages of MSL (here 35-day mean), 
    % tidal amplitude, residual, and R2% during the TWL exceedances
    nYrs = floor(days(stObs.dtH(end) - stObs.dtH(1))/ 365);
    iStart = 1;
    [rHs_over_yearly, rR2_over_yearly, rNTR_over_yearly, rTWL_over_yearly] = ...
        deal(NaN(nYrs,1));

    % use only corrected data
    for iYear = 1 : nYrs

        iStop = iYear * 24 * 365;

        rHH = stObs.rHs_corr(iStart:iStop);      
        rRR = stObs.rR2(iStart:iStop); 
        rNN = stObs.rNTR_nan(iStart:iStop); 
        rTT = stObs.rTWL(iStart:iStop);     
        rTWL_over_yearly(iYear) = mean(rTT(rTT > fBermEl));
        rHs_over_yearly(iYear)  = mean(rHH(rTT > fBermEl));
        rR2_over_yearly(iYear)  = mean(rRR(rTT > fBermEl));
        rNTR_over_yearly(iYear) = mean(rNN(rTT > fBermEl));

        iStart = iStop + 1;

    end

    % identify Hs threshold to qualify as a storm event, round nearest 0.05 m
    %rHs_over_yearly(28) = NaN; % Remove year 2007 (anonymously low)?
    nHs_min = min(rHs_over_yearly);
    nHs_threshold = floor(nHs_min / 0.05) * 0.05;     

    % Alternative method: Define Hs from lower 2 sigma of all TWL exceedances
    % nHs_min = mean(rHH(rTT > rBermEl)) - 2 * std(rHH(rTT > rBermEl));
    % nHs_threshold = floor(nHs_min / 0.05) * 0.05

    % visual check of threshold and drivers (this is hard coded to 1980, will 
    % need to be changed for other locations)
    if bPlot
        figure
        subplot(2,1,1)
        plot(nStartYear : 1 : nStartYear+nYrs-1, [rTWL_over_yearly, rR2_over_yearly, rNTR_over_yearly], '-o')
        ylabel('Contribution to TWL')
        legend('TWL', 'R2/%', 'NTR')
        subplot(2,1,2)
        plot(nStartYear : 1 : nStartYear+nYrs-1, rHs_over_yearly, '-o')
        hold on
        refline(0, nHs_threshold)
        ylabel('Hs')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find storms (translated from Ian Reeves); use corrected data only
    [iStormStart, iStormStop, dtStormStart, dtStormStop, cStormHs, ...
        cStormDur, cStormTWL, cStormRlow, cStormTp, cStormNTR, cStormAT, ...
        cStormWavD, cStormNegSurgeDT, cStormNegSurgeNTR, dtYear] = deal(cell(0));

    % VCR storms are all significantly under these max thresholds
    dtH_yr = year(stObs.dtH);
    t = 1 ;

    while t <= length(stObs.dtH)

        % discard storms where simultaneous surge is negative (will want to
        % check if we omit any, note that this omits many storms b/c of nans)
        if stObs.rHs_corr(t) >= nHs_threshold && stObs.rNTR_nan(t) >=0
            stormStart = t;
            dur = 1;
            t = t + 1;

            % if Hs drops below Hs_threshold for only 24 hrs or less, 
            % exceedence is assumed part of same weather system 
            % (Wahl et al., 2016; Li et al., 2014)
            while sum(stObs.rHs_corr(t:t+24) > nHs_threshold) > 0 
                if stObs.rHs_corr(t) > nHs_threshold
                    dur = dur + 1;
                    t = t + 1;
                else
                    t = t + 1;
                end
            end

            % minimum of an 8 hr storm (Magliocca et al., 2011)
            if dur > nMinDuration
                stormStop = t;
                iStormStart{end+1} = stormStart;
                iStormStop{end+1} = stormStop;
                dtStormStart{end+1} = datenum(stObs.dtH(stormStart));
                dtStormStop{end+1} = datenum(stObs.dtH(stormStop));
                cStormDur{end+1} = dur;
                cStormTWL{end+1} = max(stObs.rTWL(stormStart:stormStop));
                cStormRlow{end+1} = max(stObs.rRlow(stormStart:stormStop));

                % need to find the max Hs and simultaneous
                % (not max Tp and WavD)
                [cStormHs{end+1}, iHs] = max(stObs.rHs_corr(stormStart:stormStop));
                tmpTp = stObs.rTp_corr(stormStart:stormStop);
                tmpWavD = stObs.rWavD(stormStart:stormStop);
                cStormTp{end+1} = tmpTp(iHs);
                cStormWavD{end+1} = tmpWavD(iHs);

                % find the max NTR and simultaneous AT,
                % otherwise you will only have positive tidal values
                [cStormNTR{end+1}, iNTR] = max(stObs.rNTR_nan(stormStart:stormStop));
                tmpAT = stObs.rAT_nan(stormStart:stormStop);
                cStormAT{end+1} = tmpAT(iNTR);  
                dtYear{end+1} = dtH_yr(stormStart);

            end

            t = t + 1;  
        else

            % for debugging - see if there are any negative surge values
            % during large wave events
            if stObs.rHs_corr(t) >= nHs_threshold && stObs.rNTR_nan(t) <0 

                % save datetime, wave height, surge 
                cStormNegSurgeDT{end+1} = datenum(stObs.dtH(t));
                cStormNegSurgeNTR{end+1} = stObs.rNTR_nan(t);

                % takeaway: indeed, there are many large wave events with
                % negative surge or nearly zero surge

            end

            t = t + 1;

        end

    end         

    % for debugging
    % figure; scatter(datetime(cell2mat(cStormNegSurgeDT), 'ConvertFrom', 'datenum'), ...
    %    cell2mat(cStormNegSurgeNTR))

    % create structure of storm parameters
    stStorms = struct();
    stStorms.rStart = cell2mat(iStormStart)';
    stStorms.rStop = cell2mat(iStormStop)';
    stStorms.rStart_dt = cell2mat(dtStormStart)';
    stStorms.rStop_dt = cell2mat(dtStormStop)';
    stStorms.rHs = cell2mat(cStormHs)';
    stStorms.rDur = cell2mat(cStormDur)';
    stStorms.rTWL = cell2mat(cStormTWL)';
    stStorms.rNTR = cell2mat(cStormNTR)';
    stStorms.rTp = cell2mat(cStormTp)';
    stStorms.rAT = cell2mat(cStormAT)';
%     stStorms.rRlow = cell2mat(cStormRlow)';
    stStorms.rYear = cell2mat(dtYear)';
    stStorms.nHs_threshold = nHs_threshold;
    stStorms.nStorms = length(stStorms.rTWL);

%     if bPlot
%         % Plot storm TWL, Hs, Dur, Tp, NTR, AT histogram
%         figure
%         subplot(2,3,1)
%         hist(stStorms.rTWL, 50)
%         ylabel('Storm TWL [m NAVD88]')
%         title(sprintf('Berm Elev = %d m NAVD88', fBermEl)) 
% 
%         subplot(2,3,2)
%         hist(stStorms.rHs, 50)
%         ylabel('Storm Hs [m]')
%         %title('1980 - 2014')
% 
%         subplot(2,3,3)
%         hist(stStorms.rDur, 50)
%         ylabel('Storm Dur [hrs]')
% 
%         subplot(2,3,4)
%         hist(stStorms.rTp, 50)
%         ylabel('Storm Tp [s]')
% 
%         subplot(2,3,5)
%         hist(stStorms.rNTR, 50)
%         ylabel('Storm \eta_{NTR} [m]')
% 
%         subplot(2,3,6)
%         hist(stStorms.rAT, 50)
%         ylabel('Storm \eta_{A} [m NAVD88]')
%     end

end

function [stStorms, stSimStorms] = mssm(sCopula, nSimStorm, stStorms)

    % From Matlab: To generate data Xsim with a distribution "just like" 
    % (in terms of marginal distributions and correlations) the dist of
    % data in the matrix X, you need to: 
    %    1) fit marginal distributions to the columns of X
    %    2) use appropriate cdf functions to transform X to U ([0,1] space) 
    %    3) use copulafit to fit a copula to U 
    %    4) generate new data Usim from the copula
    %    5) use appropriate inverse cdf functions to transform Usim to Xsim

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: fit marginal distributions to each variable 

    % From MvCAT: struct D of fitted distributions/parameters, and struct 
    % PD representing the fitted distributions (ProbDist class). 
    [~, stPD_U1] = allfitdist(stStorms.rNTR, 'PDF'); % Fit a dist to NTR
    [~, stPD_U2] = allfitdist(stStorms.rHs, 'PDF');  % Fit a dist to Hs
    [~, stPD_U3] = allfitdist(stStorms.rTp, 'PDF');  % Fit a dist to Tp
    [~, stPD_U4] = allfitdist(stStorms.rDur, 'PDF'); % Fit a dist Duration
%     [~, stPD_U5] = allfitdist(stStorms.rAT, 'PDF');  % Fit a dist to AT

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: use appropriate cdf functions to transform X->U ([0,1] space) 

    % Matlab example transforms the data to the copula scale (unit square)
    % using a kernal estimator, which provide a smooth estimate of the CDF,
    % however my gut (and I think what Wahl did) says to use the empirical 
    % CDF (below) -- try both?
%     rEP1 = cdf(stPD_U1{1}, stStorms.rNTR);  % Rayleigh
%     rEP2 = cdf(stPD_U2{1}, stStorms.rHs);   % Generalized Pareto
%     rEP3 = cdf(stPD_U3{1}, stStorms.rTp);   % Generalized Extreme Value
%     rEP4 = cdf(stPD_U4{1}, stStorms.rDur);  % Generalized Pareto
%     rEP5 = cdf(stPD_U5{1}, stStorms.rAT);   % Nakagami Distribution? 

    % The following code produces the same result as ecdf() function! 
    % Find data ranks
    n    = length(stStorms.rNTR); 
    data = [stStorms.rNTR, stStorms.rHs, stStorms.rTp, stStorms.rDur];
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

%     % for debugging, prove that ecdf() is the same as above
%     if bPlot
%         figure
%         ecdf(stStorms.rNTR)
%         hold on
%         scatter(stStorms.rNTR,rU1)
%         hold on
%         scatter(stStorms.rNTR, rEP1) % the CDF estimate from marginal dist
%     end

    % Compute Kendal's Corelation Coefficient for each pair using the ECDF
    % NOTE: if pval(a,b) is small (less than 0.05), then the correlation 
    % rho(a,b) is significantly different from zero
    [stStorms.rTau, ~] = corr([rU1, rU2, rU3, rU4], 'type', 'kendall');
    %[rTau, pval] = corr(data, 'type', 'kendall'); % same as above

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3: use copulafit to fit a copula to U 

    % To account for interdependencies, fit copulas to the transformed 
    % four-dimensional data sets. Copulas are great b/c we can mix various 
    % marginal distributions. Wahl used elliptical copulas (Gaussian and 
    % t-student); these are not capable of modeling tail dependence (that's 
    % why we should try Vine, Archimedian, or EV copulas)...but they are 
    % easily transformed to 4 dimensions (i.e., NTR, Hs, Tp, and Dur)

    % for fitting vine copula models in R
    if sCopula == "c-vine" || sCopula == "d-vine"
        dlmwrite('U_mssm.txt',[rU1, rU2, rU3, rU4],'delimiter','\t',...
            'precision',12)
        dlmwrite('Inputs_mssm.txt',nSimStorm,'delimiter',...
            '\t','precision',12)
        system('/Users/KatherineAnardeWheels/PyCharmProjects/Barrier3D/Tools/Multivariate_Sea_Storm_Model/mssmVines.R');

    elseif sCopula == "gaussian"
        % returns an estimate, rhohat, of the matrix of linear correlation 
        % parameters for a gaussian copula (elliptical), and an estimate of 
        % the dof parameter, nuhat, given the data in [0,1] space
        %rRhoHat = copulafit('Gaussian',[rEP1 rEP2 rEP3 rEP4]); 
        rRhoHat = copulafit('Gaussian',[rU1, rU2, rU3, rU4]);  
        % the tau values look better using the ecdf (what Wahl did)

    elseif sCopula == "t" 
        % now t-student copula, also elliptical
        [rRhoHat, rNuHat, ~] = copulafit('t',[rU1, rU2, rU3, rU4]); 
        % the tau values look best for t-student, vs gaussian, for VCR data

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4: generate new data Usim from the copula

    % Generate random samples from the copula ([0,1] space, sampled 
    % from a continuous uniform distribution)

    if sCopula == "c-vine"
        % load outputs from vine model in R
        rU = load('Usim_mssm-Cvine.mat', 'uSim');
        rU = rU.uSim;  % best tau values (closest to empirical for VCR)

    elseif sCopula == "d-vine"
        rU = load('Usim_mssm-Dvine.mat', 'uSim');
        rU = rU.uSim;  

    elseif sCopula == "gaussian"
        rU = copularnd('Gaussian', rRhoHat, nSimStorm); 

    elseif sCopula == "t"
        rU = copularnd('t', rRhoHat, rNuHat, nSimStorm); 
    end

    % kendall's coefficient for the simulated data
    rTauSim = corr([rU(:,1), rU(:,2), rU(:,3), rU(:,4)],'type','Kendall');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 5: use appropriate inverse cdf functions to transform Usim->Xsim

    % use the inverse of the fitted marginal CDFs to transform simulated 
    % data from the unit hypercube space back to original scale of data
    rSimNTR = icdf(stPD_U1{1}, rU(:,1));
    rSimHs  = icdf(stPD_U2{1}, rU(:,2));
    rSimTp  = icdf(stPD_U3{1}, rU(:,3));
    rSimDur = icdf(stPD_U4{1}, rU(:,4));

    % this workflow simulated 10,000 quadruplets of NTR, Hs, Tp, and Dur in 
    % the unit hypercube (preserves the interdependencies btwn variables)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Now, simulate the tide randomly from its ecdf and save everything to
    % a structure
    [~,rU5] = ecdf(stStorms.rAT);
    randAT  = randi(length(rU5), nSimStorm, 1);
    rSimAT  = rU5(randAT);

    stSimStorms = struct();
    stSimStorms.rAT = rSimAT;
    stSimStorms.rNTR = rSimNTR;   
    stSimStorms.rHs = rSimHs;
    stSimStorms.rTp = rSimTp;    
    stSimStorms.rDur = rSimDur;
    stSimStorms.rTau = rTauSim;

end

function stSimStorms = calculate_simulated_TWL(fBeta, stSimStorms)

    rSimL0   = (9.8 * stSimStorms.rTp.^2) / (2 * pi); % wavelength       
    rSimSetup = 0.35 * fBeta * sqrt(stSimStorms.rHs .* rSimL0); 
    rSimSin   = 0.75 * fBeta * sqrt(stSimStorms.rHs .* rSimL0); % incident 
    rSimSig   = 0.06 * sqrt(stSimStorms.rHs .* rSimL0) ;  % infragravity
    rSimSwash = sqrt((rSimSin.^2) + (rSimSig.^2)); % total swash

    stSimStorms.rR2    = 1.1 * (rSimSetup + (rSimSwash/2)); % R2%
    stSimStorms.rTWL  = stSimStorms.rNTR + stSimStorms.rR2 + stSimStorms.rAT;       
    stSimStorms.rRlow = (stSimStorms.rTWL - (rSimSwash/2));  

end

function stSimStorms = apply_max_thresholds(fHs_max, nTp_max, ...
        nDur_max, stSimStorms, stStorms)

    % only save synthetic storms below thresholds
    [cSimHs, cSimDur, cSimTp, cSimNTR, cSimAT, cSimRlow, cSimTWL, ...
        cSimR2] = deal(cell(0));

    for iSim = 1 : length(stSimStorms.rDur)
        if stSimStorms.rHs(iSim) < fHs_max && ...
                stSimStorms.rTp(iSim) < nTp_max && ...
                stSimStorms.rDur(iSim) < nDur_max

            cSimHs{end+1} = stSimStorms.rHs(iSim);
            cSimTp{end+1} = stSimStorms.rTp(iSim);
            cSimDur{end+1} = stSimStorms.rDur(iSim);
            cSimNTR{end+1} = stSimStorms.rNTR(iSim);
            cSimAT{end+1}  = stSimStorms.rAT(iSim);
            cSimRlow{end+1} = stSimStorms.rRlow(iSim);
            cSimTWL{end+1}  = stSimStorms.rTWL(iSim);
            cSimR2{end+1}  = stSimStorms.rR2(iSim);

        end
    end

    % save new variables
    stSimStorms.rHs = cell2mat(cSimHs)';
    stSimStorms.rDur = cell2mat(cSimDur)';
    stSimStorms.rTWL = cell2mat(cSimTWL)';
    stSimStorms.rNTR = cell2mat(cSimNTR)';
    stSimStorms.rTp = cell2mat(cSimTp)';
    stSimStorms.rAT = cell2mat(cSimAT)';
    stSimStorms.rRlow = cell2mat(cSimRlow)';
    stSimStorms.rR2 = cell2mat(cSimR2)';

    if bPlot
        % plot Wahl Figure 6
        figure
        subplot(5,6,1)
        hist(stStorms.rNTR, 50)
        ylabel('\eta_{NTR} [m]')
        subplot(5,6,2)
        hist(stSimStorms.rNTR, 50, 'FaceColor', 'green')
        ylabel('\eta_{NTR} [m]')
        legend('sim')
        title(sprintf('%s copula', sCopula))

        subplot(5,6,7)
        scatter(stSimStorms.rNTR, stSimStorms.rAT)
        hold on
        scatter(stStorms.rNTR, stStorms.rAT)
        ylabel('\eta_{AT} [m]')

        subplot(5,6,13)
        scatter(stSimStorms.rNTR, stSimStorms.rHs)
        hold on
        scatter(stStorms.rNTR, stStorms.rHs)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(2,1)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(2,1)))
        ylabel('Hs [m]')

        subplot(5,6,19)
        scatter(stSimStorms.rNTR, stSimStorms.rTp)
        hold on
        scatter(stStorms.rNTR, stStorms.rTp)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(3,1)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(3,1)))
        ylabel('Tp [s]')

        subplot(5,6,25)
        scatter(stSimStorms.rNTR, stSimStorms.rDur)
        hold on
        scatter(stStorms.rNTR, stStorms.rDur)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(4,1)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(4,1)))
        ylabel('D [h]')
        xlabel('\eta_{NTR} [m]')

        subplot(5,6,8)
        hist(stStorms.rAT, 50)
        ylabel('\eta_{AT} [m]')
        subplot(5,6,9)
        hist(stSimStorms.rAT, 50)
        ylabel('\eta_{AT} [m]')
        legend('sim')

        subplot(5,6,14)
        scatter(stSimStorms.rAT, stSimStorms.rHs)
        hold on
        scatter(stStorms.rAT, stStorms.rHs)
        ylabel('Hs [m]')

        subplot(5,6,20)
        scatter(stSimStorms.rAT, stSimStorms.rTp)
        hold on
        scatter(stStorms.rAT, stStorms.rTp)
        ylabel('Tp [s]')

        subplot(5,6,26)
        scatter(stSimStorms.rAT, stSimStorms.rDur)
        hold on
        scatter(stStorms.rAT, stStorms.rDur)
        ylabel('D [h]')
        xlabel('\eta_{AT} [m]')

        subplot(5,6,15)
        hist(stStorms.rHs, 50)
        ylabel('Hs [m]')
        subplot(5,6,16)
        hist(stSimStorms.rHs, 50)
        ylabel('Hs [m]')
        legend('sim')

        subplot(5,6,21)
        scatter(stSimStorms.rHs, stSimStorms.rTp)
        hold on
        scatter(stStorms.rHs, stStorms.rTp)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(2,3)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(2,3)))
        ylabel('Tp [s]')

        subplot(5,6,27)
        scatter(stSimStorms.rHs, stSimStorms.rDur)
        hold on
        scatter(stStorms.rHs, stStorms.rDur)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(2,4)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(2,4)))
        ylabel('D [h]')
        xlabel('Hs [m]')

        subplot(5,6,22)
        hist(stStorms.rTp, 50)
        ylabel('Tp [s]')
        subplot(5,6,23)
        hist(stSimStorms.rTp, 50)
        ylabel('Tp [s]')
        legend('sim')

        subplot(5,6,28)
        scatter(stSimStorms.rTp, stSimStorms.rDur)
        hold on
        scatter(stStorms.rTp, stStorms.rDur)
        text(0.5,1,sprintf('t = %0.2f', stStorms.rTau(3,4)))
        text(1.5,1, sprintf('t_{sim} = %0.2f', stSimStorms.rTau(3,4)))
        ylabel('D [h]')
        xlabel('Tp [s]')

        subplot(5,6,29)
        hist(stStorms.rDur, 50)
        ylabel('D [h]')
        xlabel('D [h]')
        subplot(5,6,30)
        hist(stSimStorms.rDur, 50)
        ylabel('D [h]')
        xlabel('D [h]')
        legend('sim')
    end
end

end