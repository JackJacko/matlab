%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION A4 - Size and RI binning: concentration corrections, PSDs and 
%%% PRIDs.

%%% 4.1 Size and refractive index binning
            
%%% Size and refractive index are assigned to each particle based on 
%%% their proximity to the nodes in the FC model grid. Before binning 
%%% proper, a few correction terms are needed.

%%% 4.1.1 Sample pump flow rate correction

%%% The first correction term accounts for particle concentration 
%%% underestimation at low sample pump flow rates. This correction will 
%%% need the parameter pumpSpeed to be set by the user to be equal to 
%%% that used during the measurement. The function coeffCorrFind simply 
%%% finds the correction curve value corresponding to pumpSpeed - see 
%%% pump flow rate correction results, Chapter 3.

platCorr = 1.087628356;                                                     % Correction factor derived from flow rate test data - see pump flow rate correction results, Chapter 3
pumpSpeed = 0.5;                                                            % Pump speed in uL/s.
cCorr = coeffCorrFind(pumpSpeed);                                           % Sets the pump speed correction coefficient - more on this in the function definition and further on in this script.
coeffCorr = (1/cCorr)*platCorr;                                             % Defines the final flow rate correction coefficient

%%% 4.1.2 Sensitivity run volume correction

%%% A further correction is required to account for data transfer 
%%% overhead time. While the instrument software transfers particle 
%%% data, it cannot analyse the sample. This means not only that the 
%%% analysed volume is always less than the processed volume, but that 
%%% the analysed volume depends on the amount of detected particles 
%%% too. The more particles there are, the larger the difference 
%%% between analysed and pumped volume. In turn, this means that the 
%%% analysed volume will be different for different sensitivity 
%%% settings (high sensitivity generally means that more particles are 
%%% detected by the sensors). Individual analysed volumes are imported 
%%% and compared to produce a correction that is applied to each 
%%% sensitivity setting.

%%% The following line will require input from the user (sens_vols). It 
%%% is meant to identify the largest analysed volume among those 
%%% corresponding to each sensitivity run. These are provided directly 
%%% by the flow cytometer data analysis software (CytoClus).                        

sens_vols_max = max(sens_vols);                                             % Gets the largest analysed volume.

%%% 4.1.3 Binning

FSClin = reshape(FSC,[1,NV*RV]);                                            % Reshaping FWS array.
SSClin = reshape(SSC,[1,NV*RV]);                                            % Reshaping SWS array .
SIZElin = reshape(SIZEs,[1,NV*RV]);                                         % Reshaping diameters array.
RIslin = reshape(RIs,[1,NV*RV]);                                            % Reshaping RI array.

dataN = size(main_data_total,1);                                            % Gets the number of data points in main_data_total.
 
particleSIZEs = ones(dataN,1);                                              % Preallocates the array that will contain the sizes of each particle.
particleRIs = ones(dataN,1);                                                % Preallocates the array that will contain the RIs of each particle.
NUMb = ones(dataN,1);                                                       % A simple array of ones. Summing along this array effectively counts particles, and will be the basis of the concentration calculations.
 
for q = 1:1:dataN                                                           % This loop assigns size & RI of the closest grid point to each particle.
    [~,ind] = min(abs(FSClin-main_data_total(q,1)) + ...                             
        abs(SSClin-main_data_total(q,2)));                                  % Finds closest grid point.

    particleSIZEs(q) = SIZElin(ind);                                        % Assigns the size to the particle.
        particleRIs(q) = RIslin(ind);                                       % Assigns the RI to the particle.
end

nrsort = sort(nr);                                                          % Sorts the real refractive index values in ascending order
B1 = ones(1,numel(nrsort)-1);                                               % Preallocates the array of values used as RI bin boundaries
            
for h = 1:1:numel(nrsort)-1
    B1(h) = (nrsort(h)+nrsort(h+1))/2;                                      % Defines the values used as RI bin boundaries
end         
B1 = [nrsort(1)-((nrsort(2)-nrsort(1))/2),B1,nrsort(end)+ ...
    ((nrsort(end)-nrsort(end-1))/2)];                                       % Further redefines the first and last RI bin boundary values

B1N = numel(B1);                                                            % Gets the number of values used as RI bin boundaries
            
%%% The following line will require input from the user (size_bins). 
%%% These are log-spaced size bin boundaries.

B2 = size_bins;                                                             % The array of log-spaced values used as size bin boundaries
B2N = numel(B2);                                                            % Gets the number of log-spaced values used as size bin boundaries
            
binMAT = zeros(B1N-1,B2N-1);                                                % Preallocates the matrix that will contain the particles after being binned by RI and by size 
binMAT_FLRaye = zeros(B1N-1,B2N-1);                                         % Preallocates the matrix that will contain the FL particles after being binned by RI and by size binMAT_FLRnay = zeros(B1N-1,B2N-1);           % Preallocates the matrix that will contain the non-FL particles after being binned by RI and by size
                        
loB = 1;                                                                    % Primes the low boundary of the particles (relative to sensitivity setting).
hiB = 0;                                                                    % Primes the high boundary of the particles (relative to sensitivity setting).
            
for u = 1:1:sensS_num
    if main_data_nums(u) ~= 0                                               % Checks if there are particles corresponding to the sensitivity setting.
        if u ~= 1                                                           % Updates the low boundary with each step save the first.
            loB = loB + main_data_nums(u-1);
        end
        hiB = hiB + main_data_nums(u);                                      % Updates the high boundary.
        for tR = 1:1:B1N-1
            indRI = find(particleRIs(loB:hiB) >= B1(tR) & ...
                particleRIs(loB:hiB) < B1(tR+1));                           % Indices of the particles belonging to the RI bin
            
            pinPoint = particleSIZEs(loB:hiB);                              % Singles out particles within the sensitivity boundaries
            SEL_binMAT = pinPoint(indRI);                                   % Singles out an array containing the particles within the boundaries and belonging to the RI bin, ready to be binned by size.
            NUMbSEL = NUMb(indRI);                                          % A reduced array of ones. Summing along this array effectively counts particles.
            
            pinPoint_FLR = main_data_total_FLR(loB:hiB);                    % Singles out FL values of particles within the sensitivity boundaries
            
            SEL_binMAT_FLR_pre = pinPoint_FLR(indRI);                       % Singles out an array containing FL values of particles within the boundaries and belonging to the RI bin, ready to be binned by size.
            
            indRIFLRaye = find(SEL_binMAT_FLR_pre > 10);                    % Indices of the FL particles within the boundaries and belonging to the RI bin. 10 is the value chosen to separate FL particles from the FL noise background
            
            indRIFLRnay = find(SEL_binMAT_FLR_pre <= 10);                   % Indices of the non-FL particles within the boundaries and belonging to the RI bin.
            
            SEL_binMAT_FLRaye = SEL_binMAT(indRIFLRaye);                    % Singles out an array containing FL particles within the boundaries and belonging to the RI bin, ready to be binned by size.
            SEL_binMAT_FLRnay = SEL_binMAT(indRIFLRnay);                    % Singles out an array containing non-FL particles within the boundaries and belonging to the RI bin, ready to be binned by size.
            
            NUMbSEL_FLRaye = NUMbSEL(indRIFLRaye);                          % A reduced array of ones. Summing along this array effectively counts FL particles.
            
            NUMbSEL_FLRnay = NUMbSEL(indRIFLRnay);                          % A reduced array of ones. Summing along this array effectively counts non-FL particles.
            
            for tS = 1:1:B2N-1
                indSIZE= find(SEL_binMAT >= B2(tS) & ...
                    SEL_binMAT < B2(tS+1));                                 % Indices of the particles belonging to the size bin.
                
                indSIZE_FLRaye= find(SEL_binMAT_FLRaye >= ...
                    B2(tS) & SEL_binMAT_FLRaye < B2(tS+1));                 % Indices of FL particles belonging to the size bin.
                
                indSIZE_FLRnay= find(SEL_binMAT_FLRnay >= ...
                    B2(tS) & SEL_binMAT_FLRnay < B2(tS+1));                 % Indices of non-FL particles belonging to the size bin.
                
                binMAT_mult = sens_vols_max/sens_vols(u);                   % This is the analysed volume multiplier.
                
                binMAT(tR,tS) = binMAT(tR,tS) + ...
                    (sum(NUMbSEL(indSIZE))*binMAT_mult);                    % Modified number of elements inside the bin. Adds number of elements from previous sensitivity settings to itself.
                
                binMAT_FLRaye(tR,tS) = binMAT_FLRaye(tR,tS) + ...
                    (sum(NUMbSEL_FLRaye(indSIZE_FLRaye))*binMAT_mult);      % Modified number of FL elements inside the bin. Adds number of elements from previous sensitivity settings to itself.
                
                binMAT_FLRnay(tR,tS) = binMAT_FLRnay(tR,tS) + ...
                    (sum(NUMbSEL_FLRnay(indSIZE_FLRnay))*binMAT_mult);      % Modified number of non-FL elements inside the bin. Adds number of elements from previous sensitivity settings to itself.
            end
        end
    end
end

binMAT = (binMAT.*coeffCorr)./(sens_vols_max*(1e-3));                       % Final particle concentrations for the RI and size bins, after application of the flow rate and volume corrections. Converted to number concentration per millilitre
binMAT_FLRaye = (binMAT_FLRaye.*coeffCorr)./(sens_vols_max*(1e-3));         % Final FL particle concentrations for the RI and size bins, after application of the flow rate and volume corrections. Converted to number concentration per millilitre
binMAT_FLRnay = (binMAT_FLRnay.*coeffCorr)./(sens_vols_max*(1e-3));         % Final non-FL particle concentrations for the RI and size bins, after application of the flow rate and volume corrections. Converted to number concentration per millilitre

PSD = sum(binMAT,1);                                                        % Final total PSD
PRID = sum(binMAT,2);                                                       % Final total PRID 

FL_PSD = sum(binMAT_FLRaye,1);                                              % Final total FL PSD
FL_PRID = sum(binMAT_FLRaye,2);                                             % Final total FL PRID

nFL_PSD = sum(binMAT_FLRnay,1);                                             % Final total non-FL PSD
nFL_PRID = sum(binMAT_FLRnay,2);                                            % Final total non-FL PRID
