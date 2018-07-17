%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION B1 - PSD power law best fit and PSD extrapolation.

%%% 1.1 PSD extrapolation

%%% This section of the code defines the PSD extrapolations and 
%%% calculates the slope of the size distribution. It requires the 
%%% user to input the total PSD and the median particle size of each 
%%% size bin (plotSIZEs). Note that plotSIZEs must already contain 
%%% the bins that will correspond to the extrapolations. PSD must 
%%% accordingly contain empty (0) values in correspondence with 
%%% these bins.

[~,indminA] = max(PSD);                                                     % This index defines the particle size at which the PSD peaks before falling off due to size detection limit
[~,indmaxA] = min(abs(plotSIZEs-2e-5));                                     % This index corresponds to 20 micron, the limit of statistical significance in a majority of UKCW PSDs        

%%% The next two lines will define derived indices. These further 
%%% neglect two size bins on each end of the PSD to avoid boundary 
%%% effects     
indminAf = indminA+2;
indmaxAf = indmaxA-2;

%%% The next two lines will define a reduced PSD using the two 
%%% indices just defined.
reducedPSD = PSD*0;
reducedPSD(indminAf:indmaxAf) = PSD(indminAf:indmaxAf);

%%% The next lines define the arrays over which fitting takes place. 
%%% A further input from the user is required (plotBINs). This 
%%% contains the width of the size bins, and it will be required to 
%%% calculate the density function of the PSD. This is the one for 
%%% which distribution slope is calculated. 
fitmeXa = plotSIZEs(indminAf:indmaxAf);
fitmeBa = plotBINs(indminAf:indmaxAf);
fitmeYa = PSD(indminAf:indmaxAf);

%%% The next line defines the number of bins considered
binNum = numel(plotBINs);

expungeInd = find(fitmeYa == 0);                                            % Finds the zeroes corresponding to the yet-to-be-calculated extrapolations in preparation for the best fit calculations 

%%% The zeroes are subsequently eliminated from the fitting arrays
fitmeXa(expungeInd) = [];
fitmeBa(expungeInd) = [];
fitmeYa(expungeInd) = [];
        
pFITA = coeffvalues(fit(log(fitmeXa)',log(fitmeYa)','poly1'));              % Power law best fit of the PSD. This one will be used for the PSD extrapolations proper. Note that the fit is actually carried out as a linear fit of the logarithm of the PSD and bin values as to limit the excessive influence of the largest values on the fit
        
PDF = fitmeYa./fitmeBa;                                                     % Density function of the PSD
pFITnA = coeffvalues(fit(log(fitmeXa)',log(PDF)','poly1'));                 % Power law best fit of the density function, similarly carried out as a linear fit of log values

slope = -pFITnA(1);                                                         % The slope of the PSD, given as the slope of the underlying density function

%%% The following lines calculate the PSD extrapolations proper.
        
N = ((pFITA(1).*log(plotSIZEs))+pFITA(2));                                  % Linear extrapolation
jungePSD = exp(N);                                                          % Resulting power law extrapolation

extPSD = reducedPSD;                                                        % Preallocates extended PSD array

%%% The next line preallocates the complementary PSD array – this 
%%% contains ONLY the extensions. This is needed because IOP 
%%% calculation for the PSD will be made using the size & RI  
%%% concentration matrix rather than the single PSD. This avoids 
%%% averaging and maintains full details of the organic/inorganic 
%%% content of the particulate. On the other hand, the extensions 
%%% are NOT a matrix, but an array (which uses rRI approximations). 
%%% Calculations must be done separately and then combined.
complPSD = PSD*0;		 
        
for v = 1:1:binNum                                                          % This loop fills in complPSD and extPSD
    if reducedPSD(v) == 0	
        complPSD(v) = jungePSD(v);
        extPSD(v) = jungePSD(v);
    end
end
