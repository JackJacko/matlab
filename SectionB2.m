%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION B2 - Total particulate IOPs calculations.

%%% 2.1 FC PSD

%%% This section calculates the IOPs limitedly to the FC PSD, that 
%%% is the unextended size distribution. In practice, calculations 
%%% are made for the size & RI concentration matrix (PSDmat), which 
%%% maintains all details about composition of the particulate and 
%%% allows for organic and inorganic IOP calculations

N = PSDmat;                                                                 % The size & RI concentration matrix
N = N*1e6;                                                                  % This is a cm3 to m3 conversion

D = plotSIZEs;                                                              % Median size of the bins
RIs = plotRIs;                                                              % Median RI of the bins. This needs to be input by the user, consistently with PSDmat

%%% The following contains the RI approximations in the PSD
%%% extrapolations. The lines given here are for Mode B (see Chapter 
%%% 7). The code checks for the first and last non-zero elements of 
%%% the PSD and uses those as starting and ending point for a 4-bin 
%%% average of RI (i.e. average over first non-zero element and next 
%%% three bins, and average over last non-zero element and previous 
%%% three bins)

multRIs1 = repmat(RIs',1,65);                                               % Constructs a RI matrix to superimpose the size & RI matrix
RI_avgMe = multRIs1.*N;                                                     % Multiplies matrices in preparation for weighted average
avgRI = sum(RI_avgMe(:))/sum(N(:));                                         % Weighted average over entire PSD

check = 0;                                                                  % Primes non-zero element locator
for nz = 1:1:size(N,2)                                                      % Cycles through PSD matrix
    check = check+sum(N(:,nz));                                             % Sums over PSD matrix column
    if check ~= 0                                                           % If “check” becomes non-zero...
        nz1 = nz;
        break                                                               % ... break from loop
    end
end
check = 0;                                                                  % Primes non-zero element locator
for nz = 1:1:size(N,2)                                                      % Cycles through PSD matrix
    check = check+sum(N(:,end-nz+1));                                       % Sums over PSD matrix column
    if check ~= 0                                                           % If “check” becomes non-zero...
        nz2 = size(N,2)-nz+1;
        break                                                               % ... break from loop
    end
end

multRIs2 = repmat(RIs',1,4);                                                % Constructs a RI matrix to superimpose the size & RI matrix
N1 = N(:,nz1:nz1+3);                                                        % Isolates the section of the size & RI matrix needed for the first averaging
N2 = N(:,nz2-3:nz2);                                                        % Isolates the section of the size & RI matrix needed for the second averaging
RI_avgMe1 = multRIs2.*N1;                                                   % Multiplies matrices in preparation for weighted average
RI_avgMe2 = multRIs2.*N2;                                                   % Multiplies matrices in preparation for weighted average
avgRI1 = sum(RI_avgMe1(:))/sum(N1(:));                                      % Weighted average
avgRI2 = sum(RI_avgMe2(:))/sum(N2(:));                                      % Weighted average

RV = numel(D);                                                              % Number of elements in the diameter array
RIV = numel(RIs);                                                           % Number of elements in the RI array

%%% 2.1.1 FC PSD Mie calculations

%%% This section will manage the Mie forward modelling calculations. 
%%% The first step is to calculate the single-particle VSF and
%%% absorption corresponding to each size/RI bin of the size &
%%% refractive index matrix. Although in this appendix the
%%% corresponding code lines are given within the rest of the code, 
%%% it is advisable to conduct these calculations in a separate
%%% script and to store the resulting data in a file, since these 
%%% don’t change if the RI and size arrays remain the same. The
%%% resulting VSFs and absorptions can then be simply multiplied by 
%%% the number of elements within each bin to determine the total 
%%% particulate VSF and absorption. The angular resolution of the 
%%% VSF is defined here similarly to Appendix A.

%%% Angular ranges

%%% The next few lines define the scattering angles array. The
%%% resolution is higher at small angles.

theta1 = (0:0.01:1);                                                        % Angle range 1.
theta2 = (1.1:0.1:10);                                                      % Angle range 2.
theta3 = (11:1:180);                                                        % Angle range 3.

theta = [theta1 theta2 theta3];                                             % Scattering angles array.
thetalength = numel(theta);                                                 % Gets the number of scattering angles.
theta_rad = theta./180*pi;                                                  % Converts the scattering angles array into radians.

[~,index1] = min(abs(theta-90));                                            % Singles out the element in the scattering angles array corresponding to pi/2 for integration purposes.
[~,index2] = min(abs(theta-180));                                           % Singles out the element in the scattering angles array corresponding to pi for integration purposes.

dtheta = ones(1,thetalength);                                               % Preallocates the array which contains the increments between angles in the scattering angles array. For integration purposes.

for S1 = 1:1:(thetalength-2)                                                % This small loop fills the angle increments array.
    dtheta(S1+1) = (theta_rad(S1+2)-theta_rad(S1))/2;
end

dtheta(1) = dtheta(2);                                                      % Defines the first element of dtheta.
dtheta(thetalength) = dtheta(thetalength-1);                                % Defines the last element of dtheta.

sintheta = sin(theta_rad);                                                  % The sines of the scattering angles.

r = D/2;                                                                    % Array of log-spaced virtual particle radii.

Beta_ind_wvl = cell(RIV,RV);                                                % Preallocates the cell matrix which will contain all single-particle VSFs corresponding to each size/RI bin of the size & refractive index matrix for a certain wavelength. In this study this was chosen as 532 nm (see Section 7.3)
abs_ind_wvl = ones(RIV,RV);                                                 % Preallocates the matrix which will contain all single-particle absorption corresponding to each size/RI bin of the size & refractive index matrix for the selected wavelength.

%%% Wavelength is defined here
lambda_air = 532*1e-9;                                                      % Chosen wavelength in air (nanometres).
wnr = 1.333;                                                                % Absolute refractive index of water
lambda = lambda_air/wnr;                                                    % Chosen wavelength in water.

for RI1 = 1:1:RIV                                                           % Iterates through RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        
        %%% This part will require the user to input the imaginary
        %%% refractive indices for the selected wavelength. In this study
        %%% these were adapted from Babin et al. (2003).
        if RIs(RI1) < 1.1
            ni = niOrg_wvl;
        elseif RIs(RI1) >= 1.1
            ni = niMin_wvl;
        end
        m=RIs(RI1)+(1i*ni);                                                 % Complex refractive index.
        w=2*pi/lambda;                                                      % Wavenumber.
        x=w*r(R1);                                                          % Size parameter.
        [S1,S2,Qb,Qc,Qbb] = fastmie(x,m,[],theta_rad);                      % The core of the calculation is handled by the FASTMie script (by W. H. Slade, 2006).
        i1=abs(S1).^2;                                                      % Scattered intensity functions.                                                                  
        i2=abs(S2).^2;                                                      % Scattered intensity functions.
        Beta_ind_wvl{RI1,R1}= ((1/w)^2)*0.5*((i1+i2));                      % Single-particle VSF for diameter r(R1) and refractive index RIs(RI1)
        abs_ind_wvl(RI1,R1) = (pi/4)*(D(R1).^2).*(Qc-Qb);                   % Single-particle absorption for diameter r(R1) and refractive index RIs(RI1)
    end
end

Beta_ind = cell(RIV,RV);                                                    % Preallocates the cell matrix which will contain all total particulate VSFs corresponding to each size/RI bin of the size & refractive index matrix for the selected wavelength.
abs_ind = ones(RIV,RV);                                                     % Preallocates the cell matrix which will contain all total particulate absorptions corresponding to each size/RI bin of the size & refractive index matrix for the selected wavelength.

for RI1 = 1:1:RIV                                                           % Iterates through RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        Beta_ind{RI1,R1}= Beta_ind_wvl{RI1,R1}.*N(RI1,R1);                  % total particulate VSF for diameter r(R1) and refractive index RIs(RI1)
        abs_ind(RI1,R1) = abs_ind_wvl(RI1,R1).*N(RI1,R1);                   % total particulate absorption for diameter r(R1) and refractive index RIs(RI1)
    end
end

Beta_tot = 0;                                                               % Primes the total particulate VSF (all sizes and RIs)
for RI1 = 1:1:RIV
    for R1 = 1:1:RV
        Beta_tot = Beta_tot + Beta_ind{RI1,R1};                             % Total VSF for this wavelength
    end
end
tsd = 2*pi*dtheta.*sintheta.*Beta_tot';                                     % Total particulate VSF integral

%%% 2.1.2 FC PSD IOP calculations
%%% This section handles the calculations of the IOPs proper. The 
%%% first part will calculate total FC PSDs, with cumulative,
%%% organic and inorganic IOPs following.

FCPSD_b = sum(tsd);                                                         % FC PSD total scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
FCPSD_bb = sum(tsda);                                                       % FC PSD total backscattering coefficient
FCPSD_bbr = (FCPSD_bb/FCPSD_b)*100;                                         % FC PSD total backscattering ratio
FCPSD_a_Qa = sum(sum(abs_ind));                                             % FC PSD total absorption coefficient

%%% This small section calculates the FC PSD component of the
%%% cumulative IOPs.

cmlPSD1_b = zeros(RV);                                                      % Preallocates cumulative scattering array
cmlPSD1_bb = zeros(RV);                                                     % Preallocates cumulative backscattering array
cmlPSD1_a_Qa = zeros(RV);                                                   % Preallocates cumulative absorption array

for c = 1:1:binNum                                                          % Iterates through size bins
    Beta_tot_cml = 0;                                                       % Primes the total particulate VSF
    for RI1 = 1:1:RIV                                                       % Iterates through RIs
        for R1 = 1:1:c                                                      % Iterates through radii up to bin c
            Beta_tot_cml = Beta_tot_cml + Beta_ind{RI1,R1};                 % Total VSF for this wavelength
        end
    end
    tsd = 2*pi*dtheta.*sintheta.*Beta_tot_cml';                             % Total particulate VSF integral
    cmlPSD1_b(c) = sum(tsd);                                                % FC PSD component of cumulative scattering
    tsda = tsd(index1:index2);                                              % Backscattering VSF integral
    cmlPSD1_bb(c) = sum(tsda);                                              % FC PSD component of cumulative backscattering
    cmlPSD1_a_Qa(c) = sum(sum(abs_ind(:,1:c),2),1);                         % FC PSD component of cumulative absorption
end

%%% The following section determines organic and inorganic IOPs.

%%% This next line finds the RI element closest to 1.1.
[~,RIind] = min(abs(RIs-1.1));

Beta_tot = 0;                                                               % Primes the total particulate VSF
orgPSD_a_Qa = 0;                                                            % Primes the total particulate organic absorption
for RI1 = 1:1:RIind-1                                                       % Iterates through organic RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        Beta_tot = Beta_tot + Beta_ind{RI1,R1};                             % Total VSF for this wavelength
        orgPSD_a_Qa = orgPSD_a_Qa + abs_ind(RI1,R1);                        % Total organic absorption for this wavelength
    end
end
tsd = 2*pi*dtheta.*sintheta.*Beta_tot';                                     % Total particulate VSF integral
orgPSD_b = sum(tsd);                                                        % FC PSD total organic scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
orgPSD_bb = sum(tsda);                                                      % FC PSD total organic backscattering coefficient
orgPSD_bbr = (orgPSD_bb/orgPSD_b)*100;                                      % FC PSD total organic backscattering ratio

Beta_tot = 0;                                                               % Primes the total particulate VSF
minPSD_a_Qa = 0;                                                            % Primes the total particulate inorganic absorption
for RI1 = RIind:1:RIV                                                       % Iterates through inorganic RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        Beta_tot = Beta_tot + Beta_ind{RI1,R1};                             % Total VSF for this wavelength
        minPSD_a_Qa = minPSD_a_Qa + abs_ind(RI1,R1);                        % Total inorganic absorption for this wavelength
    end
end
tsd = 2*pi*dtheta.*sintheta.*Beta_tot';                                     % Total particulate VSF integral
minPSD_b = sum(tsd);                                                        % FC PSD total inorganic scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
minPSD_bb = sum(tsda);                                                      % FC PSD total inorganic backscattering coefficient
minPSD_bbr = (minPSD_bb/ minPSD_b)*100;                                     % FC PSD total inorganic backscattering ratio

%%% 2.2 Extended PSD

%%% This section calculates the IOPs of the extended PSD by
%%% integrating the results of the FC PSD sections with IOP
%%% calculations for the sole extensions to the PSD. The first step 
%%% calculates total particulate IOPs, followed by total cumulative 
%%% IOPs and IOPs corresponding to the upper and lower extensions of 
%%% the PSD.

N = complPSD;                                                               % The complementary PSD array
N = N*1e6;                                                                  % This is a cm3 to m3 conversion

extRIs = D*0;                                                               % Primes RI array in the extensions

%%% The next loop assigns the average RIs to the extension as they 
%%% have been calculated earlier
for p = 1:1:RV                                                              % Iterates through radii
    if plotSIZEs(p) < 1e-6                                                  % First average assigned to particles smaller than 1 micron
        extRIs(p) = avgRI1;
    elseif plotSIZEs(p) > 1e-5                                              % Second average assigned to particles larger than 10 microns
        extRIs(p) = avgRI2;
    else                                                                    % Overall average assigned to particles between 1 and 10 microns
        extRIs(p) = avgRI;
    end
    %%% Imaginary RIs are assigned similarly as before.
    if extRIs(R1) < 1.1
        ni = niOrg_wvl;
    elseif extRIs(R1) >= 1.1
        ni = niMin_wvl;
    end
end

Beta_ind = ones(RV,thetalength);                                            % Primes total VSF
abs_ind = ones(RV,1);                                                       % Primes total absorption

for R1 = 1:1:RV                                                             % Iterates through radii
    m=extRIs(R1)+(1i*ni);                                                   % Complex refractive index
    w=2*pi/lambda;                                                          % Wavenumber
    x=w*r(R1);                                                              % Size Parameter
    [S1,S2,Qb,Qc,Qbb] = fastmie(x,m,[],theta_rad);                          % The core of the calculation is handled by the FASTMie script (by W. H. Slade, 2006).
    i1=abs(S1).^2;                                                          % Scattered intensity functions.
    i2=abs(S2).^2;                                                          % Scattered intensity functions.
    
    Beta_ind_core = ((1/w)^2)*0.5*((i1+i2));                                % Core VSF calculation for r(R1)
    abs_ind_core = (pi/4)*(D(R1).^2).*(Qc-Qb);                              % Core absorption calculation for r(R1)
    
    Beta_ind(R1,:)= Beta_ind_core.*N(R1);                                   % Total particulate VSF for r(R1)
    abs_ind(R1) = abs_ind_core.*N(R1);                                      % Total particulate absorption for r(R1)
end

Beta_tot = sum(Beta_ind,1);                                                 % Total VSF for the selected wavelength
tsd = 2*pi*dtheta.*sintheta.*Beta_tot;                                      % Total particulate VSF integral
complPSD_b = sum(tsd);                                                      % Total scattering coefficient in the extensions
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
complPSD_bb = sum(tsda);                                                    % Total backscattering coefficient in the extensions
complPSD_a_Qa = sum(abs_ind);                                               % Total absorption coefficient in the extensions

extPSD_b = FCPSD_b+complPSD_b;                                              % Total particulate scattering coefficient
extPSD_bb = FCPSD_bb+complPSD_bb;                                           % Total particulate backscattering coefficient
extPSD_bbr = (extPSD_bb/extPSD_b)*100;                                      % Total particulate scattering ratio
extPSD_a_Qa = FCPSD_a_Qa+complPSD_a_Qa;                                     % Total particulate absorption coefficient

%%% This small section calculates the extension component of the
%%% cumulative IOPs.

cmlPSD2_b = zeros(RV);                                                      % Preallocates cumulative scattering array
cmlPSD2_bb = zeros(RV);                                                     % Preallocates cumulative backscattering array
cmlPSD2_a_Qa = zeros(RV);                                                   % Preallocates cumulative absorption array

for c = 1:1:binNum                                                          % Iterates through size bins
    Beta_tot_cml = sum(Beta_ind(1:c,:),1);                                  % Total VSF for this wavelength
    tsd = 2*pi*dtheta.*sintheta.*Beta_tot_cml;                              % Total particulate VSF integral
    cmlPSD2_b(c) = sum(tsd);                                                % Extension component of cumulative scattering
    tsda = tsd(index1:index2);                                              % Backscattering VSF integral
    cmlPSD2_bb(c) = sum(tsda);                                              % Extension component of cumulative backscattering
    cmlPSD2_a_Qa(c) = sum(abs_ind(c));                                      % Extension component of cumulative absorption
end

cmlPSD_b = cmlPSD1_b+cmlPSD2_b;                                             % Total cumulative particulate scattering coefficient
cmlPSD_bb = cmlPSD1_bb+cmlPSD2_bb;                                          % Total cumulative particulate backscattering coefficient
cmlPSD_a_Qa = cmlPSD1_a_Qa+cmlPSD2_a_Qa;                                    % Total cumulative particulate absorption coefficient

%%% The following section determines IOPs in the upper and lower
%%% extensions

Beta_tot_upx = sum(Beta_ind(1:nz1-1,:),1);                                  % Total upper extension VSF for the selected wavelength
tsd = 2*pi*dtheta.*sintheta.*Beta_tot_upx;                                  % Total particulate VSF integral
upxPSD_b = sum(tsd);                                                        % Total upper extension particulate scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
upxPSD_bb = sum(tsda);                                                      % Total upper extension particulate backscattering coefficient
upxPSD_bbr = (upxPSD_bb/upxPSD_b)*100;                                      % Total upper extension particulate backscattering ratio
upxPSD_a_Qa = sum(abs_ind(1:nz1-1));                                        % Total upper extension particulate absorption coefficient

Beta_tot_lox = sum(Beta_ind(nz2+1:end,:),1);                                % Total lower extension VSF for the selected wavelength
tsd = 2*pi*dtheta.*sintheta.*Beta_tot_lox;                                  % Total particulate VSF integral
loxPSD_b = sum(tsd);                                                        % Total lower extension particulate scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
loxPSD_bb = sum(tsda);                                                      % Total lower extension particulate backscattering coefficient
loxPSD_bbr = (loxPSD_bb/loxPSD_b)*100;                                      % Total lower extension particulate backscattering ratio
loxPSD_a_Qa = sum(abs_ind(nz2+1:end));                                      % Total lower extension particulate absorption coefficient

%%% 2.3 Fluorescing PSD

%%% This section repeats the calculations made in section 2.1 but %%% this time only on the fluorescent fraction of the FC PSD.

N = FLPSDmat;                                                               % The size & RI concentration matrix
N = N*1e6;                                                                  % This is a cm3 to m3 conversion

Beta_ind = cell(RIV,RV);                                                    % Re-preallocates the cell matrix which will contain total fluorescent particulate VSFs
abs_ind = ones(RIV,RV);                                                     % Re-preallocates the matrix which will contain total fluorescent particulate absorptions

for RI1 = 1:1:RIV                                                           % Iterates through RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        Beta_ind{RI1,R1}= Beta_ind_wvl{RI1,R1}.*N(RI1,R1);                  % total particulate VSF for diameter r(R1) and refractive index RIs(RI1)
        abs_ind(RI1,R1) = abs_ind_wvl(RI1,R1).*N(RI1,R1);                   % total particulate absorption for diameter r(R1) and refractive index RIs(RI1)
    end
end

Beta_tot = 0;                                                               % Primes the total particulate VSF
FLPSD_a_Qa = 0;                                                             % Primes the total particulate fluorescent absorption
for RI1 = 1:1:RIind-1                                                       % Iterates through organic RIs
    for R1 = 1:1:RV                                                         % Iterates through radii
        Beta_tot = Beta_tot + Beta_ind{RI1,R1};                             % Total VSF for this wavelength
        FLPSD_a_Qa = FLPSD_a_Qa + abs_ind(RI1,R1);                          % Total fluorescent absorption for this wavelength
    end
end
tsd = 2*pi*dtheta.*sintheta.*Beta_tot';                                     % Total particulate VSF integral
FLPSD_b = sum(tsd);                                                         % FC PSD total fluorescent scattering coefficient
tsda = tsd(index1:index2);                                                  % Backscattering VSF integral
FLPSD_bb = sum(tsda);                                                       % FC PSD total fluorescent backscattering coefficient
FLPSD_bbr = (FLPSD_bb/FLPSD_b)*100;                                         % FC PSD total fluorescent backscattering ratio

