%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION A2 - FC model grid calculations, scaling and mapping.

%%% 2.1 Refractive indices & wavelength

%%% The refractive indices used in the model are defined here, as well
%%% as the wavelength of the laser light incident on the particles.
 
wnr = 1.333;                                                                % Water refractive index
 
lambda_air = 488e-9;                                                        % Laser light wavelength in air.
lambda = lambda_air/wnr;                                                    % Laser light wavelength in water.
 
%%% The refractive index of the standard is put first to simplify 
%%% scaling later on
 
nr = [1.595 (1.335:0.01:1.585) (1.605:0.01:1.725)]./wnr;                    % Relative real refractive indices array.
NV = numel(nr);                                                             % Gets the number of elements in the real refractive index array.
ni = 0;                                                                     % Imaginary refractive index - set at negligible absorption.
    
%%% 2.2 Size & refractive index isolines i.e. Look-up table

%%% This section of code calculates scattering value nodes, forming 
%%% curves of constant size and refractive index. These nodes 
%%% constitute a look-up table of scattering values which will relate 
%%% particle scattering to particle physical properties. Calculations
%%% are based on Mie theory and handled by an implementation of the 
%%% FASTMie code by W. H. Slade, http://www.scattport.org/index.php/
%%% light-scattering-software/mie-type-codes/list/264-fastmie.html
%%% The curves are then mapped to the dataset by using reference bead
%%% scattering averages.
    
%%% Two lines of code reassigning the two columns in the modelBeadSet
%%% array to two different variables.
    
beadsFSC = modelBeadSet(:,1);
beadsSSC = modelBeadSet(:,2);
 
%%% An array of particle sizes is now generated.

D = (logspace(-8,-4,300));                                                  % Array of log-spaced virtual particle diameters.
r = D/2;                                                                    % Array of log-spaced virtual particle radii.
RV = numel(D);                                                              % Gets the number of virtual particles.
 
sizes = [0.498;0.994;4.993;10.12;50.2;100]*1e-6;                            % The actual NIST diameters for the 6 bead groups used as standard, in metres.
 
[~,Bindex1] = min(abs(D-sizes(1)));                                         % Singles out the element in the diameters array corresponding to sizes(1).
[~,Bindex2] = min(abs(D-sizes(2)));                                         % Singles out the element in the diameters array corresponding to sizes(2).
[~,Bindex3] = min(abs(D-sizes(3)));                                         % Singles out the element in the diameters array corresponding to sizes(3).
[~,Bindex4] = min(abs(D-sizes(4)));                                         % Singles out the element in the diameters array corresponding to sizes(4).
[~,Bindex5] = min(abs(D-sizes(5)));                                         % Singles out the element in the diameters array corresponding to sizes(5).
[~,Bindex6] = min(abs(D-sizes(6)));                                         % Singles out the element in the diameters array corresponding to sizes(6).
    
If = ones(1,RV);                                                            % Preallocates the array which will contain the computed total forward scattered intensity of the virtual particles.
Is = ones(1,RV);                                                            % Preallocates the array which will contain the computed total side scattered intensity of the virtual particles.
 
FSC = ones(NV,RV);                                                          % Preallocates the array which will store modeled forward scatter for each refractive index/diameter pair.
SSC = ones(NV,RV);                                                          % Preallocates the array which will store modeled side scatter for each refractive index/diameter pair.
SIZEs = repmat(D,NV,1);                                                     % Preallocates the array which will store the virtual particle diameter for each refractive index/diameter pair.
RIs = repmat(nr',1,RV);                                                     % Preallocates the array which will store the virtual particle real RI for each refractive index/diameter pair.

%%% The following loop calculates the total forward and side scatter
%%% for each virtual particle.
    
for N1=1:1:NV                                                               % Iterates through refractive indices.
    for R1=1:1:RV                                                           % Iterates through radii.
        m = nr(N1)+(1i*ni);                                                 % Complex refractive index.
        k = 2*pi/lambda;                                                    % Wavenumber.
        x = k*r(R1);                                                        % Size parameter.
        [S1,S2,Qb,Qc,Qbb] = fastmie(x,m,[],theta_rad);                      % The core of the calculation is handled by the FASTMie script (by W. H. Slade, 2006).
        i1=abs(S1).^2;                                                      % Scattered intensity functions.
        i2=abs(S2).^2;                                                      % Scattered intensity functions.

        %%% The total forward and side scatter are now calculated using
        %%% an extremely stripped down version of the VSF integral normally 
        %%% used to calculate scattering - all missing factors and 
        %%% coefficients are taken care of by the subsequent mapping. 

        Iratio=i1+i2;                                                       % Simplified total scattered intensity (angular).

        tsd = 2*pi*dtheta'.*sintheta'.*Iratio;                              % The stripped down version of the VSF integral.
        tsdf = tsd(Iindex1:Iindex2);                                        % Preparing the forward scatter sensor angle integration.
        tsdfMOD = shapecorr1.*tsdf';                                        % Applying the sensor shape correction.  
        If(:,R1) = sum(tsdfMOD);                                            % Integrates and computes the total forward scattering by the virtual particle (N1,R1).

        tsds = tsd(Iindex3:Iindex4);                                        % Preparing the side scatter sensor angle integration.
        tsdsMOD = shapecorr2.*tsds';                                        % Applying the sensor shape correction.
        Is(:,R1) = sum(tsdsMOD);                                            % Integrates and computes the total forward scattering by the virtual particle (N1,R1).
    end
        
%%% 2.2.1 FC model grid/look-up table mapping

    %%% In the following section the script takes the refractive index
    %%% isolines it has computed and maps them to the standard beadset. To 
    %%% do so, it first calculates the X-axis and Y-axis displacement 
    %%% between the averages for the 0.5 um bead group and its virtual 
    %%% counterpart. It then translates the whole isoline so that the 
    %%% positions of the two match before rescaling the whole isoline by 
    %%% matching the X-axis and Y-axis distance between the 0.5 and 1 um 
    %%% bead groups and its virtual counterpart.

    if N1 == 1                                                              % The mapping is done with the RI nodes corresponding to the latex beads of the standard.

        DisplX = If(:,Bindex1)-beadsFSC(1);                                 % 0.5 um bead group X-axis displacement.
        DisplY = Is(:,Bindex1)-beadsSSC(1);                                 % 0.5 um bead group Y-axis displacement.

        FSCtemp = If-DisplX;                                                % Applies the X-axis displacement to the whole isoline.
        SSCtemp = Is-DisplY;                                                % Applies the Y-axis displacement to the whole isoline.

        ScaleX = (beadsFSC(2)-beadsFSC(1))/(FSCtemp(:,Bindex2)- ...     
            beadsFSC(1));                                                   % X-Axis scaling coefficient.
        ScaleY = (beadsSSC(2)-beadsSSC(1))/(SSCtemp(:,Bindex2)- ...     
            beadsSSC(1));                                                   % Y-Axis scaling coefficient.

        FSC(N1,:)=((FSCtemp-beadsFSC(1))*ScaleX)+beadsFSC(1);               % FWS values of the mapped isoline.
        SSC(N1,:)=((SSCtemp-beadsSSC(1))*ScaleY)+beadsSSC(1);               % SWS values of the mapped isoline.            
    end
 
%%% The process is then applied to all other refractive indices.
        
FSCtemp = If-DisplX;                                                        % Applies the X-axis displacement to the whole isoline.
SSCtemp = Is-DisplY;                                                        % Applies the Y-axis displacement to the whole isoline.
 
FSC(N1,:)=((FSCtemp-beadsFSC(1))*ScaleX)+beadsFSC(1);                       % FWS values of the mapped isoline.
SSC(N1,:)=((SSCtemp-beadsSSC(1))*ScaleY)+beadsSSC(1);                       % SWS values of the mapped isoline.
end
