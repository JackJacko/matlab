%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION A1 - Angular ranges, sensor shape and PMT sensitivity scaling
%%% coefficients.

%%% 1.1 Angular ranges    

%%% The next few lines define the scattering angles array. The 
%%% resolution is higher at small angles.
 
theta1 = (0:0.01:1);                                                        % Angle range 1.
theta2 = (1.1:0.1:10);                                                      % Angle range 2.
theta3 = (11:1:180);                                                        % Angle range 3.

theta = [theta1 theta2 theta3];                                             % Scattering angles array.
thetalength = numel(theta);                                                 % Gets the number of scattering angles.
theta_rad = theta./180*pi;                                                  % Converts the scattering angles array into radians.
 
%%% The collection angle ranges of the forward and side scattering 
%%% sensors are defined here. These fine-tune the shape of the model 
%%% curves.
    
ang1 = 2;                                                                   % Forward scattering sensor MIN angle.
ang2 = 9.7;                                                                 % Forward scattering sensor MAX angle.
    
ang3 = 45;                                                                  % Side scattering sensor MIN angle.
ang4 = 135;                                                                 % Side scattering sensor MAX angle.
 
[~,Iindex1] = min(abs(theta-ang1));                                         % Singles out the element in the scattering angles array corresponding to ang1 for integration purposes.
[~,Iindex2] = min(abs(theta-ang2));                                         % Singles out the element in the scattering angles array corresponding to ang2 for integration purposes.
[~,Iindex3] = min(abs(theta-ang3));                                         % Singles out the element in the scattering angles array corresponding to ang3 for integration purposes.
[~,Iindex4] = min(abs(theta-ang4));                                         % Singles out the element in the scattering angles array corresponding to ang4 for integration purposes.
 
dtheta = ones(1,thetalength);                                               % Preallocates the array which contains the increments between angles in the scattering angles array. For integration purposes.

for S1 = 1:1:(thetalength-2)                                                % This small loop fills the angle increments array.
    dtheta(S1+1) = (theta_rad(S1+2)-theta_rad(S1))/2;
end
 
dtheta(1) = dtheta(2);                                                      % Defines the first element of dtheta.
dtheta(thetalength) = dtheta(thetalength-1);                                % Defines the last element of dtheta.
 
sintheta = sin(theta_rad);                                                  % The sines of the scattering angles.

%%% 1.2 Sensor shape weighting functions

%%% This section defines the weighting functions which account for the 
%%% collecting area of the sensors.
 
% Forward scattering weighting function
shapecorr1=((.5*pi)-asin(sin(theta_rad(Iindex1))./ ... 
    sin(theta_rad(Iindex1:Iindex2))))/(.5*pi);
 
% Side scattering weighting function element 1 
shapecorr2space=sin(linspace(0,pi,numel(theta_rad(Iindex3:Iindex4))));

% Side scattering weighting function element 2 
Kprime = sin(theta_rad(Iindex4)-(.5*pi)); 

% Side scattering weighting function element 3                              
shapecorr2core = Kprime*shapecorr2space;

% Side scattering weighting function    
shapecorr2=asin(shapecorr2core)/(.5*pi);

%%% 1.3 Reference polymer beads and PMT sensitivity settings

%%% This section will define reference polymer bead scattering averages 
%%% and PMT sensitivity settings which will be used both to scale and 
%%% centre the scattering values of the FC model grid/look-up table and 
%%% to combine data from multiple sensitivity runs into a single sample 
%%% dataset. PMT setting 60 was the one used as reference.
    
sensSetting = [50;60;70;80];                                                % Sensitivity settings array.
sensS_num = numel(sensSetting);                                             % Gets the number of elements in the sensitivity settings array.

%%% The following line will require input from the user. It is meant to 
%%% collect average scattering values for 1 micron reference beads into 
%%% a 4x2 matrix, with each row corresponding to a PMT setting and the 
%%% two columns corresponding to forward and side scattering 
%%% respectively.     

set1um = [avgs_1um{1};avgs_1um{2};avgs_1um{3};avgs_1um{4}];                 % The set of scattering averages corresponding to 1 micron beads

%%% The following line will also require input from the user. It is 
%%% meant to collect average scattering values for 5 and 10 micron 
%%% reference beads into a 2x2 matrix, with each row corresponding to a 
%%% bead diameter and the two columns corresponding to forward and side 
%%% scattering respectively. 
    
setMed60 = [avgs_5um;avgs_10um];                                            % The set of scattering averages corresponding to 5 and 10 micron beads

%%% The following line will require a final input from the user. It is 
%%% meant to collect average scattering values for 50 and 100 micron 
%%% reference beads into a 2x2 matrix, with each row corresponding to a 
%%% bead diameter and the two columns corresponding to forward and side 
%%% scattering respectively. PMT setting 50 was the only one for which 
%%% side scattering of 50 and 100 micron reference beads didn’t 
%%% saturate.

setLarge50 = [avgs_50um;avgs_100um];                                        % The set of scattering averages corresponding to 50 and 100 micron beads                                        

%%% This next section uses the 1 micron set of bead scattering averages 
%%% to work out scaling coefficients necessary to compose data 
%%% corresponding to each sensitivity setting into one single dataset.
%%% The type of the following fit models was determined after direct
%%% examination of the data.
    
f1 = coeffvalues(fit(sensSetting,set1um(:,1),'poly1'));                     % Calculates the linear fit of the FWS averages of the 1 um set.
f2 = coeffvalues(fit(sensSetting,set1um(:,2),'power1'));                    % Calculates the power model fit of the SWS averages of the 1 um set.
    
fws50coeff = set1um(2,1)/(f1(1)*sensSetting(1)+f1(2));                      % FWS scaling coefficient for sensitivity setting 60.
sws50coeff = set1um(2,2)/(f2(1)*(sensSetting(1)^f2(2)));                    % SWS scaling coefficient for sensitivity setting 60.
fws70coeff = set1um(2,1)/(f1(1)*sensSetting(3)+f1(2));                      % FWS scaling coefficient for sensitivity setting 70.
sws70coeff = set1um(2,2)/(f2(1)*(sensSetting(3)^f2(2)));                    % SWS scaling coefficient for sensitivity setting 70.
fws80coeff = set1um(2,1)/(f1(1)*sensSetting(4)+f1(2));                      % FWS scaling coefficient for sensitivity setting 80.
sws80coeff = set1um(2,2)/(f2(1)*(sensSetting(4)^f2(2)));                    % SWS scaling coefficient for sensitivity setting 80.
    
sensCorrCoeff = [[fws50coeff;1;fws70coeff;fws80coeff], ...              
    [sws50coeff;1;sws70coeff;sws80coeff]];                                  % Scaling coefficient matrix
 
%%% The next two lines scale the scattering averages corresponding to 
%%% 50 and 100 micron beads to PMT setting 60
    
fwsLarge60 = setLarge50(:,1)*sensCorrCoeff(1,1);
swsLarge60 = setLarge50(:,2)*sensCorrCoeff(1,2);
    
%%% The next line compiles the set of bead scattering averages 
%%% corresponding to reference PMT setting 60. It also requires input 
%%% of the 0.5 micron bead scattering averages at PMT setting 60, which 
%%% was the lowest setting for which 0.5 micron beads were detectable.

modelBeadSet = [[avgs_05um(1);avgs_1um{2}(1);setMed60(:,1); ... 
    fwsLarge60],[avgs_05um(2);avgs_1um{2}(2);setMed60(:,2);swsLarge60]];
