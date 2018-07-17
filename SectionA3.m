%%% Jacopo Agagliate
%%% University of Strathclyde
%%% Feb 2014 - Nov 2017
%%%
%%% ---- Changelog:
%%% 17 July 2018: Formatted and organized for GitHub storage.

%%% SECTION A3 - PMT sensitivity scaling of particle data and sample 
%%% dataset composition.

%%% 3.1 Rescaling and composition of measured particle scattering data

%%% Measured particle scattering data (main_data) is input here. In the 
%%% present iteration of the code it is meant to be input as a cell 
%%% array of four Nx3 matrices, each corresponding to a separate PMT 
%%% sensitivity run for the same sample. N is the number of particles 
%%% within each sensitivity run, and the three columns represent 
%%% forward scattering, side scattering and red fluorescence signal 
%%% respectively.

main_data_corr = main_data;                                                 % Preallocates cell array which will contain rescaled particle scattering values
                                
for q = 1:sensS_num
    if not(isempty(main_data{q}))
        main_data_corr{q}(:,1) = main_data{q}(:,1).* ...    			
            sensCorrCoeff(q,1);                                             % Sensitivity setting correction coefficient applied to FWS.
        main_data_corr{q}(:,2) = main_data{q}(:,2).* ...         
      		sensCorrCoeff(q,2);                                             % Sensitivity setting correction coefficient applied to SWS.
    end
end

%%% The following part of the code will select “strips” of particle 
%%% data to eliminate overlap between rescaled PMT sensitivity run data 
%%% and avoid inflating particle counts.
            
main_data_slice = cell(sensS_num,1);                                        % Preallocates cell array to contain the data “strips”
            
cutoff = [1e7,3500,700,200,1];                                              % SWS cutoff values for each sensitivity setting. User determined by direct examination of the data
cutoff_ind = cell(sensS_num,1);                                             % Prepares an array of cells to contain the arrays of indices of particles contained within the cut-off boundaries. 

for q = 1:sensS_num                                                         % This loop operates all the data cuts
    if not(isempty(main_data_corr{q}))                                      % Checks if the cell isn't empty (to prevent the code from stopping).
        cutoff_ind{q} = find(main_data_corr{q}(:,2) ...   
            < cutoff(q) & main_data_corr{q}(:,3) >= cutoff(q+1));           % Gets the indices of data points in main_data_corr{q} between cutoff(q) and cutoff(q+1).
        main_data_slice{q} = main_data_corr{q}(cutoff_ind{q},:);            % Assigns those data points to main_data_slice{q}.
    end
end
            
%%% The data “strips” are stitched back together into one single large 
%%% matrix. Fluorescence data is separated into a dedicated matrix.
            
main_data_total = [main_data_slice{1};main_data_slice{2}; ...                  
    main_data_slice{3};main_data_slice{4}];                                 % All cells from main_data_slice are compiled together into a single large array.

if not(isempty(main_data_total))
    main_data_total_FLR = main_data_total(:,3);                             % Extracts the Max FLR column from the main_data_total matrix.
    main_data_total = main_data_total(:,1:2);                               % Reduces the main_data_total matrix to the FWS and SWS columns only.
end 

%%% The number of particles left into each of the cells inside 
%%% main_data_slice are counted and kept track of. 
            
main_data_nums = [size(main_data_slice{1},1); ...  
    size(main_data_slice{2},1);size(main_data_slice{3},1);...
    size(main_data_slice{4},1)];
