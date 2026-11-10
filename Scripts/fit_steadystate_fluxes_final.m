% Fit steady state fluxes to GIT metabolites
function fit_steadystate_fluxes_final(inputFileName, outputFileName)
%addpath(genpath('.\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'
% Figures:
% 'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to cluster ions based on their GI profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file
%inputFileName = '.\Input\modeled_data_annDC.csv';
metabolite_data = readtable(inputFileName);

%outputFileName = '.\Output\model_results_annDC';
% file name to store solutions
filenameAllModel = [outputFileName '_all'];
filenameBestModel = [outputFileName '_best'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};
% separate data into a matrix of normalized intensities
% and get information on diet, tissue and mouse group
allColumns = metabolite_data.Properties.VariableNames;
% get columns that have git tissue suffix
dataColumns = zeros(size(allColumns));
for i=1:length(allColumns)
    dataColumns(i) = sum(arrayfun(@(x) contains(allColumns{i}, x), git_labels))>0;
end

% extract numeric data
meanMatrix = table2array(metabolite_data(:, dataColumns==1));
annColumns = allColumns(dataColumns==0);
dataColumns = allColumns(dataColumns==1);

% extract info about tissue, diet and mouse group
combinedTissues = cell(size(dataColumns));
combinedType = cell(size(dataColumns));
combinedDiet = cell(size(dataColumns));
for i=1:length(dataColumns)
    curcol = strsplit(dataColumns{i}, '_');
    combinedTissues{i} = curcol{3};
    combinedType{i} = curcol{2};
    combinedDiet{i} = curcol{1};
end
sampleType_unique = unique(combinedType);
sampleDiet_unique = unique(combinedDiet);
sampleTissue_unique = unique(combinedTissues);

% check if all tissues are present
if ~isempty(setdiff(git_labels, sampleTissue_unique))
    missing_tissues = setdiff(git_labels, sampleTissue_unique);
    for i=1:length(missing_tissues)
        sprintf('Warning: missing tissue: %s',missing_tissues{i});
    end
end

% check if GF group is present
if length(sampleType_unique) > 2
     sprintf('Warning: expect two groups, current group number is %d', length(sampleType_unique));
end
if sum(contains(sampleType_unique, {'GF'}))==0
   sprintf('Warning: missing GF group');
end
if ~isequal(sampleType_unique{end}, {'GF'})
    sampleType_unique = circshift(sampleType_unique,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selected_mets are metabolites detected along the GI tract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the mean normalized value across replicates

% define which metabolites to model
%%%% all metabolites
selected_mets = 1:size(meanMatrix,1);
met_gitfits = cell(length(selected_mets),1);
met_bestsols = cell(length(selected_mets),1);

for met_i = 1:length(selected_mets)
    cmpd_interest_idx = selected_mets(met_i);
    % calculate mean profiles            
    idx=1;
    kmeanMatrix_joint = zeros(4,6);
    kmeanMatrix_joint_names = cell(4,6);
    for diet_i = 1:length(sampleDiet_unique)
        for type_i = 1:length(sampleType_unique)
            selectDiet = sampleDiet_unique{diet_i};
            selectMouse = sampleType_unique{type_i};
            kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),dataColumns));
            kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
            kmeanMatrix_joint(idx,:) = kmeanMatrix;
            kmeanMatrix_joint_names(idx,:) = strcat(selectDiet,...
                                                    '_',...
                                                    selectMouse,...
                                                    '_',...
                                                    git_labels);
            idx = idx+1;
        end
    end
   
    % normalize by max intensity
    kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
    kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;

    if nnz(kmeanMatrix_joint)
        [gitfit] = fitGITmodel(kmeanMatrix_joint, kmeanMatrix_joint_names, 4, 1);
        % get best solution
        [bestsol] = select_gitfit_sol(gitfit);
        
        % save gitfit and best solution for the current metabolite
        met_gitfits{met_i} = gitfit;
        met_bestsols{met_i} = bestsol;
        
    end
    % display status every 100 metabolites
    if mod(met_i,100)==0
        fprintf('Calculated parameters for %d metabolites\n',met_i);
    end
end

% print solutions to file
% print best solutions to files
% create met_info object needed for the printing function
met_info = metabolite_data(selected_mets,annColumns);
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols, filenameBestModel);

% print solutions to files
print_allsols_to_files(met_info, met_gitfits, filenameAllModel);
