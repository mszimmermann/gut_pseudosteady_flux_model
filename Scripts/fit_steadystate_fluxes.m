% Fit steady state fluxes to GIT metabolites
function fit_steadystate_fluxes(inputFileName, outputFolder, outputFilePrefix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% inputFileName - name of input file with metabolite abundances and annotation information
% outputFolder - foldr name where the modelling results should be stored
% outputFilePrefix - prefix of the output file to be used for output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files: in the outputFolder wit houtputFilePrefix describing all and best
% solutions with differen optimization algorithms
% IP - interior point
% TRR = trust region reflective
% file types: 
% model_results*_all - all solutions
% model_results*_all_random_coefs - random coefficients for comparison
% model_results*_all_reciprocal - solution to the reverse problem
% model_results*_best# - best solution using each of the selection criteria #
% model_results*_best_reciprocal# - corresponding reciprocal solution using 
% each of the selection criteria#
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file
metabolite_data = readtable(inputFileName);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

% file name to store solutions
% check if outputFolder has prefix at the end
if ~(outputFolder(end)==filesep)
    outputFolder = [outputFolder filesep];
end

filenameAllModel = [outputFolder outputFilePrefix '_all'];
filenameBestModel = [outputFolder outputFilePrefix '_best'];

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
% expect GF group to be second of two
if ~isequal(upper(sampleType_unique(end)), {'GF'})
    sampleType_unique = circshift(sampleType_unique,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
