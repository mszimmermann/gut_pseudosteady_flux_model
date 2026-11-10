% Plot and analyse modeling results
function plot_and_analyse_modeling_results(modelresultsFileName,...
                                            sel_crit1, sel_crit2,...
                                            figurePrefix, outputFilePrefix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'
% Matlab file with manually selected clustergrams
% 'cgo_clustergrams_of_model_coefficients.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
% 'table_hierarchical_clustering_groups.csv'
% Figures:
% 'fig4b_clustergram_2LIhos1LIbact_model_coefs_annotated_ions_in_GIT_Rmodel_corr_0_7'
% 'fig4b_barh_2LIhos1LIbact_model_cluster_sizes_percent.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modelresultsFileName = '.\Output\model_results_annDCgut_test';
filename = modelresultsFileName;
corrthreshold = 0.7;
% try combining solutions
% sel_crit1 = 'IP';
% sel_crit2 = 'LI_PCC_within_high_total';
% figurePrefix = '.\Output\fig_annDCgut_test';
% outputFilePrefix = '.\Output\table_model_results_annDCgut_test';

[met_info_combined, met_bestsol_combined] = ...
                combine_bestsols_from_file([filename '_best'], ...
                sel_crit1, sel_crit2,...
                corrthreshold);
% get model coefficients
x_met_smooth = met_bestsol_combined.x;
% get correlations calculated with reverse problem
x_data_corr = met_bestsol_combined.x_sel_CorrRev;
coefvalues = met_bestsol_combined.coefvalues;

% get all solutions from file
[met_info_read, met_gitfits_read] = read_allsols_from_files([filename '_all']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% calculate differentce in corr distrbutions
selected_mets = find((met_info_combined.SumGITclusters>0) &...
                     (met_info_combined.MetaboliteFilter==1));

%filename = [figurePrefix, '_histogram_MAXcorr_', sel_crit1, '_', sel_crit2];
%plot_gitfit_model_corr(met_gitfits_read(selected_mets), filename)

% plot correlation of restored and original data for the best solution
% calculate differentce in corr distrbutions
filename = [figurePrefix '_histogram_MAXcorr_best_combined_', sel_crit1, '_', sel_crit2];
gcfs = plot_bestfit_model_corr(filter_bestsols_by_index(met_bestsol_combined, selected_mets),...
                        met_gitfits_read(selected_mets), filename);
% close gcfs
for i=1:length(gcfs)
    close(gcfs(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save combined best solution to file
outputfilename = [outputFilePrefix, '_combined_',...
            sel_crit1, '_', sel_crit2];
% add information about the method and mode from the annotation table
met_info_method = cell(length(met_info_combined.MZ),1);
met_info_mode = zeros(length(met_info_combined.MZ),1);
met_info_combined.Method = met_info_method;
met_info_combined.Mode = met_info_mode;

print_combined_bestsol_to_files(met_info_combined, met_bestsol_combined,...
                                outputfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster only metabolites with corr>0.7
clustidx = 1:size(x_met_smooth,1); % for testing purposed plot all 

clustdata = x_met_smooth(clustidx,2:end)';
clustrows = coefvalues(2:end);
% filter by resiprocal corr
clustcorr = x_data_corr(clustidx);
clustdata = clustdata(:, clustcorr>=0.7);
clustidx = clustidx(clustcorr>=0.7);

for i=1:size(clustdata,2)
    clustdata(:,i) = clustdata(:,i)/max(abs(clustdata(:,i)));
end
clustnan = isnan(sum(clustdata,1));
clustdata(:,clustnan) = [];
clustidx(clustnan)=[];

colLabels = met_info_combined.CompoundID(clustidx);
colLabels = cellfun(@(x) strrep(x,'_','-'), colLabels, 'unif', 0);

clustdist ='cityblock';%'correlation';%'euclidean';%
cgo = clustergram(clustdata,...
            'RowLabels', clustrows,...
            'ColumnLabels', colLabels,...
            'ColumnPdist',clustdist,...
            'RowPdist', clustdist,...
            'DisplayRange', 1,...
            'colormap', redbluecmap);
% solution to turn on colormap programmatically from https://stackoverflow.com/questions/20648627/turn-on-colorbar-programmatically-in-clustergram
% combined with solution from https://de.mathworks.com/matlabcentral/answers/305274-reduce-font-size-of-column-labels-in-clustergram?#answer_533054
%cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
%cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle

cbButton = findall(0,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manually add colorbar to clustergram and print to figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = cgo.plot;
orient landscape
% add title to colobar
C = findall(gcf,'type','ColorBar');                         
C.Label.String = 'Relative model coefficient value';

print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
    [figurePrefix,...
    '_clustergram_model_combined_'...
    sel_crit1, '_', sel_crit2,'corr_0_7'])
close(gcf)

% close all 
set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');
for i=1:length(allhnds)
    close(allhnds(i))
end
