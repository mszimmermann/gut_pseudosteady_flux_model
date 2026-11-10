% defined community data demo
inputFileName = '.\Input\modeled_data_annDCgut_demo.csv';
outputFolder = '.\Output\model_demo';
intermediateOutputFilePrefix = 'model_results_annDCgut_demo';
fit_steadystate_fluxes(inputFileName, outputFolder, ...
    intermediateOutputFilePrefix)


sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';
figurePrefix = '.\Output\fig_annDCgut_demo';
modelOutputFilePrefix = '.\Output\table_model_results_annDCgut_demo';

plot_and_analyse_modeling_results([outputFolder filesep intermediateOutputFilePrefix],...
                                    sel_crit1, sel_crit2,...
                                    figurePrefix, modelOutputFilePrefix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drug data demo
inputFileName = '.\Input\modeled_data_drugs.csv';
outputFolder = '.\Output\model_drugs';
intermediateOutputFilePrefix = 'model_results_drugs';
fit_steadystate_fluxes(inputFileName, outputFolder, ...
    intermediateOutputFilePrefix)


sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';
figurePrefix = '.\Output\fig_drugs';
modelOutputFilePrefix = '.\Output\table_model_results_drugs';

plot_and_analyse_modeling_results([outputFolder filesep intermediateOutputFilePrefix],...
                                    sel_crit1, sel_crit2,...
                                    figurePrefix, modelOutputFilePrefix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meier data demo
inputFileName = '.\Input\modeled_data_Meier.csv';
outputFolder = '.\Output\model_Meier';
intermediateOutputFilePrefix = 'model_results_Meier';
fit_steadystate_fluxes(inputFileName, outputFolder, ...
    intermediateOutputFilePrefix)


sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';
figurePrefix = '.\Output\fig_Meier';
modelOutputFilePrefix = '.\Output\table_model_results_Meier';

plot_and_analyse_modeling_results([outputFolder filesep intermediateOutputFilePrefix],...
                                    sel_crit1, sel_crit2,...
                                    figurePrefix, modelOutputFilePrefix)