% ########################################################################
%
% Various plots and analyses for EIB paper
% 
% Comparison of different spectral alignment methods.
%
% stefano.tambalo@unitn.it - 20240502
% ########################################################################

clear all
clc
addpath('gramm');

% Hard-coded values from previous analyses.
resultspath = './RESULTS_dynamic_analysis';
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
keyword = 'batch1';
fitpackage = {'OSPREY', 'GANNET'};
fittingmethod = 2;

% Load Spectral Quality data
align_none = readtable('./RESULTS_dynamic_analysis/20231205_102512_batch1_GANNET_DynamicSpectralQuality.csv');
align_SR = readtable('./RESULTS_dynamic_analysis/20230208_233743_batch1_GANNET_DynamicSpectralQuality.csv');
align_rSR = readtable('./RESULTS_dynamic_analysis/20230307_111406_batch1_GANNET_DynamicSpectralQuality.csv');

% Define Spectral Alignment ID
align_IDnone = repelem(1, size(align_none, 1))';
align_IDSR = repelem(2, size(align_SR, 1))';
align_IDrSR = repelem(3, size(align_rSR, 1))';
align_ID = array2table([align_IDnone; align_IDSR; align_IDrSR], 'VariableNames', {'Align'});

% align_param = [align_none(:, [param 10]); align_SR(:, [param 10]); align_rSR(:, [param 10])];
align_param = [align_none; align_SR; align_rSR];
% Merge Data and ID to compose input dataset
aligndata = [align_param, align_ID];

% Select column to test
% Column 1/4: GABA/Glx FWHM
% Column 2/5: GABA/Glx SNR
% Column 3/6: GABA/Glx FitError
% Column 10: N-Back Load
varnames = align_none.Properties.VariableNames;
param = 5;
varname = varnames{param};

% Print summary of data
summary = groupsummary(aligndata, ["Align"], "all", varname);
writetable(summary,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_Summary_', varname, '_x_Align.csv']);

% Two-way anova
model = [varname, ' ~ Align'];
aov = anova(aligndata, model);
anovatable = stats(aov);
posthoc = multcompare(aov, "Align");

% Output ANOVA results (table and posthoc test)
writetable(aligndata,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_SpectralAlignment_anova_inputdata.csv']);
writetable(anovatable(1:end, :),...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_SpectralAlignment_', varname,'_anova_table.csv']);
writetable(posthoc,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_SpectralAlignment_', varname,'_anova_posthoc.csv']);

% Plot ANOVA results
loadnames = {'Rest', '0-Back', '1-Back', '2-Back'};
alignnames = ["None", "SpecReg", "rSpecReg"];

xvar = categorical(alignnames(aligndata.Align));
yvar = aligndata{:, param};
cvar = categorical(xvar);
% cvar = loadnames(aligndata.Load);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
% graph.set_point_options('base_size',3);
% graph.set_color_options();
% graph.geom_jitter('width', 0.4, 'dodge',0.3, 'alpha', 0.05);
% graph.no_legend();
%Add individual points
% graph.update('color', cvar);
graph.stat_boxplot();
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', varname, 'column', '', 'row', '', 'color', 'Load');
graph.set_title(['Comparison of Spectral Alignment - ', varname]);
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'Ylim', [5 30]);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
    '_SpectralAlignment_', varname,'_anova_boxplot'], 'file_type', 'png');
clear graph
% close(gcf)
