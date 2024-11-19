% ########################################################################
%
% Various plots and analyses for EIB paper
%
% These are auxiliary commands to create figures and process data
% obtained with lnifmri_mrs_<dynamic|static>_analysis.m
% Inputs heavily rely on the availability of the particular outputs
% hardcoded in the script.
%
% stefano.tambalo@unitn.it - 20240502
% ########################################################################
clc
close all
curdir = pwd;
%addpath([curdir, filesep, 'gramm']);
%load([pwd, filesep, 'qualityReportMRS_workspace.mat']);
nameload = {'Rest', '0-Back', '1-Back', '2-Back'};

% Build a Quality x Subject Report (n=12 subj, dynamic MRS)
% ########################################################################
qualityReportMRSdir = 'MRS_struct_vox1';
qualityReportMRSlist = dir([pwd, filesep, qualityReportMRSdir, filesep, '*MRS_struct_vox1.csv']);
qualityReportMRSraw = table();
q=0;
for i = 1:size(qualityReportMRSlist, 1)
    dims(i) = size(readtable([qualityReportMRSlist(i).folder, filesep, qualityReportMRSlist(i).name]), 2);
    if dims(i) == 44
        q=q+1;
        qualityReportMRSraw(q, :) = readtable([qualityReportMRSlist(i).folder, filesep, qualityReportMRSlist(i).name]);
    end
end%for

% Data filtering
qualityReportMRSfilter = (qualityReportMRSraw.Date_of_analysis == '2023-12-05');
qualityReportMRStable = qualityReportMRSraw(qualityReportMRSfilter, :);

% Load subjects labels
resultspath = [pwd, filesep, 'RESULTS_dynamic_analysis',...
    filesep, '20230208_233743_batch1_GANNET'];
labelsfile = [resultspath, filesep, '20230208_233743_batch1_GANNET_LoadTimeSubjectLabels.csv'];
labelstable = readtable(labelsfile);

% Merge info
qualityReportMRStable = [qualityReportMRStable labelstable];
writetable(qualityReportMRStable,...
    [resultspath, filesep, '20230208_233743_batch1_GANNET_MRS_Quality_RAW.csv']);

% Compute summary
qualityReportMRSsummary = groupsummary(qualityReportMRStable(:, 5:47), "Subject", "all");
writetable(qualityReportMRSsummary,...
    [resultspath, filesep, '20230208_233743_batch1_GANNET_MRS_Quality_x_Subject.csv']);

save([pwd, filesep, 'qualityReportMRS_workspace.mat']);

%% Plot representative spectrum of a Subject to compare pre-post alignment
% Sub 9, Load 3 from static MRS results:
% 20231023_152702_MRS_034_19940705EETR_1_s25004_fMRS_GANNET
% ########################################################################

load([pwd, filesep, 'Gannet_MRS_struct_s9_l3.mat']);
edited(1, :) = real(MRS_struct.spec.vox1.GABAGlx.diff_noalign)+3000;
edited(2, :) = real(MRS_struct.spec.vox1.GABAGlx.diff);

editedV2 = [real(MRS_struct.spec.AllFramesFT)'; real(MRS_struct.spec.AllFramesFTrealign)'];
spectype = [repelem(1, size(MRS_struct.spec.AllFramesFT, 2), 1); repelem(2, size(MRS_struct.spec.AllFramesFT, 2), 1)];

% ########################################################################
% Only for retrospective reconstruction of plots.
% This line is not necessary for analysis run from scratch after 30/10/23
specfreq = readmatrix([pwd, filesep, 'specfreq.csv']);
% stefano.tambalo@unitn.it - 30/10/23
% ########################################################################

speclimits = find(specfreq>0 & specfreq<5);
gabalimits = find(specfreq>2.75 & specfreq<3.25);
glxlimits = find(specfreq>3.45 & specfreq<4.1);

%Plot individual spectra
individualspeclimits = find(specfreq>2.5 & specfreq<4.2);
%Color mapping for the polygons
cmap = [1 0.5 0.5; % red (Glx)
    0.5 0.5 1];% blue (GABA)
line = {'--', '-'};
linecol = {'noalign', 'align'};

% First figure
g = gramm('x', specfreq(individualspeclimits), 'y', edited(:, individualspeclimits), 'lightness', linecol);
g.geom_line();
g.geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');
g.geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');
g.no_legend();
g.set_names('x', 'Freq.(ppm)', 'y', '', 'color', '');
%Draw figure
g.set_title('');
g.set_color_options('lightness_range',[0 55],'chroma_range',[0 0]);
g.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {}, 'Ylim', [-2000 10000],...
                    'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [pwd, filesep, 'IndividualEditedSpectrum_alignment_before-after'],...
                'file_type', 'png', 'width', 18, 'unit', 'centimeters');
clear g
close(gcf)

% % Second figure
% graph = gramm('x', specfreq(individualspeclimits), 'y', edited(2, individualspeclimits), 'color', linecol);
% graph.geom_line();
% graph.set_color_options('map', [0 0 0; 0 0 0]);
% graph.geom_polygon('x',{[specfreq(glxlimits(1)), specfreq(glxlimits(end))];  [specfreq(gabalimits(1)), specfreq(gabalimits(end))]},'color', cmap, 'alpha', 0.3);
% graph.geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');
% graph.geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');
% graph.no_legend();
% graph.set_names('x', 'Freq.(ppm)', 'y', '', 'color', '');
% %Draw figure
% graph.set_title('Individual Edited Spectrum');
% graph.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {}, 'Ylim', [-2000 10000],...
%                     'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
% graph.draw();
% graph.export('file_name', [pwd, filesep, 'IndividualEditedSpectrum'],...
%                 'file_type', 'png', 'width', 18, 'unit', 'centimeters');
% clear graph
% close(gcf)

% Third figure
g = gramm('x', specfreq(individualspeclimits), 'y', editedV2(:, individualspeclimits), 'color', spectype);
g.stat_summary();
%graph.set_color_options('map', [0 0 0; 0 0 0]);
%g.geom_polygon('x',{[specfreq(glxlimits(1)), specfreq(glxlimits(end))];  [specfreq(gabalimits(1)), specfreq(gabalimits(end))]},'color', cmap, 'alpha', 0.3);
g.geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');
g.geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');
g.no_legend();
g.set_names('x', 'Freq.(ppm)', 'y', '', 'color', '');
%Draw figure
g.set_title('Individual Edited Spectrum');
g.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {},...%, 'Ylim', [-2000 10000],...
                    'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
%g.export('file_name', [pwd, filesep, 'IndividualEditedSpectrum], 'file_type', 'png');
%close(gcf)

%% Plot Water data (FWHM, F0shifts)
% From all MRS data (n=24, static MRS)
% to compare before/after QC exclusion criteria

staticQA_rejectfilterfile = [pwd, filesep, 'RESULTS_static_analysis', filesep,'20220628_090756_batch1_GANNET',filesep,'20220628_090756_batch1_GANNET_MRS_rejectfilter.csv'];
staticQA_rejectfilter = readtable(rejectfilterfile);
staticQA_labelsubj = repelem(1:24, 4)';
staticQA_labelload = repmat(1:4, 1, 24)';

% Path to individual static fMRS results
staticpath = [pwd, filesep, 'RESULTS_static_analysis', filesep, '20220628_090756_batch1_GANNET'];
% Get list of all datasets
staticlist = dir([staticpath, filesep, '*_fMRS_GANNET', filesep, 'Gannet_MRS_struct.mat']);
% Get QA data on water reference
for i=1:size(staticlist, 1)
   cursub = load([staticlist(i).folder, filesep, staticlist(i).name]);
   staticQA_waterFWHM(i, 1) = cursub.MRS_struct.out.vox1.water.FWHM;
   staticQA_GlxFWHM(i, 1) = cursub.MRS_struct.out.vox1.Glx.FWHM;
   staticQA_GABAFWHM(i, 1) = cursub.MRS_struct.out.vox1.GABA.FWHM;
   staticQA_avgDeltaF0(i, 1) = cursub.MRS_struct.out.AvgDeltaF0;
   staticQA_driftF0{i, 1} = 1:size(cursub.MRS_struct.spec.F0freq{1, 1}, 2);
   staticQA_driftF0{i, 2} = cursub.MRS_struct.spec.F0freq{1, 1};
   staticQA_GABAGlx{i, 1} = real(cursub.MRS_struct.spec.vox1.GABAGlx.diff_noalign(:, individualspeclimits));
   staticQA_GABAGlx{i, 2} = real(cursub.MRS_struct.spec.vox1.GABAGlx.diff(:, individualspeclimits));

end
%clear cursub i
% Compose datatable
staticQA_varnames = ["Subject", "Load", "DriftF0", "avgDeltaF0", "waterFWHM", "GABAGlx", "rejectFilter"];
staticQA_tableraw = table(staticQA_labelsubj, staticQA_labelload,...
                            staticQA_driftF0, staticQA_avgDeltaF0, staticQA_waterFWHM,...
                            staticQA_GABAGlx, staticQA_rejectfilter{:, 1},...
                            'VariableNames', staticQA_varnames);

%% Draw figures
% Boxplot of avgDeltaF0 x Load
xvar=nameload(categorical(staticQA_tableraw.Load));
yvar=staticQA_tableraw.avgDeltaF0;
cvar=staticQA_tableraw.Load;

staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);
xvar2=nameload(categorical(staticQA_passeddata.Load));
yvar2=staticQA_passeddata.avgDeltaF0;
cvar2=staticQA_passeddata.Load;

g = gramm('x', xvar, 'y', yvar, 'color', cvar);
g.stat_boxplot('width', 2);
g.no_legend();
g.set_names('x', 'Load', 'y', 'Avg. Delta F0 (ppm)', 'color', '');
g.set_title('');
g.set_color_options('chroma', 0, 'lightness', 95); %We make it light grey
g.set_order_options('color', 0, 'x', 0);

g.update('x', xvar2, 'y', yvar2, 'color', cvar2);
g.stat_boxplot('width', 1.5);   
g.set_color_options();
g.set_order_options('color', 0, 'x', 0);
g.no_legend();

g.axe_property('Ylim', [-0.15 0.15],'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_avgDeltaF0'],...
                'file_type', 'png', 'width', 18, 'unit', 'centimeters');
clear g
close(gcf)

%% Boxplot of waterFWHM x Load
xvar=nameload(categorical(staticQA_tableraw.Load));
yvar=staticQA_tableraw.waterFWHM;
cvar=staticQA_tableraw.Load;

staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);
xvar2=nameload(categorical(staticQA_passeddata.Load));
yvar2=staticQA_passeddata.waterFWHM;
cvar2=staticQA_passeddata.Load;

g = gramm('x', xvar, 'y', yvar, 'color', cvar);
g.stat_boxplot('width', 2);
g.no_legend();
g.set_names('x', 'Load', 'y', 'Water Peak (FWHM)', 'color', '');
g.set_title('');
g.set_color_options('chroma', 0, 'lightness', 95); %We make it light grey
g.set_order_options('color', 0, 'column', nameload);

g.update('x', xvar2, 'y', yvar2, 'color', cvar2);
g.stat_boxplot('width', 1.5);
g.set_color_options();
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();

g.axe_property('Ylim', [7 17],'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_waterFWHM'],...
            'file_type', 'png', 'width', 18, 'unit', 'centimeters');

close(gcf)
clear g

% Compute one-way np ANOVA to test for FWHM differences across loadings
% Test for normality
ksFWHM = kstest(yvar2);

% Distribution is not normal, then compute
% Kruskall-Wallis nonparametric ANOVA
[pKW_FWHM, tblKW_FWHM, statsKW_FWHM] =  kruskalwallis(yvar2, cvar2);
multcompare(statsKW_FWHM)


%% Line plot of FO drift over time x Load

staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);

g = gramm('x', staticQA_tableraw.DriftF0(:, 1), 'y', staticQA_tableraw.DriftF0(:, 2),...
                'color', staticQA_tableraw.Load);
g.facet_wrap(nameload(staticQA_tableraw.Load), 'ncols', 2);
g.stat_summary('type', 'sem');
g.set_color_options('chroma', 0, 'lightness', 75); %We make it light grey
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();
g.set_names('x', 'Averages', 'y', 'Water Freq. Drift (ppm)', 'color', '', 'column', '');
g.set_title('');

g.update('x', staticQA_passeddata.DriftF0(:, 1), 'y', staticQA_passeddata.DriftF0(:, 2),...
                 'color', staticQA_passeddata.Load);
g.facet_wrap(nameload(staticQA_passeddata.Load), 'ncols', 2);
g.stat_summary('type', 'sem');
g.geom_hline('yintercept', 4.68, 'style', 'k--');
g.set_color_options();
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();

g.axe_property('Ylim', [4.65 4.75],'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_waterDriftF0_summary'],...
            'file_type', 'png', 'width', 18, 'unit', 'centimeters');
close(gcf)
clear g

% Line plot of FO drift over time x Load - VERSION 2

staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);

g = gramm('x', staticQA_tableraw.DriftF0(:, 1), 'y', staticQA_tableraw.DriftF0(:, 2),...
                'color', staticQA_tableraw.Load);
g.facet_wrap(nameload(staticQA_tableraw.Load), 'ncols', 2);
g.stat_glm();
g.set_color_options('chroma', 0, 'lightness', 75); %We make it light grey
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();
g.set_names('x', 'Averages', 'y', 'Water Freq. Drift (ppm)', 'color', '', 'column', '');
g.set_title('');

g.update('x', staticQA_passeddata.DriftF0(:, 1), 'y', staticQA_passeddata.DriftF0(:, 2),...
                 'color', staticQA_passeddata.Load);
g.facet_wrap(nameload(staticQA_passeddata.Load), 'ncols', 2);
g.stat_glm();
g.geom_hline('yintercept', 4.68, 'style', 'k--');
g.set_color_options();
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();

g.axe_property('Ylim', [4.65 4.75],'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_waterDriftF0_glm'],...
            'file_type', 'png', 'width', 18, 'unit', 'centimeters');
close(gcf)
clear g

% Line plot of FO drift over time x Load - VERSION 3

staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);
yvar = staticQA_tableraw.DriftF0(:, 2);
yvar = yvar

g = gramm('x', staticQA_tableraw.DriftF0(:, 1), 'y', yvar,...
                'color', staticQA_tableraw.Load);
g.facet_wrap(nameload(staticQA_tableraw.Load), 'ncols', 2);
g.stat_glm();
g.set_color_options('chroma', 0, 'lightness', 75); %We make it light grey
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();
g.set_names('x', 'Averages', 'y', 'Water Freq. Drift (ppm)', 'color', '', 'column', '');
g.set_title('');

g.update('x', staticQA_passeddata.DriftF0(:, 1), 'y', staticQA_passeddata.DriftF0(:, 2),...
                 'color', staticQA_passeddata.Load);
g.facet_wrap(nameload(staticQA_passeddata.Load), 'ncols', 2);
g.stat_glm();
g.geom_hline('yintercept', 4.68, 'style', 'k--');
g.set_color_options();
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();

g.axe_property('Ylim', [4.65 4.75],'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_waterDriftF0_glm'],...
            'file_type', 'png', 'width', 18, 'unit', 'centimeters');
close(gcf)
clear g

%% Line plot of edited spectra
clear g
staticQA_passeddata = staticQA_tableraw(staticQA_tableraw.rejectFilter==0, :);

g = gramm('x', specfreq(individualspeclimits), 'y', staticQA_tableraw.GABAGlx(:, 1),...
                'color', staticQA_tableraw.Load);
g.facet_wrap(nameload(staticQA_tableraw.Load), 'ncols', 2);
g.stat_summary();
g.set_color_options('chroma', 0, 'lightness', 75); %We make it light grey
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();  
g.set_names('x', 'Frequency (ppm)', 'y', '', 'color', '', 'column', '');
g.set_title('');

g.update('x', specfreq(individualspeclimits), 'y', staticQA_passeddata.GABAGlx(:, 2),...
                 'color', staticQA_passeddata.Load);
g.facet_wrap(nameload(staticQA_passeddata.Load), 'ncols', 2);
g.stat_summary();
g.geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');
g.geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');
g.set_color_options();
g.set_order_options('color', 0, 'column', nameload);
g.no_legend();

g.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', '', 'Ylim', [-3000 6000],...
                'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
g.draw();
g.export('file_name', [staticpath, filesep, 'staticQA_GABAGlx'],...
            'file_type', 'png', 'width', 18, 'unit', 'centimeters');
close(gcf)
clear g
%% Save current state of matlab workspace
save([pwd, filesep, 'qualityReportMRS_workspace.mat']);

