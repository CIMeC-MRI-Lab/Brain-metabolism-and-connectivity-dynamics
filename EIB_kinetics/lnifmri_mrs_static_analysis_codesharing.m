% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20220504
% # REV. 20220630
% # REV. 20240503
% ########################################################################

function lnifmri_mrs_static_analysis(fittingmethod, keyword)

%clear variables
%close all
clc

compThread = maxNumCompThreads(8);

curdir = pwd;
ospreypath=[curdir, filesep, 'osprey'];
gannetpath=[curdir, filesep, 'Gannet'];
exportfigpath=[gannetpath, filesep, 'export_fig'];
addpath([curdir, filesep, 'spm12']);
addpath([curdir, filesep, 'gramm']);

% Set output folder name
resultsdirname = 'RESULTS_static_analysis';
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if ~isfolder([curdir, filesep, resultsdirname])
    createresultsdir = ['mkdir ', curdir, filesep, resultsdirname];
    system(createresultsdir);
end
resultspath = [curdir, filesep, resultsdirname];

% ########################################################################
% ##### EDIT SECTION #####################################################
% ########################################################################
% Select fitting method:
% 1 = OSPREY
% 2 = GANNET
%
% fittingmethod = 1;
%
% If GANNET, then specify the fitting model:
model = 'unused';
% model = '1lorentz';
%
% If OSPREY, then optionally add a custom basis set:
% basissetfolder = [ospreypath, filesep, 'MRS_basis_set', filesep, 'all_metabolites', filesep, 'basis'];
% BASIS = fit_makeBasis(basissetfolder, 1, 'MEGA', 'GABA');
% basissetfolder = [ospreypath, filesep, 'MRS_basis_set', filesep, 'all_metabolites', filesep, 'GABAMM'];
% BASIS = load([basissetfolder, filesep, 'BASIS_Siemens_MEGA_PRESS_GABA68_wMM.mat']);
%
% Specify a substring of relevant data folders to process
% e.g.: 'sub' to process all folders whose name starts with 'sub'
%
% keyword = 'batch2';
%
%
% ########################################################################
% ##### END OF EDIT SECTION ##############################################
% ########################################################################


% Set fitting model
fitpackage = {'OSPREY', 'GANNET'};

% Clean output directory
outdir = [timestamp, '_', keyword, '_', fitpackage{fittingmethod}];
mkoutdir = ['mkdir ', resultspath, filesep, outdir];
system(mkoutdir);
% Update resultspath
resultspath = [resultspath, filesep, outdir];

%Start analysis
listsub = dir([curdir, filesep, keyword, '*']);
listdataset = dir([curdir, filesep, keyword, '*', filesep, 'MRI_MRS', filesep, '*_fMRS']);
metabdiff = nan(size(listdataset, 1), 32);
metaboff = nan(size(listdataset, 1), 32);
metabqual = nan(size(listdataset, 1), 9);
specqual = zeros(size(listdataset, 1), 32769);

namesmetab = cell(1, 31);
namesrow = cell(size(listdataset, 1), 1);
for k =1:size(namesrow, 1)
    namesrow(k) = {[num2str(k, '%04d'), '_NaN']};
end

%initialize row counter to store results
r=0;
for s = 1:size(listsub, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag to enable coreg and segment only on the first run of a new
    % subject - CURRENTLY UNUSED
    same = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    listruns = dir([listsub(s).folder, filesep, listsub(s).name, filesep, 'MRI_MRS', filesep, '*_fMRS']);
    % No need to edit beyond this line
    % ########################################################################

    for n = 1:size(listruns, 1)
        % Save current state of matlab workspace
        save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);
        % ########################################################################
        % Set paths
        datapath = [pwd, filesep, listsub(s).name, filesep, 'MRI_MRS'];
        % define path for GLOBAL DATASET FOLDER
        rootpath = listruns(n).name;
        % get available folders
        availdata = dir([datapath, filesep, rootpath, filesep, '*']);
        dirfilter = contains({availdata.name}, '.');
        % filter '.' and '..' directories
        availdata(dirfilter) = [];
        % define path for EDITED MRS DATA
        mrsfilter = contains({availdata.name}, 'EDIT');
        mrspath = availdata(mrsfilter).name;
        % define path for WATER REFERENCE DATA
        reffilter = contains({availdata.name}, 'WREF');
        refpath = availdata(reffilter).name;
        % define path for ANATOMICAL REFERENCE DATA
        anatpath = 'anat';

        % set number of averages
        onoffspectra = dir([datapath, filesep, rootpath, filesep, mrspath, filesep, '*.IMA']);
        number = length(onoffspectra); %this number must be even

        % define basename of metabolite data files
        onoffnames = {onoffspectra.name};
        firstmetabfile = onoffnames{1};
        tokens = strsplit(firstmetabfile, '.');
        tokenID = size(tokens, 2)-10;
        imaname = [strjoin(tokens(1:tokenID), '.'), '.'];
        imaID = tokens{1};

        % define filename of first water-reference DICOM file
        wrefspectra = dir([datapath, filesep, rootpath, filesep, refpath, filesep, '*.IMA']);
        wrefnames = {wrefspectra.name};
        refname = wrefnames{1};

        % define the path for your data folder
        path_metab = [datapath, filesep, rootpath, filesep, mrspath];%, filesep, imaname];
        % define the path to the first of your water reference files
        path_ref = [datapath, filesep, rootpath, filesep, refpath];
        % define the path to the anatomical image
        path_t1 = [datapath, filesep, anatpath];

        % create output folder
        unique_name = strrep(rootpath, filesep, '_');
        outputfoldername = [timestamp, '_MRS_', num2str(r, '%03d'), '_', imaID, '_', unique_name, '_', fitpackage{fittingmethod}];
        outputfolder = [resultspath, filesep, outputfoldername];
        createoutputfolder = ['mkdir ', outputfolder];
        system(createoutputfolder);

        % copy or convert anatomical image
        anatlist = dir([path_t1, filesep, '*.nii']);
        if ~isempty(anatlist) && contains(anatlist.name, 'anat')
            disp('Copying existing structural image to analysis folder');
            copyanat = ['cp ', path_t1, filesep, 'anat.nii ', outputfolder];
            system(copyanat);
        elseif ~isempty(anatlist) && contains(anatlist.name, 'mprage')
            disp('Copying existing structural image (raw MEmprage) to analysis folder');
            copyanat = ['cp ', path_t1, filesep, anatlist.name, ' ', outputfolder, filesep, 'anat.nii'];
            system(copyanat);
        else
            disp('Converting and copying existing structural image to analysis folder');
            nifticonversion = ['dcm2niix -b y -o ' outputfolder, ' -f anat ', path_t1];
            system(nifticonversion);
        end%if

        switch fittingmethod
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 1 % Run OSPREY
                addpath(genpath(ospreypath));
                jobname = 'lnifmri_mrs_osprey_MM_job';

                MRS_struct = lnifmri_mrs_static_osprey(jobname, path_metab, path_ref, number, same, outputfolder);

                save([outputfolder, filesep, 'MRS_struct.mat']);
                movejobmat = ['mv ', jobname, '.mat ', outputfolder];
                system(movejobmat);
                % increase row counter
                r = r+1;
                % append load index (rest, 0Back, 1Back, 2Back)
                metabdiff(r, :) = [table2array(MRS_struct.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 2}), n];
                metaboff(r, :) = [table2array(MRS_struct.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 1}), NaN, n];
                metabqual(r, 1:numel(MRS_struct.QM.tables{1, :})) = MRS_struct.QM.tables{1, :};
                metabqual(r, end) = n;

                namesrow{r} = [num2str(n, '%03d'), '_', num2str(r, '%03d'), '_', listdataset(r).name];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % flag to enable coreg and segment only on the first run of a new
                % subject - CURRENTLY UNUSED
                same = 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                namesmetab = [MRS_struct.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 2}.Properties.VariableNames, 'Load'];
                namesqual = [MRS_struct.QM.tables.Properties.VariableNames, 'x001', 'Load'];

                specsize = numel(MRS_struct.processed.metab{1,1}.ppm);
                specqual(r, 1:specsize) = real(MRS_struct.processed.metab{1,1}.specs(:,3));
                specqual(r, end) = n;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 2 % Run GANNET
                addpath(genpath(gannetpath));
                MRS_struct = lnifmri_mrs_static_gannet(model, path_metab, path_ref, refname, number, same, gannetpath, outputfolder);
                r = r+1;
                metabdiff(r, 1) = MRS_struct.out.vox1.GABA.ConcIU_TissCorr;
                metabdiff(r, 2) = MRS_struct.out.vox1.GABA.ConcIU_AlphaTissCorr;
                metabdiff(r, 3) = MRS_struct.out.vox1.Glx.ConcIU_TissCorr;
                metabdiff(r, 4) = MRS_struct.out.vox1.Glx.ConcIU_AlphaTissCorr;
                metabdiff(r, end) = n;

                metaboff(r, 1) = n;
                metaboff(r, 2) = n;
                metaboff(r, 3) = n;
                metaboff(r, 4) = n;
                metaboff(r, 5) = n;
                metaboff(r, end) = n;

                namesrow{r} = [num2str(n, '%03d'), '_', num2str(r, '%03d'), '_', listdataset(r).name];
                namesmetab = cell(1, 31);
                namesmetab(1:2) = {'GABA_ConcIU_TissCorr', 'GABA_ConcIU_AlphaTissCorr'};
                namesmetab(3:4) = {'Glx_ConcIU_TissCorr', 'Glx_ConcIU_AlphaTissCorr'};
                emptyplaces = numel(namesmetab(5:end));
                namesmetab(5:end) = strsplit(deblank(sprintf('x%03d ', 1:emptyplaces)));
                namesmetab = [namesmetab, 'Load'];

                metabqual(r, 1) = MRS_struct.out.vox1.GABA.FWHM;
                metabqual(r, 2) = MRS_struct.out.vox1.GABA.SNR;
                metabqual(r, 3) = MRS_struct.out.vox1.GABA.FitError;
                metabqual(r, 4) = MRS_struct.out.vox1.Glx.FWHM;
                metabqual(r, 5) = MRS_struct.out.vox1.Glx.SNR;
                metabqual(r, 6) = MRS_struct.out.vox1.Glx.FitError;
                metabqual(r, 7) = MRS_struct.out.vox1.tissue.fGM;
                metabqual(r, 8) = MRS_struct.out.vox1.ResidWater.SuppressionFactor;
                metabqual(r, 9) = n;
                namesqual = cell(1, 8);
                namesqual(1:8) = {'GABA_FWHM', 'GABA_SNR', 'GABA_FitError', 'Glx_FWHM', 'Glx_SNR', 'Glx_FitError', 'fGM', 'WaterSuppFactor'};
                %namesqual = strsplit(deblank(sprintf('x%03d ', 1:emptyplaces)));
                namesqual = [namesqual, 'Load'];

                specsize = numel(MRS_struct.spec.vox1.GABAGlx.diff_scaled);
                specqual(r, 1:specsize) = real(MRS_struct.spec.vox1.GABAGlx.diff_scaled);
                specqual(r, end) = n;
        end%switch

        save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.mat'], 'metabdiff');
        save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metaboff.mat'], 'metaboff');
        save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabqual.mat'], 'metabqual');
        save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specqual.mat'], 'specqual');

        difftable = array2table(metabdiff, 'VariableNames', namesmetab, 'RowNames', namesrow);
        writetable(difftable,[resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.csv'],'WriteRowNames',true)

        offtable = array2table(metaboff, 'VariableNames', namesmetab, 'RowNames', namesrow);
        writetable(offtable,[resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metaboff.csv'],'WriteRowNames',true)

        qualtable = array2table(metabqual, 'VariableNames', namesqual, 'RowNames', namesrow);
        writetable(qualtable,[resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabqual.csv'],'WriteRowNames',true)

    end%listnback
end%listsub

% Save only once the frequency array (used later for plotting)
switch fittingmethod
    case 1
        specfreq = MRS_struct.processed.metab{1,1}.ppm;
    case 2
        specfreq = MRS_struct.spec.freq;
end
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specfreq.mat'], 'specfreq');

% Save current state of matlab workspace
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Glx/GABA ratio

%Load metrics for OSPREY(1) or GANNET(2)
switch fittingmethod
    case 1
        fGlx = {'Glx'};
        fGABA = {'GABAplus'};
        namevars = {'GlxGABAplus', 'Load'};
    case 2
        fGlx = {'Glx_ConcIU_TissCorr', 'Glx_ConcIU_AlphaTissCorr'};
        fGABA = {'GABA_ConcIU_TissCorr', 'GABA_ConcIU_AlphaTissCorr'};
        namevars = {'GlxGABA_ConcIU_TissCorr', 'GlxGABA_ConcIU_AlphaTissCorr', 'Load'};
end%switch

diffcsv = dir([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.csv']);
glxgabadata = readtable([diffcsv.folder, filesep, diffcsv.name], 'Delimiter', 'comma');
glx = glxgabadata{:, fGlx};
gaba = glxgabadata{:, fGABA};
glxgabaratio = glx./gaba;
glxgabaratio = [glxgabaratio, glxgabadata{:, {'Load'}}];
glxgabatable = array2table(glxgabaratio, 'VariableNames', namevars, 'RowNames', namesrow);
writetable(glxgabatable,[resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabGlxGABAratio.csv'],'WriteRowNames',true)

% Save current state of matlab workspace
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quality-based rejection filter

%Load metrics for OSPREY(1) or GANNET(2)
switch fittingmethod
    case 1
        QCmetrics = {'NAA_SNR', 'NAA_FWHM', 'water_FWHM'};
    case 2
        QCmetrics = {'GABA_SNR', 'GABA_FitError', 'Glx_SNR', 'Glx_FitError'};
end%switch

%Load quality parameters
listcsv = dir([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabqual.csv']);
dataqual = readtable([listcsv.folder, filesep, listcsv.name], 'Delimiter', 'comma');
loadindex = unique(dataqual.Load);

%Initialize evaluation structure
qualitytable = 2.*ones([size(dataqual, 1)/numel(loadindex) numel(QCmetrics) numel(loadindex)]);

%Z-score of quality metrics
for i = 1:numel(loadindex)
    filter = dataqual.Load == loadindex(i);
    qualitytable(:,:, i) = abs(zscore(dataqual{filter, QCmetrics}));
end

%Exclude data with z-score > threshold
zth = 1;
qualitytable(qualitytable < zth) = 0;
qualitytable(qualitytable >= zth) = 1;
% If the sum is greater than one, it means a z-score>2 in more than one
% dimension (quality parameter or run)
exclusionmatrix = sum(qualitytable, [2 3]);
%Build a subject-rejection filter
rejectfilter = repelem(exclusionmatrix>1, numel(loadindex));
%Store rejected subjects in a csv file
rejectedlist = dataqual(rejectfilter, 'Row');
writetable(rejectedlist, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_rejectedList.csv'],'WriteRowNames',true);
writematrix(rejectfilter, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_rejectfilter.csv']);

% Save current state of matlab workspace
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANOVA Stats and Graphical Output
listcsv = dir([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.csv']);

figuresize = get(0,'ScreenSize');
%figuresize(3:4) = figuresize(3:4)./2;
namerun = {'Rest', '0-Back', '1-Back', '2-Back'};
figuresize(4) = figuresize(4)/2;
figuresize(3) = (figuresize(3)/8)*numel(namerun);

for p = 1:size(listcsv, 1)
    %find(specfreq>2.75 & specfreq<3.25);
    % load data file
    datatoplot = readtable([resultspath, filesep, listcsv(p).name], 'Delimiter', 'comma');
    datatoplot = [datatoplot(:, "GABA_ConcIU_AlphaTissCorr"),...
        datatoplot(:, "Glx_ConcIU_AlphaTissCorr"),...
        glxgabatable(:, "GlxGABA_ConcIU_AlphaTissCorr"),...
        datatoplot(:, "Load")];
    % exclude subjects
    datatoplot(rejectfilter, :) = [];
    datatoplot_outlierremoval = filloutliers(datatoplot(:, 2:end-1), 'linear');
    writetable(datatoplot_outlierremoval, [resultspath, filesep, [strtok(listcsv(p).name, '.'), '_nooutliers.txt']]);
    filename = strtok(listcsv(p).name, '.');
    varnames = datatoplot.Properties.VariableNames;

    % skip rowname, skip Load
    for r = 1:size(varnames, 2)-1
        vartoplot = varnames{r};
        % Plot only sensible data
        % if ~contains(vartoplot, 'x0')
        if contains(vartoplot, 'AlphaTissCorr')
            % Filter out NaN values
            nanfilter = ~isnan(datatoplot.(vartoplot));
            filteredX = datatoplot.Load(nanfilter);
            filteredY = datatoplot.(vartoplot)(nanfilter);
            % Scale all values considering Rest run as subject-specific baseline
            baseline = repelem(filteredY(filteredX==1), 4);
            filteredY = (filteredY-baseline)./baseline;
            % This is to improve readability of the code
            % Replace numeric Load code with readable literal labels
            groupvar = namerun(filteredX);
            % Derive number of subjects
            subnum = nnz(rejectfilter)/numel(namerun);
            % Define grouping factor
            factors = repmat(string(namerun), subnum, 1);

            % Stats - One-way anova on AlphaTissueCorrected concentration
            % Run ANOVA 1-way
            aov = anova(factors(:), filteredY(:), FactorNames="Load");
            multcompare(aov)
            anovatable = stats(aov);
            posthoc = multcompare(aov);
            writetable(array2table([factors(:), filteredY(:)]),...
                [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
                '_ANOVA_', vartoplot,'_inputdata.csv']);
            writetable(anovatable(1:end, :),...
                [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
                '_ANOVA_', vartoplot,'_table.csv']);
            writetable(posthoc,...
                [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
                '_ANOVA_', vartoplot,'_posthoc.csv']);


            xvar = categorical(factors(:));%reshape(aov.Factors{:, 'Load'}, [], 4);
            yvar = filteredY(:);%reshape(aov.Y, [], 4);
            cvar = categorical(xvar);

            graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
            graph.stat_boxplot();
            graph.set_order_options('color', 0, 'x', 0);
            graph.set_names('x', '', 'y', vartoplot, 'column', '', 'row', '', 'color', 'Load');
            graph.set_title([vartoplot, ' 1-way ANOVA']);
            graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
            graph.no_legend();
            %Add individual points
            graph.update('color', cvar);
            graph.set_point_options('base_size',6);
            graph.set_color_options();
            graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
            graph.no_legend();
            graph.draw();
            graph.export('file_name', [resultspath, filesep, outdir, '_ANOVA_', vartoplot, '_boxplot'], 'file_type', 'png');
            clear graph
            close(gcf)

        else
        end
    end%for r
end%for p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Spectral Quality for GABA and Glx

clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/1;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

load([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specqual.mat']);
load([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specfreq.mat']);
specqual(rejectfilter, :) = [];
speclimits = find(specfreq>0 & specfreq<5);
gabalimits = find(specfreq>2.75 & specfreq<3.25);
glxlimits = find(specfreq>3.45 & specfreq<4.1);

clear graph

%Plot GABA
graph(1, 1) = gramm('x', specfreq(speclimits), 'y', specqual(:, speclimits), 'color', groupvar);
graph(1, 1).geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');

graph(1, 1).stat_summary();
graph(1, 1).facet_grid([], groupvar);
graph(1, 1).set_order_options('color', 0, 'column', 0);
graph(1, 1).no_legend();
graph(1, 1).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

graph(2, 1) = gramm('x', specfreq(gabalimits), 'y', specqual(:, gabalimits), 'color', groupvar);
graph(2, 1).stat_summary();
graph(2, 1).facet_grid([], groupvar);
graph(2, 1).set_order_options('color', 0, 'column', 0);
graph(2, 1).no_legend();
graph(2, 1).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

%Plot Glx
graph(1, 2) = gramm('x', specfreq(speclimits), 'y', specqual(:, speclimits), 'color', groupvar);
graph(1, 2).geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');

graph(1, 2).stat_summary();
graph(1, 2).facet_grid([], groupvar);
graph(1, 2).set_order_options('color', 0, 'column', 0);
graph(1, 2).no_legend();
graph(1, 2).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

graph(2, 2) = gramm('x', specfreq(glxlimits), 'y', specqual(:, glxlimits), 'color', groupvar);
graph(2, 2).stat_summary();
graph(2, 2).facet_grid([], groupvar);
graph(2, 2).set_order_options('color', 0, 'column', 0);
graph(2, 2).no_legend();
graph(2, 2).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

%Draw figure
graph.set_title('Average Edited Spectra');
graph.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {}, 'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
graph.draw();
saveas(gcf, [resultspath, filesep, filename, '_AverageEditedSpectra.png'])
close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Spectral Quality for GABA and Glx (REST only)

clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figuresize(4) = figuresize(4)/1;
figure('Position', figuresize);

load([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specqual.mat']);
load([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specfreq.mat']);
specqual(rejectfilter, :) = [];
speclimits = find(specfreq>0 & specfreq<5);
gabalimits = find(specfreq>2.75 & specfreq<3.25);
glxlimits = find(specfreq>3.45 & specfreq<4.1);

clear graph

%Plot GABA
graph(1, 1) = gramm('x', specfreq(speclimits), 'y', specqual(:, speclimits), 'subset', contains(groupvar, 'Rest'));
graph(1, 1).geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','k--');

graph(1, 1).stat_summary();
graph(1, 1).facet_grid([], groupvar);
graph(1, 1).set_order_options('color', 0, 'column', 0);
graph(1, 1).no_legend();
graph(1, 1).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

graph(2, 1) = gramm('x', specfreq(gabalimits), 'y', specqual(:, gabalimits), 'subset', contains(groupvar, 'Rest'));
graph(2, 1).stat_summary();
graph(2, 1).facet_grid([], groupvar);
graph(2, 1).set_order_options('color', 0, 'column', 0);
graph(2, 1).no_legend();
graph(2, 1).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

%Plot Glx
graph(1, 2) = gramm('x', specfreq(speclimits), 'y', specqual(:, speclimits), 'subset', contains(groupvar, 'Rest'));
graph(1, 2).geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','k--');

graph(1, 2).stat_summary();
graph(1, 2).facet_grid([], groupvar);
graph(1, 2).set_order_options('color', 0, 'column', 0);
graph(1, 2).no_legend();
graph(1, 2).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

graph(2, 2) = gramm('x', specfreq(glxlimits), 'y', specqual(:, glxlimits), 'subset', contains(groupvar, 'Rest'));
graph(2, 2).stat_summary();
graph(2, 2).facet_grid([], groupvar);
graph(2, 2).set_order_options('color', 0, 'column', 0);
graph(2, 2).no_legend();
graph(2, 2).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

%Draw figure
graph.set_title('Average Edited Spectra');
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {}, 'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
graph.draw();
saveas(gcf, [resultspath, filesep, filename, '_AverageEditedSpectraRestOnly.png'])
close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save current state of matlab workspace
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

switch fittingmethod
    case 1
        rmpath(genpath(ospreypath));
    case 2
        rmpath(genpath(gannetpath));
end

return;