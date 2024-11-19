% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20240504
% ########################################################################

function lnifmri_mrs_dynamic_analysis(fittingmethod, keyword, wsize, wstep)

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
resultsdirname = 'RESULTS_dynamic_analysis';
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
% model = 'unused';
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

inputpars = table(string(fitpackage{fittingmethod}), string(keyword), wsize, wstep, 'VariableNames', ["package", "keyword", "wsize", "wstep"]);
writetable(inputpars, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_InputParameters.csv']);

listsub = dir([curdir, filesep, keyword, '*']);
listdataset = dir([curdir, filesep, keyword, '*', filesep, 'MRI_MRS', filesep, '*_fMRS']);

% Load demographic data to be used as covariates
demographicdata = readtable([curdir, filesep, 'fMRI_fMRS_dataset_', keyword, '.csv']);
covariates = demographicdata(:, {'Sex', 'Age'});

namesmetab = cell(1, 31);
namesrow = cell(size(listdataset, 1), 1);
for k =1:size(namesrow, 1)
    namesrow(k) = {[num2str(k, '%04d'), '_NaN']};
end
namerun = {'Rest', '0Back', '1Back', '2Back'};

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
                slidingwindow = [wsize, wstep];
                swstruc = lnifmri_mrs_dynamic_osprey(jobname, path_metab, path_ref, slidingwindow, outputfolder);
                % Append subject code and load
                subjcode = repmat(s,   [size(swstruc.diff, 1),   1]);
                loadcode = repmat(n,   [size(swstruc.diff, 1),   1]);
                
                codes = [subjcode loadcode];
                TableCodes = array2table(codes, 'VariableNames', {'SubCode' 'Load'});
                                
                swstruc.diff = [swstruc.diff, codes];
                swstruc.off = [swstruc.off, codes];
                swstruc.qual = [swstruc.qual, codes];
                
                swstruc.TableDiff = [swstruc.TableDiff, TableCodes];
                swstruc.TableOff = [swstruc.TableOff, TableCodes];
                swstruc.TableQual = [swstruc.TableQual, TableCodes];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 2 % Run GANNET
                addpath(genpath(gannetpath));
                addpath(genpath(exportfigpath));
                slidingwindow = [wsize, wstep];
                swstruc = lnifmri_mrs_dynamic_gannet(path_metab, path_ref, slidingwindow, gannetpath, outputfolder);
                % Append subject code and load
                subjcode = repmat(s,   [size(swstruc.diff, 1),   1]);
                loadcode = repmat(n,   [size(swstruc.diff, 1),   1]);
                
                codes = [subjcode loadcode];
                TableCodes = array2table(codes, 'VariableNames', {'SubCode' 'Load'});
                                
                swstruc.diff = [swstruc.diff, codes];
                swstruc.off = [swstruc.off, codes];
                swstruc.qual = [swstruc.qual, codes];
                
                swstruc.TableDiff = [swstruc.TableDiff, TableCodes];
                swstruc.TableOff = [swstruc.TableOff, TableCodes];
                swstruc.TableQual = [swstruc.TableQual, TableCodes];
        end%switch
        
        if exist('swmetabstruc', 'var')
            swmetabstruc = [swmetabstruc, swstruc];
        else
            swmetabstruc = swstruc;
        end
        
        clear swstruc
        
    end%listruns
   save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);
       
end%listsub
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract data for stats and plots

%Load data structure to plot
figurediff = swmetabstruc(1).diff;
for s = 2:size(swmetabstruc, 2)
    figurediff = [figurediff; swmetabstruc(s).diff];
end
%Extract number of subjects
subnum = numel(unique(figurediff(:, 5)));
%Extract number of conditions (levels: Rest - 0B, 1B, 2B)
[condnum, ~, ic] = unique(figurediff(:, 6));
%Extract number of sliding window frames for each condition
framenum = accumarray(ic, 1)./subnum;
frame = [];
%Define condition labels
for r = 1:numel(condnum)
    frame = [frame, 1:framenum(r)];
end
%Define sliding window frame vector
timeframe = repmat(frame', [subnum 1]);

%Extract baseline value from each condition, first window frame
baselineFilter = false(size(timeframe));
baselineFilter(1:size(frame, 2):end) = true;
for i=1:size(framenum, 1)
    baselineFilter(sum(framenum(1:i))+1:size(frame, 2):end) = true;
end
baselineRef = figurediff(baselineFilter, :);
%Expand reference values for subsequent normalization
%baselineRef = repelem(baselineRef, size(frame, 2), 1);
baselineRef = repelem(baselineRef, repmat(framenum, subnum, 1), 1);

%Store an un-normalized version to compute EIB
raw_figurediff = figurediff;

%Normalize all subsequent values by baseline: Adjust baseline value to 0 for each condition
figurediff(:, 1:4) = (figurediff(:, 1:4)-baselineRef(:, 1:4))./baselineRef(:, 1:4);

%Concatenate metabolite data with sliding window frame vector
figurediff = [figurediff, timeframe];

%Repeat for EIB
%Compute EIB (Excitation / Inhibition Ratio)
EIBdiff = raw_figurediff(:, 4)./raw_figurediff(:, 2); 
EIBRef = EIBdiff(baselineFilter);
EIBRef = repelem(EIBRef, repmat(framenum, subnum, 1), 1);

%Normalize all subsequent values by baseline
EIBdiff = (EIBdiff - EIBRef)./EIBRef;

%Concatenate EIB data with condition index and sliding window frame vector
condition = figurediff(:, 6);
EIBdiff = [EIBdiff, condition, timeframe];

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%% Compute temporal features of GABA curve
GABAkinetic = array2table(figurediff(:, [2, 6, 7]), 'VariableNames', {'GABA', 'Load', 'Time'});
GABAsubjlabel = array2table(figurediff(:, 5), 'VariableNames', {'Subject'});
GABAkinetic = [GABAkinetic, GABAsubjlabel];

checkfilter = zeros(size(GABAkinetic.GABA, 1), 1);
for s = 1:size(unique(GABAkinetic.Subject), 1)
    for l = 1:size(unique(GABAkinetic.Load), 1)
    
    rowfilter = logical((GABAkinetic.Subject == s) .* (GABAkinetic.Load == l));
    GABAdata = GABAkinetic{rowfilter, "GABA"};
    GABAcheckfilter = [checkfilter, rowfilter];    
    % Smooth outliers
    GABAsmooth = smoothdata(GABAdata, "rlowess", 5);
    % Remove vertical shift introduced by smoothing (baseline correction)
    GABAdata = GABAsmooth - GABAsmooth(1);
    % Store smoothed value
    GABAkinetic{rowfilter, "GABA"} = GABAdata;

    % Compute AUC
    x = 1:numel(GABAdata);
    GABA_auc = cumtrapz(x, GABAdata);
    GABA_aucFunc = @(a,b) max(GABA_auc(x<=b)) - min(GABA_auc(x>=a));
    GABA_tf_auc(s, l) = GABA_aucFunc(x(1), x(end));

    % Compute TTP
    %[M, TTP] = max(eibdata - abs(eibdata));
    %EIB_tf_ttp(s, l) = TTP;
    info = stepinfo(GABAdata);
    GABA_tf_ttp(s, l) = info.PeakTime;

    % Compute slope
    [p1,~,mu1] = polyfit(x([1 info.PeakTime])', GABAdata([1 info.PeakTime]), 1);
    %[p2,~,mu2] = polyfit(x(1:info.PeakTime)', GABAdata(1:info.PeakTime), 1);
    GABA_tf_slope(s, l) = p1(1);

    % Compute zerocrossing
    zerox = 0;
    for i=1:size(GABAdata)-1
        % Test if 0 is contained in the (sorted) interval between
        % two adjacent points
        if ~isnan(discretize(0, sort([GABAdata(i), GABAdata(i+1)])))
            % If true: crossing detected, increase the counter
            zerox = zerox+1;
        end
    end
    GABA_tf_zcross(s, l) = zerox;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get temporal properties via visibility graph
    % https://doi.org/10.1073/pnas.0709247105
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transform time series to visibility graph
    tsfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_GABAdata_s', num2str(s), '_l', num2str(l), '.dat'];
    vgfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_GABAvisgraph_s', num2str(s), '_l', num2str(l), '.dat'];
    writematrix(GABAdata, tsfilename);
    % Define ts2vg command line call - output edge list of directed graph
    call_ts2vg = ['ts2vg ', tsfilename, ' -t horizontal -d left_to_right -o ', vgfilename];
    status = system(call_ts2vg);
    % Load graph edge list
    % Column 1 -> Source Node
    % Column 2 -> Destination Node
    GABA_vg_edges = readmatrix(vgfilename);
    % Adjust node indices in the range [1 ... n];
    GABA_vg_edges = GABA_vg_edges+1;
    GABA_vg_numnodes = numel(GABAdata);

%     % Build Adjacency Matrix (UPPER TRIANGULAR)
%     EIB_vg_adjmtxUT = sparse(EIB_vg_edges(:, 1), EIB_vg_edges(:, 2), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxUT = full(EIB_vg_adjmtxUT);
%     % Build Adjacency Matrix (LOWER TRIANGULAR)
%     EIB_vg_adjmtxLT = sparse(EIB_vg_edges(:, 2), EIB_vg_edges(:, 1), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxLT = full(EIB_vg_adjmtxLT);
%     % Build Adjacency Matrix (SYMMETRIC)
%     EIB_vg_adjmtx = EIB_vg_adjmtxLT + EIB_vg_adjmtxUT;
    
    % Compute in_degree sequence
    degree_seq_in = zeros(GABA_vg_numnodes, 1);
    % Extract unique Destinations
    [N, ia, ic_in] = unique(GABA_vg_edges(:, 2));
    % Compute k_in (number of times that a node is a destination)
    ds_in = accumarray(ic_in, 1);
    degree_seq_in(N) = ds_in;
    
    % Compute in_degree distribution
    degree_dist_in = zeros(GABA_vg_numnodes, 2);
    [K, ia, ic_in] = unique(degree_seq_in);
    dd_in = accumarray(ic_in, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_in(K+1, :) = [dd_in, dd_in./GABA_vg_numnodes];

    % Compute out_degree sequence
    degree_seq_out = zeros(GABA_vg_numnodes, 1);
    % Extract unique Sources
    [N, ia, ic_out] = unique(GABA_vg_edges(:, 1));
    % Compute k_out (number of times a node is a source)
    ds_out = accumarray(ic_out, 1);
    degree_seq_out(N) = ds_out;
    GABA_vg_degree_out(s, l) = mean(degree_seq_out);

    % Compute out_degree distribution
    degree_dist_out = zeros(GABA_vg_numnodes, 2);
    [K, ia, ic_out] = unique(degree_seq_out);
    dd_out = accumarray(ic_out, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_out(K+1,:) = [dd_out, dd_out./GABA_vg_numnodes];
    
    % Compute stationarity via irreversibility of directed VG
    % (https://arxiv.org/pdf/1510.01318.pdf)
    P = degree_dist_out(:,2);
    Q = degree_dist_in(:, 2);
    ind = (P > 0 & Q > 0);
    PP = P(ind)./sum(P(ind));
    QQ = Q(ind)./sum(Q(ind));
    kld = PP.*(log10(PP)-log10(QQ));
    ind = ~isinf(kld);
    Dkld = sum(kld(ind), 'omitnan');
    GABA_vg_stationarity(s, l) = Dkld;

    end%loop on Load
end%loop on Subject

GABA_tf_AUC = array2table(GABA_tf_auc, "VariableNames", namerun);
GABA_tf_TTP = array2table(GABA_tf_ttp, "VariableNames", namerun);
GABA_tf_SLOPE = array2table(GABA_tf_slope, "VariableNames", namerun);
GABA_tf_ZCROSS = array2table(GABA_tf_zcross, "VariableNames", namerun);

% Temporal features extracted from Visibility Graph
GABA_vg_DEGREE_OUT = [subjTable array2table(GABA_vg_degree_out, "VariableNames", namerun)];
GABA_vg_STATIONARITY = [subjTable array2table(GABA_vg_stationarity, "VariableNames", namerun)];

writetable(GABA_tf_AUC, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_temporalfeatures_AUC.csv']);
writetable(GABA_tf_TTP, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_temporalfeatures_TTP.csv']);
writetable(GABA_tf_SLOPE, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_temporalfeatures_SLOPE.csv']);
writetable(GABA_tf_ZCROSS, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_temporalfeatures_ZCROSS.csv']);
% Temporal features extracted from Visibility Graph
writetable(GABA_vg_DEGREE_OUT, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_visibilitygraph_DEGREE_OUT.csv']);
writetable(GABA_vg_STATIONARITY, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_GABA_visibilitygraph_STATIONARITY.csv']);

GABAsummary = groupsummary(GABAkinetic, ["Load"], "mean", "GABA");
writetable(GABAsummary,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_Mean_GABA_x_Load.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute temporal features of Glx curve
Glxkinetic = array2table(figurediff(:, [4, 6, 7]), 'VariableNames', {'Glx', 'Load', 'Time'});
Glxsubjlabel = array2table(figurediff(:, 5), 'VariableNames', {'Subject'});
Glxkinetic = [Glxkinetic, Glxsubjlabel];

checkfilter = zeros(size(Glxkinetic.Glx, 1), 1);
for s = 1:size(unique(Glxkinetic.Subject), 1)
    for l = 1:size(unique(Glxkinetic.Load), 1)
    
    rowfilter = logical((Glxkinetic.Subject == s) .* (Glxkinetic.Load == l));
    Glxcheckfilter = [checkfilter, rowfilter];
        
    Glxdata = Glxkinetic{rowfilter, "Glx"};
    
    % Smooth outliers
    Glxsmooth = smoothdata(Glxdata, "rlowess", 5);
    % Remove vertical shift introduced by smoothing (baseline correction)
    Glxdata = Glxsmooth - Glxsmooth(1);
    % Store smoothed value
    Glxkinetic{rowfilter, "Glx"} = Glxdata;

    % Compute AUC
    x = 1:numel(Glxdata);
    Glx_auc = cumtrapz(x, Glxdata);
    Glx_aucFunc = @(a,b) max(Glx_auc(x<=b)) - min(Glx_auc(x>=a));
    Glx_tf_auc(s, l) = Glx_aucFunc(x(1), x(end));

    % Compute TTP
    %[M, TTP] = max(eibdata - abs(eibdata));
    %EIB_tf_ttp(s, l) = TTP;
    info = stepinfo(Glxdata);
    Glx_tf_ttp(s, l) = info.PeakTime;

    % Compute slope
    [p1,~,mu1] = polyfit(x([1 info.PeakTime])', Glxdata([1 info.PeakTime]), 1);
    %[p2,~,mu2] = polyfit(x(1:info.PeakTime)', Glxdata(1:info.PeakTime), 1);
    Glx_tf_slope(s, l) = p1(1);
    
    % Compute zerocrossing
    zerox = 0;
    for i=1:size(Glxdata)-1
        % Test if 0 is contained in the (sorted) interval between
        % two adjacent points
        if ~isnan(discretize(0, sort([Glxdata(i), Glxdata(i+1)])))
            % If true: crossing detected, increase the counter
            zerox = zerox+1;
        end
    end
    Glx_tf_zcross(s, l) = zerox;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get temporal properties via visibility graph
    % https://doi.org/10.1073/pnas.0709247105
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transform time series to visibility graph
    tsfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_Glxdata_s', num2str(s), '_l', num2str(l), '.dat'];
    vgfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_Glxvisgraph_s', num2str(s), '_l', num2str(l), '.dat'];
    writematrix(Glxdata, tsfilename);
    % Define ts2vg command line call - output edge list of directed graph
    call_ts2vg = ['ts2vg ', tsfilename, ' -t horizontal -d left_to_right -o ', vgfilename];
    status = system(call_ts2vg);
    % Load graph edge list
    % Column 1 -> Source Node
    % Column 2 -> Destination Node
    Glx_vg_edges = readmatrix(vgfilename);
    % Adjust node indices in the range [1 ... n];
    Glx_vg_edges = Glx_vg_edges+1;
    Glx_vg_numnodes = numel(Glxdata);

%     % Build Adjacency Matrix (UPPER TRIANGULAR)
%     EIB_vg_adjmtxUT = sparse(EIB_vg_edges(:, 1), EIB_vg_edges(:, 2), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxUT = full(EIB_vg_adjmtxUT);
%     % Build Adjacency Matrix (LOWER TRIANGULAR)
%     EIB_vg_adjmtxLT = sparse(EIB_vg_edges(:, 2), EIB_vg_edges(:, 1), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxLT = full(EIB_vg_adjmtxLT);
%     % Build Adjacency Matrix (SYMMETRIC)
%     EIB_vg_adjmtx = EIB_vg_adjmtxLT + EIB_vg_adjmtxUT;
    
    % Compute in_degree sequence
    degree_seq_in = zeros(Glx_vg_numnodes, 1);
    % Extract unique Destinations
    [N, ia, ic_in] = unique(Glx_vg_edges(:, 2));
    % Compute k_in (number of times that a node is a destination)
    ds_in = accumarray(ic_in, 1);
    degree_seq_in(N) = ds_in;
    
    % Compute in_degree distribution
    degree_dist_in = zeros(Glx_vg_numnodes, 2);
    [K, ia, ic_in] = unique(degree_seq_in);
    dd_in = accumarray(ic_in, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_in(K+1, :) = [dd_in, dd_in./Glx_vg_numnodes];

    % Compute out_degree sequence
    degree_seq_out = zeros(Glx_vg_numnodes, 1);
    % Extract unique Sources
    [N, ia, ic_out] = unique(Glx_vg_edges(:, 1));
    % Compute k_out (number of times a node is a source)
    ds_out = accumarray(ic_out, 1);
    degree_seq_out(N) = ds_out;
    Glx_vg_degree_out(s, l) = mean(degree_seq_out);

    % Compute out_degree distribution
    degree_dist_out = zeros(Glx_vg_numnodes, 2);
    [K, ia, ic_out] = unique(degree_seq_out);
    dd_out = accumarray(ic_out, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_out(K+1,:) = [dd_out, dd_out./Glx_vg_numnodes];
    
    % Compute stationarity via irreversibility of directed VG
    % (https://arxiv.org/pdf/1510.01318.pdf)
    P = degree_dist_out(:,2);
    Q = degree_dist_in(:, 2);
    ind = (P > 0 & Q > 0);
    PP = P(ind)./sum(P(ind));
    QQ = Q(ind)./sum(Q(ind));
    kld = PP.*(log10(PP)-log10(QQ));
    ind = ~isinf(kld);
    Dkld = sum(kld(ind), 'omitnan');
    Glx_vg_stationarity(s, l) = Dkld;

    end%loop on Load
end%loop on Subject

Glx_tf_AUC = array2table(Glx_tf_auc, "VariableNames", namerun);
Glx_tf_TTP = array2table(Glx_tf_ttp, "VariableNames", namerun);
Glx_tf_SLOPE = array2table(Glx_tf_slope, "VariableNames", namerun);
Glx_tf_ZCROSS = array2table(Glx_tf_zcross, "VariableNames", namerun);

% Temporal features extracted from Visibility Graph
Glx_vg_DEGREE_OUT = [subjTable array2table(Glx_vg_degree_out, "VariableNames", namerun)];
Glx_vg_STATIONARITY = [subjTable array2table(Glx_vg_stationarity, "VariableNames", namerun)];

writetable(Glx_tf_AUC, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_temporalfeatures_AUC.csv']);
writetable(Glx_tf_TTP, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_temporalfeatures_TTP.csv']);
writetable(Glx_tf_SLOPE, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_temporalfeatures_SLOPE.csv']);
writetable(Glx_tf_ZCROSS, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_temporalfeatures_ZCROSS.csv']);
% Temporal features extracted from Visibility Graph
writetable(Glx_vg_DEGREE_OUT, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_visibilitygraph_DEGREE_OUT.csv']);
writetable(Glx_vg_STATIONARITY, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_Glx_visibilitygraph_STATIONARITY.csv']);

Glxsummary = groupsummary(Glxkinetic, ["Load"], "mean", "Glx");
writetable(Glxsummary,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_Mean_Glx_x_Load.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute temporal features of EIB curve
EIBkinetic = array2table(EIBdiff, 'VariableNames', {'EIB', 'Load', 'Time'});
EIBsubjlabel = array2table(figurediff(:, 5), 'VariableNames', {'Subject'});
EIBkinetic = [EIBkinetic, EIBsubjlabel];

checkfilter = zeros(size(EIBkinetic.EIB, 1), 1);
for s = 1:size(unique(EIBkinetic.Subject), 1)
    for l = 1:size(unique(EIBkinetic.Load), 1)
    
    rowfilter = logical((EIBkinetic.Subject == s) .* (EIBkinetic.Load == l));
    EIBcheckfilter = [checkfilter, rowfilter];
        
    EIBdata = EIBkinetic{rowfilter, "EIB"};
    
    % Smooth outliers
    eibsmooth = smoothdata(EIBdata, "rlowess", 5);
    % Remove vertical shift introduced by smoothing (baseline correction)
    EIBdata = eibsmooth - eibsmooth(1);
    % Store smoothed value
    EIBkinetic{rowfilter, "EIB"} = EIBdata;

    % Compute AUC
    x = 1:numel(EIBdata);
    EIB_auc = cumtrapz(x, EIBdata);
    EIB_tf_auc_std(s, l) = trapz(diff(EIB_auc));
    EIB_aucFunc = @(a,b) max(EIB_auc(x<=b)) - min(EIB_auc(x>=a));
    EIB_tf_auc(s, l) = EIB_aucFunc(x(1), x(end));

    % Compute different temporal features
    info = lsiminfo(EIBdata, x);    
    EIB_tf_ttp(s, l) = info.MaxTime/numel(EIBdata);
    %EIB_tf_transient(s, l) = info.TransientTime;
    %EIB_tf_settling(s, l) = info.SettlingTime;

    % Compute slope
    [p1,~,mu1] = polyfit(x([1 info.MaxTime])', EIBdata([1 info.MaxTime]), 1);
    %[p2,~,mu2] = polyfit(x(1:info.MaxTime)', EIBdata(1:info.MaxTime), 1);
    EIB_tf_slope(s, l) = p1(1);
    f1 = polyval(p1,x',[],mu1);
    %f2 = polyval(p2,x',[],mu2);
%     hold on
%     plot(x',EIBdata,'o--')
%     plot(x(2:end)',diff(EIBdata),'k-')
%     plot(x',f1, 'LineWidth', 2)
%     %plot(x',f2)
%     hold off
%     ax = gca;
%     exportgraphics(ax,['EIBslopeTEST_s_', num2str(s), '_l_', num2str(l), '.png'],'Resolution',300)
%     close (gcf);

    % Compute zerocrossing
    zerox = 0;
    for i=1:size(EIBdata)-1
        % Test if 0 is contained in the (sorted) interval between
        % two adjacent points
        if ~isnan(discretize(0, sort([EIBdata(i), EIBdata(i+1)])))
            % If true: crossing detected, increase the counter
            zerox = zerox+1;
        end
    end
    EIB_tf_zcross(s, l) = zerox;

    % Compute midcross
    mx = midcross(EIBdata, x, MidPercentReferenceLevel=75, Statelevels=[min(EIBdata) max(EIBdata)]);
    EIB_tf_midcross(s, l) = mx(1)/numel(EIBdata);
%     %Only for visualization
%     midcross(EIBdata, x, MidPercentReferenceLevel=75, Statelevels=[min(EIBdata) max(EIBdata)]);
%     ax = gca;
%     exportgraphics(ax,['EIBmidcrossTEST_s_', num2str(s), '_l_', num2str(l), '.png'],'Resolution',300)
%     close (gcf);

    % Test stationarity
    h = kpsstest(EIBdata);
    EIB_tf_stationarity(s, l) = h;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get temporal properties via visibility graph
    % https://doi.org/10.1073/pnas.0709247105
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transform time series to visibility graph
    tsfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_EIBdata_s', num2str(s), '_l', num2str(l), '.dat'];
    vgfilename = [resultspath, filesep, timestamp,...
                '_', keyword, '_', fitpackage{fittingmethod},...
                '_EIBvisgraph_s', num2str(s), '_l', num2str(l), '.dat'];
    writematrix(EIBdata, tsfilename);
    % Define ts2vg command line call - output edge list of directed graph
    call_ts2vg = ['ts2vg ', tsfilename, ' -t horizontal -d left_to_right -o ', vgfilename];
    status = system(call_ts2vg);
    % Load graph edge list
    % Column 1 -> Source Node
    % Column 2 -> Destination Node
    EIB_vg_edges = readmatrix(vgfilename);
    % Adjust node indices in the range [1 ... n];
    EIB_vg_edges = EIB_vg_edges+1;
    EIB_vg_numnodes = numel(EIBdata);

%     % Build Adjacency Matrix (UPPER TRIANGULAR)
%     EIB_vg_adjmtxUT = sparse(EIB_vg_edges(:, 1), EIB_vg_edges(:, 2), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxUT = full(EIB_vg_adjmtxUT);
%     % Build Adjacency Matrix (LOWER TRIANGULAR)
%     EIB_vg_adjmtxLT = sparse(EIB_vg_edges(:, 2), EIB_vg_edges(:, 1), 1, EIB_vg_numnodes, EIB_vg_numnodes);
%     EIB_vg_adjmtxLT = full(EIB_vg_adjmtxLT);
%     % Build Adjacency Matrix (SYMMETRIC)
%     EIB_vg_adjmtx = EIB_vg_adjmtxLT + EIB_vg_adjmtxUT;
    
    % Compute in_degree sequence
    degree_seq_in = zeros(EIB_vg_numnodes, 1);
    % Extract unique Destinations
    [N, ia, ic_in] = unique(EIB_vg_edges(:, 2));
    % Compute k_in (number of times that a node is a destination)
    ds_in = accumarray(ic_in, 1);
    degree_seq_in(N) = ds_in;
    degree_seq_in_all(N, s, l) = ds_in;

    % Compute in_degree distribution
    degree_dist_in = zeros(EIB_vg_numnodes, 2);
    [K, ia, ic_in] = unique(degree_seq_in);
    dd_in = accumarray(ic_in, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_in(K+1, :) = [dd_in, dd_in./EIB_vg_numnodes];

    % Compute out_degree sequence
    degree_seq_out = zeros(EIB_vg_numnodes, 1);
    % Extract unique Sources
    [N, ia, ic_out] = unique(EIB_vg_edges(:, 1));
    % Compute k_out (number of times a node is a source)
    ds_out = accumarray(ic_out, 1);
    degree_seq_out(N) = ds_out;
    degree_seq_out_all(N, s, l) = ds_out;

    EIB_vg_degree_out(s, l) = mean(degree_seq_out);

    % Compute out_degree distribution
    degree_dist_out = zeros(EIB_vg_numnodes, 2);
    [K, ia, ic_out] = unique(degree_seq_out);
    dd_out = accumarray(ic_out, 1);
    % Compose: unique degree (k), number of nodes with degree k, P(k)
    degree_dist_out(K+1,:) = [dd_out, dd_out./EIB_vg_numnodes];
    
    % Compute stationarity via irreversibility of directed VG
    % (https://arxiv.org/pdf/1510.01318.pdf)
    P = degree_dist_out(:,2);
    Q = degree_dist_in(:, 2);
    ind = (P > 0 & Q > 0);
    PP = P(ind)./sum(P(ind));
    QQ = Q(ind)./sum(Q(ind));
    kld = PP.*(log10(PP)-log10(QQ));
    ind = ~isinf(kld);
    Dkld = sum(kld(ind), 'omitnan');
    EIB_vg_stationarity(s, l) = Dkld;

    end%loop on Load
end%loop on Subject

subjTable = unique(subjlabel);
EIB_tf_AUC = [subjTable array2table(EIB_tf_auc, "VariableNames", namerun)];
EIB_tf_TTP = [subjTable array2table(EIB_tf_ttp, "VariableNames", namerun)];
EIB_tf_SLOPE = [subjTable array2table(EIB_tf_slope, "VariableNames", namerun)];
EIB_tf_ZCROSS = [subjTable array2table(EIB_tf_zcross, "VariableNames", namerun)];
EIB_tf_STATIONARITY = [subjTable array2table(EIB_tf_stationarity, "VariableNames", namerun)];
% Temporal features extracted from Visibility Graph
EIB_vg_DEGREE_OUT = [subjTable array2table(EIB_vg_degree_out, "VariableNames", namerun)];
EIB_vg_STATIONARITY = [subjTable array2table(EIB_vg_stationarity, "VariableNames", namerun)];

writetable(EIB_tf_AUC, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_temporalfeatures_AUC.csv']);
writetable(EIB_tf_TTP, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_temporalfeatures_TTP.csv']);
writetable(EIB_tf_SLOPE, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_temporalfeatures_SLOPE.csv']);
writetable(EIB_tf_ZCROSS, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_temporalfeatures_ZCROSS.csv']);
writetable(EIB_tf_STATIONARITY, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_temporalfeatures_STATIONARITY.csv']);
% Temporal features extracted from Visibility Graph
writetable(EIB_vg_DEGREE_OUT, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_visibilitygraph_DEGREE_OUT.csv']);
writetable(EIB_vg_STATIONARITY, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_EIB_visibilitygraph_STATIONARITY.csv']);

EIBsummary = groupsummary(EIBkinetic, ["Load"], "mean", "EIB");
writetable(EIBsummary,...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_Mean_EIB_x_Load.csv']);

writetable(EIBkinetic(:, 2:end),...
            [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_LoadTimeSubjectLabels.csv']);

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stats - One-way anova on AUC of GABA, Glx and EIB curve
factors = repmat(string(namerun), subnum, 1);

% Test for normality
ksGABA = kstest(GABA_tf_auc(:));
ksGlx = kstest(Glx_tf_auc(:));
ksEIB = kstest(EIB_tf_auc(:));

% Distribution is not normal, then compute
% Kruskall-Wallis nonparametric ANOVA
[pKW_GABA, tblKW_GABA, statsKW_GABA] =  kruskalwallis(GABA_tf_auc(:), factors(:));
[pKW_Glx, tblKW_Glx, statsKW_Glx] =     kruskalwallis(Glx_tf_auc(:), factors(:));
[pKW_EIB, tblKW_EIB, statsKW_EIB] =     kruskalwallis(EIB_tf_auc(:), factors(:));

multcompare(statsKW_GABA)
multcompare(statsKW_Glx)
multcompare(statsKW_EIB)

%% Plot ANOVA results
xvar = categorical(factors(:));
yvar = GABA_tf_auc(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'GABA (AUC)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('AUC of GABA');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_GABA_tf_AUC_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = Glx_tf_auc(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Glx (AUC)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('AUC of Glx');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_tf_AUC_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = EIB_tf_auc(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'EIB (AUC)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('AUC of EIB');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB_tf_AUC_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats - One-way anova on OUT-DEGREE of GABA, Glx and EIB curve
factors = repmat(string(namerun), subnum, 1);

ksGABA = kstest(GABA_vg_degree_out(:));
ksGlx = kstest(Glx_vg_degree_out(:));
ksEIB = kstest(EIB_vg_degree_out(:));
% Distribution is not normal, then compute
% Kruskall-Wallis nonparametric ANOVA
[pKW_GABA, tblKW_GABA, statsKW_GABA] =  kruskalwallis(GABA_vg_degree_out(:), factors(:));
[pKW_Glx, tblKW_Glx, statsKW_Glx] =     kruskalwallis(Glx_vg_degree_out(:), factors(:));
[pKW_EIB, tblKW_EIB, statsKW_EIB] =     kruskalwallis(EIB_vg_degree_out(:), factors(:));

multcompare(statsKW_GABA)
multcompare(statsKW_Glx)
multcompare(statsKW_EIB)

%% Plot ANOVA results
xvar = categorical(factors(:));
yvar = GABA_vg_degree_out(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Out Degree (a.u)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('GABA');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [0.75 1.75]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_GABA_vg_degree_out_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = Glx_vg_degree_out(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Out Degree (a.u)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('Glx');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [0.75 1.75]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_vg_degree_out_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = EIB_vg_degree_out(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Out Degree (a.u)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('EIB');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [0.75 1.75]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB_vg_degree_out_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats - One-way anova on STATIONARITY of GABA, Glx and EIB curve
factors = repmat(string(namerun), subnum, 1);

ksGABA = kstest(GABA_vg_stationarity(:));
ksGlx = kstest(Glx_vg_stationarity(:));
ksEIB = kstest(EIB_vg_stationarity(:));
% Distribution is not normal, then compute
% Kruskall-Wallis nonparametric ANOVA
[pKW_GABA, tblKW_GABA, statsKW_GABA] =  kruskalwallis(GABA_vg_stationarity(:), factors(:));
[pKW_Glx, tblKW_Glx, statsKW_Glx] =     kruskalwallis(Glx_vg_stationarity(:), factors(:));
[pKW_EIB, tblKW_EIB, statsKW_EIB] =     kruskalwallis(EIB_vg_stationarity(:), factors(:));

multcompare(statsKW_GABA)
multcompare(statsKW_Glx)
multcompare(statsKW_EIB)

%% Plot ANOVA results
xvar = categorical(factors(:));
yvar = GABA_vg_stationarity(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Stationarity (a.u.)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('GABA');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-0.02 0.2]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_GABA_vg_stationarity_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = Glx_vg_stationarity(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Stationarity (a.u.)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('Glx');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-0.02 0.2]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_vg_stationarity_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

xvar = categorical(factors(:));
yvar = EIB_vg_stationarity(:);
cvar = categorical(xvar);

graph = gramm('x', xvar, 'y', yvar, 'color', cvar);
graph.stat_boxplot();
graph.set_text_options('base_size', 14, 'label_scaling', 1.2, 'big_title_scaling', 1.4);
graph.set_order_options('color', 0, 'x', 0);
graph.set_names('x', '', 'y', 'Stationarity (a.u.)', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('EIB');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-0.02 0.2]);
graph.no_legend();
%Add individual points
graph.update('color', cvar);
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB_vg_stationarity_anova_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other plots
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figure('Position', figuresize);

% Plot GABA
graph = gramm('x', figurediff(:, 7), 'y', 100.*GABAkinetic.GABA, 'group', figurediff(:, 5));
graph.geom_point();
graph.geom_line();
graph.set_text_options('base_size', 16, 'label_scaling', 1.4, 'big_title_scaling', 1.6, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_color_options('chroma',0,'lightness',85); %We make it light grey
graph.set_names('x', 'Frame', 'y', '\Delta GABA', 'column', '', 'row', '');
graph.no_legend();
%graph.facet_grid([], namerun(EIBkinetic.Load));
graph.facet_wrap(namerun(GABAkinetic.Load), 'ncols', 2);
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-40 40]);
% figure('Position', figuresize.*5);
% graph.draw();
graph.update('color', GABAkinetic.Load, 'group', []);
graph.stat_summary('type', 'sem');
graph.set_color_options();
graph.set_line_options('base_size', 4);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_GABA'], 'file_type', 'png');
clear graph
close(gcf)

% Plot Glx
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figure('Position', figuresize);

graph = gramm('x', figurediff(:, 7), 'y', 100.*Glxkinetic.Glx, 'group', figurediff(:, 5));
graph.geom_point();
graph.geom_line();
graph.set_text_options('base_size', 16, 'label_scaling', 1.4, 'big_title_scaling', 1.6, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_color_options('chroma',0,'lightness',85); %We make it light grey
graph.no_legend();
%graph.facet_grid([], namerun(EIBkinetic.Load));
graph.facet_wrap(namerun(GABAkinetic.Load), 'ncols', 2);
graph.set_names('x', 'Frame', 'y', '\Delta Glx', 'column', '', 'row', '');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-40 40]);
% figure('Position', figuresize.*5);
% graph.draw();
graph.update('color', figurediff(:, 6), 'group', []);
graph.stat_summary('type', 'sem');
graph.set_color_options();
graph.set_line_options('base_size', 4);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx'], 'file_type', 'png');
clear graph
close(gcf)

% Plot EIB
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figure('Position', figuresize);

graph = gramm('x', EIBkinetic.Time, 'y', 100.*EIBkinetic.EIB, 'group', EIBkinetic.Subject);%, 'subset', EIBkinetic.Subject~=1 & EIBkinetic.Subject~=3 & EIBkinetic.Subject~=5);%, 'group', figurediff(:, 5));
graph.geom_point();
graph.geom_line();
graph.set_text_options('base_size', 16, 'label_scaling', 1.4, 'big_title_scaling', 1.6, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_color_options('chroma',0,'lightness',85); %We make it light grey
graph.no_legend();
%graph.facet_grid([], namerun(EIBkinetic.Load));
graph.facet_wrap(namerun(GABAkinetic.Load), 'ncols', 2);
graph.set_names('x', 'Frame', 'y', '\Delta EIB', 'column', '', 'row', '');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-40 40]);
%figure('Position', figuresize.*5);
%graph.draw();
graph.update('color', figurediff(:, 6), 'group', []);
%graph.geom_point();
graph.stat_summary('type', 'sem');
%graph.stat_smooth();
graph.set_color_options();
graph.set_line_options('base_size', 4);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB'], 'file_type', 'png');
clear graph
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic plots
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

% Plot EIB x Load
graph = gramm('x', categorical(string(namerun(EIBkinetic.Load))), 'y', 100.*EIBkinetic.EIB, 'color', namerun(EIBkinetic.Load));
graph.stat_boxplot();
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.set_order_options('color', 0, 'x', namerun);
graph.set_names('x', '', 'y', '\Delta EIB', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('\Delta EIB x Load');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5]);
graph.no_legend();
%Add individual points
graph.update('color', namerun(figurediff(:, 6)));
graph.set_point_options('base_size',6);
graph.set_color_options();
graph.geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB_boxplot'], 'file_type', 'png');
clear graph
close(gcf)

% Plot EIB x Frame
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

% % #####################################################################################################################
% % 20240802 - stefano.tambalo@unitn.it
% % Restored the time scale in frames
% % #####################################################################################################################
% % 20240725 - stefano.tambalo@unitn.it
% % Only for graphical purposes
% % Replace Number of frames in x-axis with time scale from sliding window
% % This needs to be fixed and automatically calculated in an updated
% % code refactoring - here is hard coded for semplicity:
% timescale = duration(0,[60*4:12:448]./60, 0, 'Format', 'mm:ss')';
% timerest = timescale(1:12);
% timeexp = [timerest; timescale; timescale; timescale];
% kinetictime = repmat(timeexp, 1, 12);
% graph = gramm('x', categorical(kinetictime), 'y', 100.*EIBkinetic.EIB, 'color', namerun(EIBkinetic.Load));
% % 20240725 - stefano.tambalo@unitn.it
% % #####################################################################################################################

graph = gramm('x', EIBkinetic.Time, 'y', 100.*EIBkinetic.EIB, 'color', namerun(EIBkinetic.Load));
graph.stat_summary('type', 'sem');
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_names('x', 'Frames', 'y', '\Delta EIB', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('\Delta EIB x Frame');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-10 20]);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_EIB_timeresponsecurve'], 'file_type', 'png');
clear graph
close(gcf)

% Plot Glx x Frame
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

graph = gramm('x', Glxkinetic.Time, 'y', 100.*Glxkinetic.Glx, 'color', namerun(Glxkinetic.Load));
graph.stat_summary('type', 'sem');
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_names('x', 'Frame', 'y', '\Delta Glx', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('\Delta Glx x Frame');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-15 15]);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_timeresponsecurve'], 'file_type', 'png');
clear graph
close(gcf)

% Plot GABA x Frame
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/2;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

graph = gramm('x', GABAkinetic.Time, 'y', 100.*GABAkinetic.GABA, 'color', namerun(GABAkinetic.Load));
graph.stat_summary('type', 'sem');
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.set_order_options('color', 0, 'column', 0);
graph.set_names('x', 'Frame', 'y', '\Delta GABA', 'column', '', 'row', '', 'color', 'Load');
graph.set_title('\Delta GABA x Frame');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-15 15]);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_GABA_timeresponsecurve'], 'file_type', 'png');
clear graph
close(gcf)

% Plot GABA~Glx GLM
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/1.75;
figure('Position', figuresize);

graph = gramm('x', 100.*GABAkinetic.GABA, 'y', 100.*Glxkinetic.Glx, 'color', namerun(GABAkinetic.Load));
graph.stat_glm('fullrange',1);
graph.set_text_options('base_size', 16, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.facet_wrap(namerun(GABAkinetic.Load), 'ncols', 2);
graph.set_order_options('color', 0, 'column', namerun);
graph.set_names('x', '\Delta GABA', 'y', '\Delta Glx', 'column', '', 'row', '');
graph.set_title('\Delta Glx ~ \Delta GABA');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-10 10], 'XLim', [-20 20]);
graph.no_legend();
graph.update('lightness', GABAkinetic.Time);
graph.geom_point();
graph.set_point_options('base_size', 8);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_GABA_scatterplot'], 'file_type', 'png');
clear graph
close(gcf)

% Plot GABA~Glx GLM over Frames
clear figuresize;
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/1.75;
figure('Position', figuresize);

graph = gramm('x', 100.*GABAkinetic.GABA, 'y', 100.*Glxkinetic.Glx, 'color', namerun(GABAkinetic.Load), 'lightness', GABAkinetic.Time);
graph.stat_glm("geom", "line", "fullrange", 1);
graph.facet_wrap(namerun(GABAkinetic.Load), 'ncols', 2);
graph.set_order_options('color', 0, 'column', namerun);
graph.set_line_options("base_size", 3);
graph.set_text_options('base_size', 16, 'legend_title_scaling', 1, 'legend_scaling', 1, 'label_scaling', 1.2, 'big_title_scaling', 1.4, 'Interpreter','tex');
graph.set_names('x', '\Delta GABA', 'y', ' \Delta Glx', 'color', '', 'lightness', 'Frame', 'column', '', 'row', '');
graph.set_title('\Delta Glx ~ GABA over Frames');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5], 'YLim', [-10 10], 'XLim', [-20 20]);
graph.no_legend();
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_Glx_GABA_over_Frames_lineplot'], 'file_type', 'png');
clear graph
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Spectral Quality for GABA and Glx

%Load data structure to plot
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/1.75;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

figurequal = swmetabstruc(1).qual;
for s = 2:size(swmetabstruc, 2)
    figurequal = [figurequal; swmetabstruc(s).qual];
end
%Add time frame index
figurequal = [figurequal, timeframe];

varnames = [swmetabstruc(1).TableQual.Properties.VariableNames, 'Frame'];
spectralQualityTable = array2table(figurequal, "VariableNames", varnames);
writetable(spectralQualityTable, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_DynamicSpectralQuality.csv']);
spectralQualityTableSummary = grpstats(spectralQualityTable, 'Load');
writetable(spectralQualityTableSummary, [resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod},...
            '_DynamicSpectralQualitySummary.csv']);

save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

groupvar = namerun(figurequal(:, 10));
framevar = figurequal(:,11);
colorvar = namerun(figurequal(:, 10));

%Plot GABA SNR
graph(1, 1) = gramm('x', framevar, 'y', figurequal(:,2), 'color', colorvar);
%Define boxplot
graph(1, 1).stat_summary();
graph(1, 1).set_order_options('color', 0, 'column', 0, 'row', 0);
graph(1, 1).facet_grid([], groupvar);
graph(1, 1).set_names('x', '', 'y', 'GABA SNR', 'column', '', 'row', '', 'color', 'Load');
graph(1, 1).set_title('');
graph(1, 1).axe_property('YLim', [0 30]);
graph(1, 1).no_legend();

%Add individual points
graph(1, 1).update('color', colorvar);
graph(1, 1).set_point_options('base_size',6);
graph(1, 1).set_color_options();
graph(1, 1).geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph(1, 1).no_legend();

%Plot Glx SNR
graph(2, 1) = gramm('x', framevar, 'y', figurequal(:,5), 'color', colorvar);
%Define boxplot
graph(2, 1).stat_summary();
graph(2, 1).set_order_options('color', 0, 'column', 0, 'row', 0);
graph(2, 1).facet_grid([], groupvar, 'column_labels', 0);
graph(2, 1).set_names('x', 'Frames', 'y', 'Glx SNR', 'column', '', 'row', '', 'color', 'Load');
graph(2, 1).set_title('');
graph(2, 1).axe_property('YLim', [0 30]);
graph(2, 1).no_legend();

%Add individual points
graph(2, 1).update('color', colorvar);
graph(2, 1).set_point_options('base_size',6);
graph(2, 1).set_color_options();
graph(2, 1).geom_jitter('width', 0.4, 'dodge',0.9, 'alpha', 0.25);
graph(2, 1).no_legend();

graph.set_title('Spectral Quality');
graph.axe_property('TickDir','out','XGrid','on','Ygrid','on','GridColor', [0.5 0.5 0.5]);

graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_DynamicSpectralQuality'], 'file_type', 'png');
clear graph
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot AverageEditedSpectra
clear figuresize
figuresize = get(0,'ScreenSize');
figuresize(3) = figuresize(3)/1;
figuresize(4) = figuresize(4)/1.5;
figure('Position', figuresize);

for s = 1:size(swmetabstruc, 2)
    tmpspec = mean(swmetabstruc(s).specqual, 2);
    specqual(s, :) = tmpspec';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only for retrospective reconstruction of plots.
% This line is not necessary for analysis run from scratch after 30/10/23
specfreq = readmatrix([pwd, filesep, 'specfreq.csv']);
% stefano.tambalo@unitn.it - 30/10/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speclimits = find(specfreq>0 & specfreq<5);
gabalimits = find(specfreq>2.75 & specfreq<3.25);
glxlimits = find(specfreq>3.45 & specfreq<4.1);

groupvar = namerun(repmat([1 2 3 4], 1, size(specqual, 1)/size(namerun, 2))');

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
graph(2, 1).facet_grid([], groupvar, 'column_labels', 0);
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
graph(2, 2).facet_grid([], groupvar, 'column_labels', 0);
graph(2, 2).set_order_options('color', 0, 'column', 0);
graph(2, 2).no_legend();
graph(2, 2).set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');

%Draw figure
graph.set_title('Average Edited Spectra');
graph.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {}, 'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
graph.draw();
graph.export('file_name', [resultspath, filesep, outdir, '_AverageEditedSpectra'], 'file_type', 'png');
clear graph
close(gcf)

%Plot individual spectra
individualspeclimits = find(specfreq>1 & specfreq<4.2);
%Color mapping for the polygons
cmap = [1 0.5 0.5; % red (Glx)
    0.5 0.5 1];% blue (GABA)
i=0;
for s = 1:size(unique(EIBkinetic.Subject), 1)
    for l = 1:size(unique(EIBkinetic.Load), 1)
        i = i+1;
        graph = gramm('x', specfreq(individualspeclimits), 'y', specqual(i, individualspeclimits));
        graph.geom_line();
        graph.set_color_options('map', [0 0 0]);
        graph.geom_polygon('x',{[specfreq(glxlimits(1)), specfreq(glxlimits(end))];  [specfreq(gabalimits(1)), specfreq(gabalimits(end))]},'color', cmap, 'alpha', 0.3);
        graph.geom_vline('xintercept', [specfreq(glxlimits(1)), specfreq(glxlimits(end))], 'style','r--');
        graph.geom_vline('xintercept', [specfreq(gabalimits(1)), specfreq(gabalimits(end))], 'style','b--');
        graph.no_legend();
        graph.set_names('x', 'Freq.(ppm)', 'y', '', 'column', '', 'row', '');
        %Draw figure
        graph.set_title(['Individual Edited Spectrum (s:', num2str(s), ' l:', num2str(l), ')']);
        graph.axe_property('xdir', 'reverse', 'TickDir','out', 'YTickLabel', {},...
            'Ylim', [-6e-5, 6e-5], 'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5]);
        graph.draw();
        graph.export('file_name', [resultspath, filesep, outdir,'_IndividualEditedSpectrum_s', num2str(s), '_l_', num2str(l)], 'file_type', 'png');
        clear graph
        close(gcf)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping
% Save current state of matlab workspace
save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_workspace.mat']);

% Set folder to store .dat files created by ts2vg
datdirname = 'datfiles';
if ~isfolder([resultspath, filesep, datdirname])
    createdatdir = ['mkdir ', resultspath, filesep, datdirname];
    system(createdatdir);
end
datdir = [resultspath, filesep, datdirname];
movedatfiles = ['mv *.dat ', datdir];
system(movedatfiles);

% Remove fitting package from path
% switch fittingmethod
%    case 1
%        rmpath(genpath(ospreypath));
%    case 2
%        rmpath(genpath(gannetpath));
%        rmpath(genpath(exportfigpath));
% end

end%function
