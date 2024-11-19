% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20220504
% ########################################################################

function swstruc = lnifmri_mrs_dynamic_osprey(jobname, path_metab, path_ref, window, outputfolder)

rng default;
prefixmetab = 'swmetabfile';
ext = '.ima';
% Create dir to store raw MRS data
createrawdir = ['mkdir ', outputfolder, filesep, prefixmetab];
system(createrawdir);
% Copy raw MRS data
rawfiles = dir([path_metab, filesep, '*.IMA']);
rawnumber = size(rawfiles, 1);

%define a sliding window filter
swf = lnifmri_mrs_dynamic_createswf(window, rawnumber);

swstruc = struct();

for s=1:size(swf, 2)
    swfiles = rawfiles(find(swf(:, s))); %#ok<FNDSB>
    swnumber = size(swfiles, 1);
    
    for mm=1:swnumber
        disp(['Copying METAB file: ', num2str(mm,'%04d'), ' / ', num2str(swnumber), ' step: ', num2str(s), '/', num2str(size(swf, 2))]);
        command = ['cp ', path_metab, filesep, swfiles(mm).name, ' ', outputfolder, filesep, prefixmetab, filesep, num2str(mm, '%04d'), ext];
        system(command);
    end

    % Load job file
    MRS_struc = OspreyJob([jobname, '.m']);

    % Assign relevant paths
    MRS_struc.files        = {[outputfolder, filesep, prefixmetab, filesep]};
    MRS_struc.files_ref    = {[path_ref, filesep]};
    MRS_struc.files_nii    = {[outputfolder, filesep, 'anat.nii']};
    MRS_struc.outputFolder = outputfolder;

    % Osprey pipeline
    MRS_struc = OspreyLoad(MRS_struc);
    MRS_struc = OspreyProcess(MRS_struc);
    MRS_struc = OspreyFit(MRS_struc);
    MRS_struc = OspreyCoreg(MRS_struc);
    MRS_struc = OspreySeg(MRS_struc);
    MRS_struc = OspreyQuantify(MRS_struc);
    
    swmetabdiff(s, :) = table2array(MRS_struc.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 2});
    swmetaboff(s, :) = table2array(MRS_struc.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 1});
    swmetabqual(s, 1:numel(MRS_struc.QM.tables{1, :})) = MRS_struc.QM.tables{1, :};

    swnamesrow{s} = ['Frame_', num2str(s, '%03d')];
    swnamesmetabdiff = [MRS_struc.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 2}.Properties.VariableNames];
    swnamesmetaboff = [MRS_struc.quantify.tables.metab.TissCorrWaterScaled.Voxel_1{1, 1}.Properties.VariableNames];
    swnamesqual = [MRS_struc.QM.tables.Properties.VariableNames];
    swspecsize = numel(MRS_struc.processed.metab{1,1}.ppm);
    swspecqual(s, 1:swspecsize) = real(MRS_struc.processed.metab{1,1}.specs(1,3));
    
%     save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.mat'], 'swmetabdiff');
%     save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metaboff.mat'], 'swmetaboff');
%     save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabqual.mat'], 'swmetabqual');
%     save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_specqual.mat'], 'swspecqual');
    
    swdifftable = array2table(swmetabdiff, 'VariableNames', swnamesmetabdiff, 'RowNames', swnamesrow);
    writetable(swdifftable,[outputfolder, filesep, 'Osprey_Dynamic_MRS_swmetabdiff.csv'],'WriteRowNames',true);
    
    swofftable = array2table(swmetaboff, 'VariableNames', swnamesmetaboff, 'RowNames', swnamesrow);
    writetable(swofftable,[outputfolder, filesep, 'Osprey_Dynamic_MRS_swmetaboff.csv'],'WriteRowNames',true);
    
    swqualtable = array2table(swmetabqual, 'VariableNames', swnamesqual, 'RowNames', swnamesrow);
    writetable(swqualtable,[outputfolder, filesep, 'Osprey_Dynamic_MRS_swmetabqual.csv'],'WriteRowNames',true);
    
    delete([outputfolder, filesep, prefixmetab, filesep, '*', ext]);
    delete('lnifmri_mrs_osprey_MM_job.mat');
    clear MRS_struc
end%for

swstruc.diff = swmetabdiff;
swstruc.off = swmetaboff;
swstruc.qual = swmetabqual;
swstruc.swf = swf;
swstruc.TableDiff = swdifftable;
swstruc.TableOff = swofftable;
swstruc.TableQual = swqualtable;

end%function