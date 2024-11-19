% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20220504
% ########################################################################

function MRS_struc = lnifmri_mrs_static_osprey(jobname, path_metab, path_ref, number, same, outputfolder)

rng default;
prefixmetab = 'metabfile';
ext = '.ima';
% Create dir to store raw MRS data
createrawdir = ['mkdir ', outputfolder, filesep, prefixmetab];
system(createrawdir);
% Copy raw MRS data

mfiles = dir([path_metab, filesep, '*.IMA']);
mnumber = size(mfiles, 1);
for mm=1:mnumber
    disp(['Copying METAB file: ', num2str(mm,'%04d'), ' / ', num2str(mnumber)]);
    command = ['cp ', path_metab, filesep, mfiles(mm).name, ' ', outputfolder, filesep, prefixmetab, filesep, num2str(mm, '%04d'), ext];
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
% MRS_struc.opts.savePDF = 1;
MRS_struc = OspreyLoad(MRS_struc);
MRS_struc = OspreyCoreg(MRS_struc);
MRS_struc = OspreySeg(MRS_struc);
MRS_struc = OspreyProcess(MRS_struc);
MRS_struc = OspreyFit(MRS_struc);
MRS_struc = OspreyQuantify(MRS_struc);

end%function