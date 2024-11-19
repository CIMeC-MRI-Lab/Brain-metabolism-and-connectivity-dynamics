% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20220504
% ########################################################################

function MRS_struct = lnifmri_mrs_static_gannet(model, path_metab, path_ref, refname, number, same, gannetpath, outputfolder)

rng default;
prefixmetab = 'filemetab';
prefixwater = 'filewater';
ext = '.ima';
gannetoutputpath = [fileparts(gannetpath), filesep, 'Gannet_output'];

cmd1 = ['mkdir ', gannetpath, filesep, prefixmetab];
cmd2 = ['mkdir ', gannetpath, filesep, prefixwater];
system(cmd1);
system(cmd2);

% Process METAB files
mfiles = dir([path_metab, filesep, '*.IMA']);
mnumber = size(mfiles, 1);
for mm=1:mnumber
    disp(['Copying METAB file: ', num2str(mm,'%04d'), ' / ', num2str(mnumber)]);
    command = ['cp ', path_metab, filesep, mfiles(mm).name, ' ', gannetpath, filesep, prefixmetab, filesep, num2str(mm, '%04d'), ext];
    system(command);
end

% Process WATER REFERENCE files
wfiles = dir([path_ref, filesep, '*.IMA']);
wnumber = size(wfiles, 1);
for ww=1:wnumber
    disp(['Copying WATER file: ', num2str(ww,'%04d'), ' / ', num2str(wnumber)]);
    command = ['cp ', path_ref, filesep, wfiles(ww).name, ' ', gannetpath, filesep, prefixwater, filesep, num2str(ww, '%04d'), ext];
    system(command);
end

% Run GANNET
MRS_struct = GannetLoad({[gannetpath, filesep, prefixmetab, filesep, '0001', ext]}, {[gannetpath, filesep, prefixwater, filesep, '0001', ext]});
target = cell2mat(MRS_struct.p.target);
switch target
    case 'GABAGlx'
        % Load all the files in Gannet
        MRS_struct=GannetFit(MRS_struct);
        MRS_struct=GannetCoRegister(MRS_struct, {[outputfolder, filesep, 'anat.nii']});
        MRS_struct=GannetSegment(MRS_struct);
        MRS_struct=GannetQuantify(MRS_struct);
        save([outputfolder, filesep, 'Gannet_MRS_struct.mat']);
        movepdf = ['cp ', gannetoutputpath, filesep, '*.pdf ', outputfolder, filesep];
        system(movepdf);
end %switch target

% Cleanup
delete([gannetpath, filesep, prefixmetab, filesep, '*']);
delete([gannetpath, filesep, prefixwater, filesep, '*']);
delete([gannetoutputpath, filesep, '*.pdf']);

end%function