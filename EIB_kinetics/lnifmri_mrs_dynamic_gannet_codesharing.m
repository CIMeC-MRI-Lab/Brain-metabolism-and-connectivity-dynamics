% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # REV. 20220504
% ########################################################################

function swstruc = lnifmri_mrs_dynamic_gannet(path_metab, path_ref, window, gannetpath, outputfolder)

rng default;
prefixmetab = 'filemetabdynamic';
prefixwater = 'filewaterdynamic';
ext = '.ima';
gannetoutputpath = [fileparts(gannetpath), filesep, 'Gannet_output'];

cmd1 = ['mkdir ', gannetpath, filesep, prefixmetab];
cmd2 = ['mkdir ', gannetpath, filesep, prefixwater];
system(cmd1);
system(cmd2);

% List raw MRS METAB data
rawfiles = dir([path_metab, filesep, '*.IMA']);
rawnumber = size(rawfiles, 1);

% Process WATER REFERENCE files
wfiles = dir([path_ref, filesep, '*.IMA']);
wnumber = size(wfiles, 1);
for ww=1:wnumber
    disp(['Copying WATER file: ', num2str(ww,'%04d'), ' / ', num2str(wnumber)]);
    command = ['cp ', path_ref, filesep, wfiles(ww).name, ' ', gannetpath, filesep, prefixwater, filesep, num2str(ww, '%04d'), ext];
    system(command);
end

%define a sliding window filter
swf = lnifmri_mrs_dynamic_createswf(window, rawnumber);

swstruc = struct();

for s=1:size(swf, 2)
    % Process METAB files
    swfiles = rawfiles(find(swf(:, s))); %#ok<FNDSB>
    swnumber = size(swfiles, 1);
    
    for mm=1:swnumber
        filecode = [num2str(mm,'%04d'), '_', num2str(s,'%04d')];
        disp(['Copying METAB file: ', num2str(mm,'%04d'), ' / ', num2str(swnumber), ' step: ', num2str(s), '/', num2str(size(swf, 2))]);
        %command = ['cp ', path_metab, filesep, swfiles(mm).name, ' ', gannetpath, filesep, prefixmetab, filesep, num2str(mm, '%04d'), ext];
        command = ['cp ', path_metab, filesep, swfiles(mm).name, ' ', gannetpath, filesep, prefixmetab, filesep, filecode, ext];
        system(command);
    end
    
    % Run GANNET
    MRS_struct = GannetLoad({[gannetpath, filesep, prefixmetab, filesep, '0001_', num2str(s,'%04d'), ext]}, {[gannetpath, filesep, prefixwater, filesep, '0001', ext]});
    target = cell2mat(MRS_struct.p.target);
    switch target
        case 'GABAGlx'
            % Load all the files in Gannet
            MRS_struct=GannetFit(MRS_struct);
            MRS_struct=GannetCoRegister(MRS_struct, {[outputfolder, filesep, 'anat.nii']});
            MRS_struct=GannetSegment(MRS_struct);
            MRS_struct=GannetQuantify(MRS_struct);
            save([outputfolder, filesep, 'Gannet_Dynamic_MRS_struct.mat']);
    end %switch target
    
    swmetabdiff(s, 1) = MRS_struct.out.vox1.GABA.ConcIU_TissCorr;
    swmetabdiff(s, 2) = MRS_struct.out.vox1.GABA.ConcIU_AlphaTissCorr;
    swmetabdiff(s, 3) = MRS_struct.out.vox1.Glx.ConcIU_TissCorr;
    swmetabdiff(s, 4) = MRS_struct.out.vox1.Glx.ConcIU_AlphaTissCorr;

    % UNUSED - left here for compatibility with OSPREY version
    swmetaboff(s, 1) = 0;
    swmetaboff(s, 2) = 0;
    swmetaboff(s, 3) = 0;
    swmetaboff(s, 4) = 0;
    
    swnamesrow{s} = ['Frame_', num2str(s, '%03d')];
    swnamesmetab(1:2) = {'GABA_ConcIU_TissCorr', 'GABA_ConcIU_AlphaTissCorr'};
    swnamesmetab(3:4) = {'Glx_ConcIU_TissCorr', 'Glx_ConcIU_AlphaTissCorr'};
    
    swmetabqual(s, 1) = MRS_struct.out.vox1.GABA.FWHM;
    swmetabqual(s, 2) = MRS_struct.out.vox1.GABA.SNR;
    swmetabqual(s, 3) = MRS_struct.out.vox1.GABA.FitError;
    swmetabqual(s, 4) = MRS_struct.out.vox1.Glx.FWHM;
    swmetabqual(s, 5) = MRS_struct.out.vox1.Glx.SNR;
    swmetabqual(s, 6) = MRS_struct.out.vox1.Glx.FitError;
    swmetabqual(s, 7) = MRS_struct.out.vox1.tissue.fGM;
    swmetabqual(s, 8) = MRS_struct.out.vox1.ResidWater.SuppressionFactor;
    swnamesqual(1:8) = {'GABA_FWHM', 'GABA_SNR', 'GABA_FitError', 'Glx_FWHM', 'Glx_SNR', 'Glx_FitError', 'fGM', 'WaterSuppFactor'};
    
    swspecsize(s) = numel(MRS_struct.spec.vox1.GABAGlx.diff_scaled);
    swspecqual(1:swspecsize(s), s) = real(MRS_struct.spec.vox1.GABAGlx.diff_scaled);
    specfreq = MRS_struct.spec.freq;

    % save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabdiff.mat'], 'metabdiff');
    % save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metaboff.mat'], 'metaboff');
    % save([resultspath, filesep, timestamp, '_', keyword, '_', fitpackage{fittingmethod}, '_MRS_metabqual.mat'], 'metabqual');
    % save([outputfolder, filesep, 'Gannet_Dynamic_MRS_swspecqual.mat'], 'swspecqual');
    
    swdifftable = array2table(swmetabdiff, 'VariableNames', swnamesmetab, 'RowNames', swnamesrow);
    writetable(swdifftable,[outputfolder, filesep, 'Gannet_Dynamic_MRS_swmetabdiff.csv'],'WriteRowNames',true);
    
    swofftable = array2table(swmetaboff, 'VariableNames', swnamesmetab, 'RowNames', swnamesrow);
    writetable(swofftable,[outputfolder, filesep, 'Gannet_Dynamic_MRS_swmetaboff.csv'],'WriteRowNames',true);
    
    swqualtable = array2table(swmetabqual, 'VariableNames', swnamesqual, 'RowNames', swnamesrow);
    writetable(swqualtable,[outputfolder, filesep, 'Gannet_Dynamic_MRS_swmetabqual.csv'],'WriteRowNames',true);
    
    writematrix(specfreq, [outputfolder, filesep, 'specfreq.csv'])

    % Cleanup
    delete([gannetpath, filesep, prefixmetab, filesep, '*']);
%     delete([gannetpath, filesep, prefixwater, filesep, '*']);
    movepdf = ['cp ', gannetoutputpath, filesep, '*.pdf ', outputfolder, filesep];
    system(movepdf);
    listpdf = dir([outputfolder, filesep, 'Gannet*.pdf']);
    for d=1:size(listpdf, 1)
        renamepdf = ['mv ', listpdf(d).folder, filesep, listpdf(d).name, ' ', listpdf(d).folder, filesep, filecode, '_', listpdf(d).name];
        system(renamepdf);
    end
    delete([gannetoutputpath, filesep, '*.pdf']);
    
    copygannetsetup = ['cp ', gannetpath, filesep, 'GannetPreInitialise.m ', outputfolder, filesep];
    system(copygannetsetup);

    clear MRS_struc
end%for

swstruc.diff = swmetabdiff;
swstruc.off = swmetaboff;
swstruc.qual = swmetabqual;
swstruc.specqual = swspecqual;
swstruc.specfreq = specfreq;
swstruc.swf = swf;
swstruc.TableDiff = swdifftable;
swstruc.TableOff = swofftable;
swstruc.TableQual = swqualtable;

end%function