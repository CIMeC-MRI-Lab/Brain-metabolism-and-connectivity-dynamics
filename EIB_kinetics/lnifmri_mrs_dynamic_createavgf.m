% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # create a filter to select off-on spectra for dynamic
% # spectral quantification of metabolites
% #
% #
% # Input:
% #
% #    mnumber = number of available ON-OFF pairs of MEGA-PRESS spectra
% #    
% # Output:
% #    
% #    avgf = average filter, a binary matrix whose columns can be used
% #            to select corresponding ON-OFF pairs of MEGA-PRESS spectra
% #
% # REV. 20220513
% ########################################################################

function avgf = lnifmri_mrs_dynamic_createavgf(mnumber)
    step = zeros(mnumber/2, 1);
    af = [];
    for i = 1:numel(step)
        step(1:i) = 1;
        af = cat(2, af, step);
    end
    avgf = repmat(af, 2, 1);
end%function