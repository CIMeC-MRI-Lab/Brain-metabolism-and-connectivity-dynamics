% ########################################################################
% #
% # LNiF MRI Lab Pipeline for preprocessing of
% # Spectroscopy MEGA-PRESS GABA-edited dataset
% # stefano.tambalo@unitn.it - 20211126
% #
% # Screate a sliding window filter to select off-on spectra for dynamic
% # spectral quantification of metabolites
% #
% #
% # Input:
% #
% #    window = two-element vector containing window length and sliding step size
% #    mnumber = number of available ON-OFF pairs of MEGA-PRESS spectra
% #    
% # Output:
% #    
% #    swf = sliding window filter, a binary matrix whose columns can be used
% #            to select corresponding ON-OFF pairs of MEGA-PRESS spectra
% #
% # REV. 20220513
% ########################################################################

function swf = lnifmri_mrs_dynamic_createswf(window, mnumber)
    wl = window(1);
    ws = window(2);
    init = zeros(mnumber/2, 1);
    init(1:wl) = 1;
    nshift = 1;
    wf(:, nshift) = init;
    rem = numel(init)-wl;
    while rem > ws
        nshift = nshift+1;
        wf(:, nshift) = circshift(wf(:, nshift-1), ws);
        rem = numel(init) - nnz(sum(wf, 2));
    end%while
    swf = repmat(wf, 2, 1);
end%function