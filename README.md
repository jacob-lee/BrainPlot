# Brain Plot


Creates a plot of an fMRI dataset à la J.D. Power's "The Plot", 
 https://doi.org/10.1016/j.neuroimage.2016.08.009

# Use

To use, instantiate an object of the class by passing in a scalar
`opts` structure with plot parameters, and call the public method
`make`. An empty opts structure can be obtained by calling the static
method `defaults`, and modifying the structure suitably. E.g.,

Assume that in the current working directory, we have three 
masks: gray_mask.nii, white_mask.nii, and 
csf_mask.nii; a realigment parameters file rp.txt; and our
SPM.mat.

```matlab
    % Set up plot parameters
    load SPM.mat
    rp = load('rp.txt');
    opts = BrainPlot.defaults(); % Get the defaults structure.
    opts.brain.vols = spm_vol(SPM.xY.VY);
    opts.mask.vols = spm_vol(spm_select('FPList', pwd, '^.*_mask.nii$'));
    opts.mask.labels = cellstr(spm_select('List', pwd, '^.*_mask.nii$'));
    opts.mask.labels = cellfun(@(s)strrep(s, '_', ' '), opts.mask.labels, ...
                            'UniformOutput',false);
    opts.filter.filter = 1; % We will apply the high pass filter
    opts.whiten.whiten = 1; % We will apply the whitening matrix
    opts.adjust.adjust = 1; % We will adjust for confounds (like
                          motion parameters
    opts.adjust.contrast = 5;  % We will use an F-contrast indexed 5
                             to adjust.
    opts.spmpath = '/path/to/SPM.mat';
    opts.extra(1).data = rp(:,1:3);
    opts.extra(1).xlabel = [];
    opts.extra(1).ylabel = 'translation rp';
    opts.extra(2).data = rp(:,4:6);
    opts.extra(2).xlabel = [];
    opts.extra(2).ylabel = 'rotation rp';

    % Instantiate the object, make the plot, and save.
    bp = BrainPlot(opts);
    bp.make(); % May take a while.
    bp.save(fullfile(pwd, '1234.jpg'));
```

For more information about the opts structure, call: `help
BrainPlot.defaults`, or look at the code.