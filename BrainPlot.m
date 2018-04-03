classdef BrainPlot < handle
    % Creates a plot of an fMRI dataset à la J.D. Power's "The Plot",
    % https://doi.org/10.1016/j.neuroimage.2016.08.009
    %
    %
    % To use, instantiate an object of the class by passing in a scalar
    % `opts` structure with plot parameters, and call the public method
    % `make`. An empty opts structure can be obtained by calling the static
    % method `defaults`, and modifying the structure suitably. E.g.,
    %
    %     % Assume that in the current working directory, we have three
    %     % masks: gray_mask.nii, white_mask.nii, and
    %     % csf_mask.nii; a realigment parameters file rp.txt; and our
    %     % SPM.mat.
    %
    %     % Set up plot parameters
    %     load SPM.mat
    %     rp = load('rp.txt');
    %     opts = BrainPlot.defaults(); % Get the defaults structure.
    %     opts.brain.vols = spm_vol(SPM.xY.VY);
    %     opts.mask.vols = spm_vol(spm_select('FPList', pwd, '^.*_mask.nii$'));
    %     opts.mask.labels = cellstr(spm_select('List', pwd, '^.*_mask.nii$'));
    %     opts.mask.labels = cellfun(@(s)strrep(s, '_', ' '), opts.mask.labels, ...
    %                            'UniformOutput',false);
    %     opts.filter.filter = 1; % We will apply the high pass filter
    %     opts.whiten.whiten = 1; % We will apply the whitening matrix
    %     opts.adjust.adjust = 1; % We will adjust for confounds (like
    %                             % motion parameters
    %     opts.adjust.contrast = 5;  % We will use an F-contrast indexed 5
    %                                % to adjust.
    %     opts.spmpath = '/path/to/SPM.mat';
    %     opts.extra(1).data = rp(:,1:3);
    %     opts.extra(1).ylabel = 'translation rp';
    %     opts.extra(2).data = rp(:,4:6);
    %     opts.extra(2).ylabel = 'rotation rp';
    %
    %     % Instantiate the object, make the plot, and save.
    %     bp = BrainPlot(opts);
    %     bp.make(); % May take a while.
    %     bp.save(fullfile(pwd, '1234.jpg'));
    %
    
    
    % Static methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static = true)
        function d = defaults()
            % Static method to generate a default option structure.
            % -----------------------------------------------------
            % The scalar opts structure is passed to the constructor, and
            % has the following fields.
            %
            % brain_vols: fMRI volumes of interest; the data structure
            %             returned by `spm_vol`.
            %
            % mask:       mask = struct('vols',[],'labels',[])
            %             and is scalar.
            % mask.vols:  A structure of nifti volumes, as returned from
            %             spm_vol, e.g.
            %             opts.mask.vols = spm_vol(char({'gray_mask.nii'
            %                                            'white_mask.nii'
            %                                            'csf_mask.nii'}))
            % mask.labels A cell array of labels for the masks, e.g.
            %             opts.mask.labels = {'Gray', 'White', 'CSF'};
            %
            % spmpath:    Leave as [] if not using. Path to your SPM.mat file.
            %             The SPM.mat file is used to adjust data, filter
            %             the data with the SPM's specification, and to
            %             obtain the SPM's whitening matrix.
            %
            % mask_threshold: Required. The exact value will depend on the
            %             type of mask. In most cases, just leave as is.
            %
            % adjust:     adjust = struct('adjust', false, contrast, [])
            %             Whether to adjust to a F-contrast of interest
            %             from the SPM. If not adjusting, the data will be
            %             detrended and demeaned.
            % adjust.adjust: truthy value if adjusting for contrast.
            % adjust.contrast: index of F-contrast to adjust for. Columns
            %             of the F-contrast which are all zeros will be
            %             treated as confounds and removed. If left empty,
            %             then all effects will be modeled out, and the
            %             residuals of the model returned.
            %
            % whiten:     whiten = struct('whiten', false, 'W', [])
            %             Whether to whiten the data before plotting.
            % whiten.whiten: truthy value if wanting to whiten data.
            % whiten.W:   The whitening matrix to whiten data with. If left
            %             as [], then W = SPM.xX.W
            %
            % filter:     filter = struct('filter', false, 'K', [])
            %             Whether to apply a high pass filter to the data.
            % filter.filter: truthy value if wanting to filter before
            %             plotting.
            % filter.K    Parameters to specify the high pass filter, as
            %             used by spm_filter. If left as [], then SPM.xX.K
            %             is used.
            %
            % extra:      extra = struct('data',{}, 'ylabel',{})
            %             A possibly zero-length structure array specifying
            %             additional data to plot, e.g. frame-displacement,
            %             or the realignment parameters. For N extra plots,
            %             the structure should have length N.
            % extra.data: nxm matrix, with n records, for m features.
            % extra.ylabel: ylabel for item
            %
            % display:    A structure containing parameters for appearance
            %             of plot, such as background color.
            % display.figure: Is passed to figure(), and contains figure
            %            Properties, such as Name, NumberTitle, etc.
            %            Crucially, it turns off InvertHardCopy, so that
            %            saving the figure saves it to file as it appears
            %            on screen.
            % display.text.Color: white by default, since background is
            %            black by default.
            % display.text.Rotation: 45 by default
            % display.text.FontSize: 8 by deafult
            % display.scaling: if not empty, specifies a CLIM for the image
            %            map produced by imagesc.
            % display.colormap: colormap for brain plot. 'gray' by default.
            %
            %
            % save:      save = structure('savefile',[], 'format',[])
            %            If not empty, will automatically save file using
            %            path `savefile`, optionally with format `format`.
            %
            
            d = struct();
            d.brain.vols = [];
            d.mask.vols = [];
            d.mask.labels = [];
            d.spmpath = [];
            d.mask_threshold = 0.95;
            d.adjust.adjust = false;
            d.adjust.contrast = [];
            d.filter.filter = false;
            d.filter.K = [];
            d.whiten.whiten = false;
            d.whiten.W = [];
            % the {} forces a zero length structure:
            d.extra = struct('data', {}, 'ylabel', {});
            d.display.figure.Name = 'Compartment Plot';
            d.display.figure.NumberTitle = 'on';
            d.display.figure.InvertHardCopy = 'off';
            d.display.figure.Color = 'black';
            d.display.text.Color = 'white';
            d.display.text.Rotation = 45;
            d.display.text.FontSize = 8;
            d.display.scaling = [];
            d.display.colormap = 'gray';
            d.save.savefile = [];
            d.save.format = [];
        end
        
        function xyz = get_xyz(dim)
            % Static method to construct 3xN matrix of voxel coordinates.
            % -----------------------------------------------------------
            % Args:
            %    dim: 1X3 matrix of dimensions of source image (voxels)
            % Returns:
            %    xyz: 3xN matrix of voxel coordinates.
            xyz = zeros(3, dim(1)*dim(2)*dim(3));
            index = 1;
            for z = 1:dim(3)
                for y = 1:dim(2),
                    for x = 1:dim(1)
                        xyz(1:3, index) = [x;y;z];
                        index=index+1;
                    end;
                end;
            end
        end
        
        function prop_copy(source, target_handle)
            % Static method to copy properties from one figure to another.
            % ------------------------------------------------------------
            % source_handle: handle of the source figure, or its
            %                properties as retrieved by get(gcf).
            % target_handle: handle of the target figure to copy into.
            if isstruct(source)
                children = source.Children;
            else
                children = get(source, 'children');
            end
            if ~isempty(children)
                h = copyobj(children, target_handle);
                for ii = 1 : numel(children)
                    prop_copy(children(ii), h(ii));
                end
            end
            return
        end
    end
    
    % Public properties
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        compartments = struct();
        SPM = [];
        beta = [];
        opts = [];
        XYZ = [];
        fig = [];
        fig_properties = [];
    end
    
    % Public methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = BrainPlot(opts)
            % Constructor
            % opts - structure containing parameters for plot. Static
            %        method BrainPlot.defaults() produces a default opts
            %        structure that can be used.
            obj.opts = opts;
        end
        
        function make(obj)
            % Call this method to load the brain data and construct the
            % plot.
            obj.parse_opts();
            obj.XYZ = BrainPlot.get_xyz(obj.opts.brain.vols(1).dim);
            obj.make_compartments();
            obj.prepare_data();
            obj.create_image();
            if ~isempty(obj.opts.save.savefile)
                obj.save();
            end
        end
        
        function save(obj, savefile, format)
            % Method to save plot to file. A relatively thin wrapper around
            % `saveas`. If plot was closed, calling this method will
            % reconstruct the plot by recalling the figure properties
            % generated by make.
            %
            % Args:
            %   savefile: filepath to save plot. If not specified, then
            %             looks for filename in opts.save.savefile.
            %
            %   format:   format argument, as used by `saveas`.
            %             If not specified, uses opts.save.format.
            %
            if ~exist('savefile', 'var')
                savefile = obj.opts.save.savefile;
            end
            if ~exist('format', 'var')
                format = obj.opts.save.format;
            end
            if isempty(savefile)
                error('Empty save file name.');
            end
            
            % Get the figure
            if isempty(obj.fig)
                % Oops, no figure
                warning('Calling save before figure made. Call make()');
                return
            end
            try
                v = get(obj.fig);
                f = figure(obj.fig);
            catch
                % Figure unavailable. But we have the
                % properties, so we'll create a new one from those.
                f = figure();
                fig_props = obj.fig_properties;
                BrainPlot.prop_copy(fig_props, f);
            end
            
            % Save it
            if ~isempty(format)
                saveas(f, savefile, format);
            else
                saveas(f, savefile);
            end
            
        end
    end
    
    
    
    % Private methods
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access = 'private')
        
        function parse_opts(obj)
            defaults = BrainPlot.defaults();
            try
                flds = fieldnames(obj.opts);
            catch
                ME = MException('InvalidInput:NotStructure', ...
                    'Arg `opts` is not a structure.');
                throw(ME);
            end
            for ii = 1 : length(flds)
                defaults.(flds{ii}) = obj.opts.(flds{ii});
            end
            
            
            % Make sure SPM is structure is loaded, if given.
            if ~isempty(obj.opts.spmpath)
                obj.SPM = load(obj.opts.spmpath);
                obj.SPM = obj.SPM.SPM;
            end
            
            % Make sure image vols image vols
            if ~isempty(obj.opts.brain.vols)
                obj.opts.brain.vols = spm_vol(obj.opts.brain.vols);
            else
                try
                    obj.opts.brain.vols = obj.SPM.xY.VY;
                catch
                    error(['No brain data: if no image volumes specified, ' ...
                        'then you must provide a SPM structure.']);
                end
            end
            if ~isempty(obj.opts.brain.vols)
                obj.opts.mask.vols = spm_vol(obj.opts.mask.vols);
            end
        end
        
        function make_compartments(obj)
            % Masks and vectorizes the data.
            num_mask_vols = numel(obj.opts.mask.vols);
            for ii = 1:num_mask_vols
                mask_vol = obj.opts.mask.vols(ii);
                trans = obj.opts.brain.vols(1).mat;
                augxyz = [obj.XYZ; ones(1, size(obj.XYZ,2))];
                j = mask_vol.mat \ trans * augxyz;
                mask = spm_get_data(mask_vol, ...
                    j);
                Q = mask ~= 0 & ~isnan(mask);
                xyz = obj.XYZ(:,Q);
                obj.compartments(ii).data = spm_get_data(...
                    obj.opts.brain.vols, xyz);
                obj.compartments(ii).XYZ = xyz;
                obj.compartments(ii).Q = Q;
            end
        end
        
        
        function prepare_data(obj)
            % Adjusts data according to options set in opts.
            for ii = 1 : numel(obj.compartments)
                
                % Whiten
                if obj.opts.whiten.whiten
                    if isempty(obj.opts.whiten.W)
                        W = obj.SPM.xX.W;
                    else
                        W = obj.opts.whiten.W;
                    end
                    whitened = W * obj.compartments(ii).data;
                    obj.compartments(ii).data = whitened;
                end
                
                % High Pass Filter
                if obj.opts.filter.filter
                    if isempty(obj.opts.filter.K)
                        K = obj.SPM.xX.K;
                    else
                        K = obj.opts.filter.K;
                    end
                    obj.compartments(ii).data = spm_filter(K, ...
                        obj.compartments(ii).data);
                end
                
                % Adjust
                if obj.opts.adjust.adjust
                    [spmdir, ~, ~] = fileparts(obj.opts.spmpath);
                    for jj = 1 : numel(obj.SPM.Vbeta)
                        [~, fname, ext] = fileparts(...
                            obj.SPM.Vbeta(jj).fname);
                        obj.SPM.Vbeta(jj).fname = fullfile(spmdir, [fname ext]);
                    end
                    obj.beta = spm_get_data(obj.SPM.Vbeta, ...
                        obj.compartments(ii).XYZ);
                    if isempty(obj.opts.adjust.contrast)
                        % Adjust for everything, i.e. get residuals
                        obj.compartments(ii).data = obj.compartments(ii).data - ...
                            obj.SPM.xX.xKXs.X * obj.beta;
                    else
                        predicted = spm_FcUtil('Y0', ...
                            obj.SPM.xCon(obj.opts.adjust.contrast), ...
                            obj.SPM.xX.xKXs, obj.beta);
                        obj.compartments(ii).data = obj.compartments(ii).data - ...
                            predicted;
                    end
                end
                
                % Detrend and demean, but only if not adjusted already
                if ~obj.opts.adjust.adjust
                    % Regression using left matrix division.
                    r0 = ones(1, size(obj.compartments(ii).data, 1));
                    r1 = linspace(0, 1, size(obj.compartments(ii).data, 1));
                    r = [r0; r1];
                    betas = r' \ obj.compartments(ii).data;
                    prediction = r' * betas;
                    obj.compartments(ii).data = obj.compartments(ii).data - ...
                        prediction;
                end
                
                % Now transpose because that's how we're going to plot it
                obj.compartments(ii).data = obj.compartments(ii).data';
            end
        end
        
        function obj = create_image(obj)
            % Creates the figure;
            
            f = figure(obj.opts.display.figure);
            hold on;
            num_subplots = numel(obj.opts.extra) + numel(obj.compartments);
            
            % Plot `extra` data
            % -----------------
            current = 0;
            for ii = 1 : numel(obj.opts.extra)
                current = current + 1;
                p(current) = subplot(num_subplots, 1, ii);
                plot(obj.opts.extra(ii).data);
                ax = ylabel(obj.opts.extra(ii).ylabel);
                flds = fields(obj.opts.display.text);
                for jj = 1 : numel(flds)
                    set(ax, flds{jj}, obj.opts.display.text.(flds{jj}));
                end
                set(p(current), 'ytick',[])
                set(p(current), 'yticklabel',[])
            end
            
            % Plot brain voxels
            % -----------------
            numrows = nan(numel(obj.compartments),1);
            for ii = 1 : numel(obj.compartments)
                numrows(ii) = size(obj.compartments(ii).data, 1);
            end
            total_numrows = sum(numrows);
            if ~isempty(obj.opts.display.scaling)
                scaling = obj.opts.display.scaling;
            else
                scaling = obj.get_scale();
            end
            colormap(obj.opts.display.colormap);
            for ii = 1 : numel(obj.compartments)
                current = current + 1;
                p(current) = subplot(num_subplots, 1, current);
                imagesc(obj.compartments(ii).data, scaling);
                set(p(current), 'xtick',[])
                set(p(current), 'ytick',[])
                set(p(current), 'yticklabel',[])
                set(p(current), 'XColor', 'black')
                set(p(current), 'YColor', 'black')
                set(p(current), 'LineWidth', 5)
                ax = ylabel(obj.opts.mask.labels{ii});
                flds = fields(obj.opts.display.text);
                for jj = 1 : numel(flds)
                    set(ax, flds{jj}, obj.opts.display.text.(flds{jj}))
                end
            end
            
            % Scale plots to better use space
            % -------------------------------
            vertical_margin = 0.025;
            horizontal_margin = 0.075;
            extra_top = 1 - vertical_margin;
            extra_bottom = 0.85;
            left = 0 + horizontal_margin;
            width = 1 - 2 * horizontal_margin;
            brain_top = extra_bottom - vertical_margin;
            brain_bottom = vertical_margin;
            brain_height = brain_top - brain_bottom;
            
            % Scale extra plots to fit top part of page
            bottom_pos = linspace(extra_top, extra_bottom, ...
                numel(obj.opts.extra)+1);
            for ii = 1 : numel(obj.opts.extra)
                set(p(ii), 'Position', [left, bottom_pos(ii+1), ...
                    width, bottom_pos(ii)-bottom_pos(ii+1)]);
            end
            
            % Scale each brain subplot, but ensure minimum height
            min_height = 0.05;
            height = max((numrows / total_numrows) * brain_height, min_height);
            height = (height ./ sum(height)) * brain_height;
            bottom = brain_top - cumsum(height);
            for ii = 1:numel(obj.compartments)
                set(p(ii+numel(obj.opts.extra)), ...
                    'Position', [left, bottom(ii), width, height(ii)]);
            end
            hold off;
            obj.fig = f;
            obj.fig_properties = get(obj.fig);
        end
        
        function scaling = get_scale(obj)
            data = [];
            for ii = 1 : numel(obj.compartments)
                data = cat(1, data, obj.compartments(ii).data(:));
            end
            scaling = [mean(data) - std(data), mean(data) + std(data)];
        end
    end
end