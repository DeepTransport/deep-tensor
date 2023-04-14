classdef TTDIRT < DIRT
    methods
        function obj = TTDIRT(func, arg, varargin)
            %{
            %
            p = inputParser;
            addRequired(p, 'func',  @(x) isa(x, 'function_handle'));
            addRequired(p, 'arg');
            addRequired(p, 'bridge');
            %
            addOptional(p, 'ref', DIRT.defaultRef);
            %
            addOptional(p, 'sirt_option',defaultSIRTOption);
            addOptional(p, 'dirt_option',DIRT.defaultDIRTOption);
            addOptional(p, 'lis_option',DIRT.defaultLISOption);
            %
            addParameter(p, 'init_samples',DIRT.defaultInitSamples);
            %
            p.KeepUnmatched = true;
            parse(p,func,arg,bridge,varargin{:});
            %
            ref = p.Results.ref;
            sirt_opt = p.Results.sirt_option;
            dirt_opt = p.Results.dirt_option;
            lis_opt = p.Results.lis_option;
            %
            init_samples = p.Results.init_samples;
            %
            obj@DIRT(func, arg, bridge, ref, sirt_opt, dirt_opt, lis_opt, init_samples);
            %}
            defaultSIRTOption = FTTOption('max_als', 2, 'tt_method', 'random');
            obj@DIRT(defaultSIRTOption, func, arg, varargin{:});
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ns = get_pre_sample_size(obj)
            ns = obj.dirt_opt.n_samples+obj.dirt_opt.n_debugs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function deb = get_debugger(obj, n_layers, samples, density)
            if obj.dirt_opt.n_debugs > 0
                ind = (1:obj.dirt_opt.n_debugs)+obj.dirt_opt.n_samples;
                switch obj.lis_opt.method
                    case {'none', 'unitary'}
                        deb = Debugger(samples(:,ind), density(ind));
                    case {'reduction'}
                        keep_ind = 1:obj.basis_r(n_layers+1);
                        deb = Debugger(samples(keep_ind,ind), density(ind));
                end
            else
                deb = Debugger();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function irt = get_new_layer(obj, func, bases, sirt_opt, n_layers, samples, density)
            ref_opt = sirt_opt;
            if n_layers <= 1
                ref_opt.max_als = max(2,sirt_opt.max_als); % use at least 2 ALS iteration at the start
            end
            % debugging
            deb = get_debugger(obj, n_layers, samples, density);
            % density
            switch obj.lis_opt.method
                case {'none'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, []);
                case {'reduction', 'unitary'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, obj.basis{n_layers+1});
            end
            % bases and samples
            ind = 1:obj.dirt_opt.n_samples;
            switch obj.lis_opt.method
                case {'none', 'unitary'}
                    %keep_ind = 1:obj.d;
                    if n_layers > 1
                        tmp_base = obj.irts{n_layers}.approx;
                        tmp_base.opt = ref_opt;
                    else
                        tmp_base = bases{n_layers+1};
                    end
                    % approximate
                    if isempty(obj.prev_approx)
                        irt = TTSIRT(newf, tmp_base, ref_opt, ...
                            'defensive', obj.dirt_opt.defensive, ...
                            'samples', samples(:,ind), 'debug', deb);
                    else
                        ip = length(obj.prev_approx);
                        ip = min(ip, n_layers+1);
                        irt = TTSIRT(newf, tmp_base, ref_opt, ...
                            'defensive', obj.dirt_opt.defensive, ...
                            'samples', samples(:,ind), 'debug', deb, ...
                            'previous', obj.prev_approx{ip});
                    end
                case {'reduction'}
                    % change the dimension of the approximation base
                    %keep_ind = 1:obj.basis_r(n_layers+1);
                    remo_ind = (obj.basis_r(n_layers+1)+1):obj.d;
                    if n_layers == 0
                        tmp_base = remove_bases(bases{1}, remo_ind);
                    else
                        tmp_base = remove_bases(bases{end}, remo_ind);
                    end
                    % approximate
                    irt = TTSIRT(newf, tmp_base, ref_opt, ...
                        'defensive', obj.dirt_opt.defensive, ...
                        'samples', samples(:,ind), 'debug', deb);
            end
        end
    end
end
