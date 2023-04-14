classdef SparseDIRT < DIRT
    
    methods
        obj = build_adapt(obj, func, bases, sirt_opt)
        
        function obj = SparseDIRT(func, arg, varargin)
            %{
            defaultSIRTOption = AdaptiveSparseTensorAlgorithm();
            defaultSIRTOption.tol = 5e-2;
            defaultSIRTOption.bulkParameter = 0.5;
            defaultSIRTOption.adaptiveSampling = true;
            defaultSIRTOption.adaptationRule = 'reducedmargin';
            defaultSIRTOption.display = true;
            defaultSIRTOption.displayIterations = true;
            defaultSIRTOption.nbSamples = 10;
            defaultSIRTOption.maxNbSamples = 1e3;
            defaultSIRTOption.maxDimBasis = 100;
            %}
            %
            %
            %{
            defaultSIRTOption   = SparseOption();
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
            init_samples = p.Results.init_samples;
            %}

            defaultSIRTOption = SparseOption();
            obj@DIRT(defaultSIRTOption, func, arg, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ns = get_pre_sample_size(obj)
            ns = obj.dirt_opt.n_samples;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function deb = get_debugger(obj, n_layers, samples, density)
            if obj.dirt_opt.n_debugs > 0
                ind = 1:min(obj.dirt_opt.n_samples,obj.dirt_opt.n_debugs);
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
            % debugging
            deb = get_debugger(obj, n_layers, samples, density);
            % target function
            switch obj.lis_opt.method
                case {'none'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, []);
                case {'reduction', 'unitary'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, obj.basis{n_layers+1});
            end
            % bases
            switch obj.lis_opt.method
                case {'none', 'unitary'}
                    if n_layers == 0
                        tmp_base = bases{1};
                    else
                        tmp_base = bases{end};
                    end
                    % approximate
                    if isempty(obj.prev_approx)
                        irt = SparseSIRT(newf, tmp_base, sirt_opt, ...
                            'defensive', obj.dirt_opt.defensive, 'debug', deb);
                    else
                        ip = length(obj.prev_approx);
                        ip = min(ip, n_layers+1);
                        irt = SparseSIRT(newf, tmp_base, sirt_opt, ...
                            'defensive', obj.dirt_opt.defensive, 'debug', deb, ...
                            'previous', obj.prev_approx{ip});
                    end
                case {'reduction'}
                    % change the dimension of the approximation base
                    remo_ind = (obj.basis_r(n_layers+1)+1):obj.d;
                    if n_layers == 0
                        tmp_base = remove_bases(bases{1}, remo_ind);
                    else
                        tmp_base = remove_bases(bases{end}, remo_ind);
                    end
                    % approximate
                    irt = SparseSIRT(newf, tmp_base, sirt_opt, ...
                        'defensive', obj.dirt_opt.defensive, 'debug', deb);
            end
            irt = clean(irt);
        end
    end
end