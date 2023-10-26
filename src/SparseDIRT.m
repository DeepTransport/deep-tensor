classdef SparseDIRT < DIRT
    
    properties (Constant = true)
        defaultSIRTOpt = SparseOption();
    end
    
    methods
        obj = build_adapt(obj, func, bases, sirt_opt)
        
        function obj = SparseDIRT(func, arg, varargin)
            obj@DIRT(func, arg, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ns = get_pre_sample_size(obj)
            ns = obj.dirt_opt.n_samples;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function deb = get_inputdata(obj, base, n_layers, samples, density)
            if obj.dirt_opt.n_debugs > 0
                ind = 1:min(obj.dirt_opt.n_samples,obj.dirt_opt.n_debugs);
                switch obj.lis_opt.method
                    case {'none', 'unitary'}
                        %deb = InputData([], samples(:,ind), density(ind));
                        tmp = SIRT.get_potential_to_density(base, density(ind), samples(:,ind));
                        deb = InputData([], samples(:,ind), tmp);
                    case {'reduction'}
                        keep_ind = 1:obj.basis_r(n_layers+1);
                        tmp = SIRT.get_potential_to_density(base, density(ind), samples(keep_ind,ind));
                        deb = InputData([], samples(keep_ind,ind), tmp);
                        %deb = InputData([], samples(keep_ind,ind));
                end
            else
                deb = InputData();
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function irt = get_new_layer(obj, func, bases, sirt_opt, n_layers, samples, density)
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
                        % debugging and initialization samples
                        deb = get_inputdata(obj, tmp_base, n_layers, samples, density);
                        irt = SparseSIRT(newf, tmp_base, sirt_opt, 'var', deb, ...
                            'defensive', obj.dirt_opt.defensive);
                    else
                        ip = length(obj.prev_approx);
                        ip = min(ip, n_layers+1);
                        % debugging and initialization samples
                        deb = get_inputdata(obj, obj.prev_approx{ip}.base, n_layers, samples, density);
                        irt = SparseSIRT(newf, obj.prev_approx{ip}, sirt_opt, ...
                            'var', deb, 'defensive', obj.dirt_opt.defensive);
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
                    % debugging and initialization samples
                    deb = get_inputdata(obj, tmp_base, n_layers, samples, density);
                    irt = SparseSIRT(newf, tmp_base, sirt_opt, 'var', deb, ...
                        'defensive', obj.dirt_opt.defensive);
            end
            irt = clean(irt);
        end
    end
end