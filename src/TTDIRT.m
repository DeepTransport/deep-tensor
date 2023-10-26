classdef TTDIRT < DIRT
    
    properties (Constant = true)
        defaultSIRTOpt = TTOption('max_als', 1, 'tt_method', 'random');
    end
    
    methods
        function obj = TTDIRT(func, arg, varargin)
            obj@DIRT(func, arg, varargin{:});
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ns = get_pre_sample_size(obj)
            ns = obj.dirt_opt.n_samples+obj.dirt_opt.n_debugs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function deb = get_inputdata(obj, base, n_layers, samples, density)
            ind1 = 1:obj.dirt_opt.n_samples;
            if obj.dirt_opt.n_debugs > 0
                ind2 = (1:obj.dirt_opt.n_debugs)+obj.dirt_opt.n_samples;
                switch obj.lis_opt.method
                    case {'none', 'unitary'}
                        tmp = SIRT.get_potential_to_density(base, density(ind2), samples(:,ind2));
                        deb = InputData(samples(:,ind1), samples(:,ind2), tmp);
                    case {'reduction'}
                        keep_ind = 1:obj.basis_r(n_layers+1);
                        tmp = SIRT.get_potential_to_density(base, density(ind2), samples(keep_ind,ind2));
                        deb = InputData(samples(keep_ind,ind1), samples(keep_ind,ind2), tmp);
                        %deb = InputData(samples(keep_ind,ind1), samples(keep_ind,ind2));
                end
            else
                deb = InputData(samples(:,ind1));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function irt = get_new_layer(obj, func, bases, sirt_opt, n_layers, samples, density)
            ref_opt = sirt_opt;
            %if n_layers <= 1
            %    ref_opt.max_als = max(2,sirt_opt.max_als); % use at least 2 ALS iteration at the start
            %end
            % density
            switch obj.lis_opt.method
                case {'none'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, []);
                case {'reduction', 'unitary'}
                    newf = @(z) ratio_fun(obj.bridge, func, obj, z, obj.basis{n_layers+1});
            end
            switch obj.lis_opt.method
                case {'none', 'unitary'}
                    if isempty(obj.prev_approx)
                        if n_layers <= 1 % start from fresh
                            % debugging and initialization samples
                            deb = get_inputdata(obj, bases{n_layers+1}, n_layers, samples, density);
                            irt = TTSIRT(newf, bases{n_layers+1}, ref_opt, 'var', deb, ...
                                'defensive', obj.dirt_opt.defensive);
                        else % start from the previous approximation
                            % debugging and initialization samples
                            deb = get_inputdata(obj, obj.irts{n_layers}.approx.base, n_layers, samples, density);
                            irt = TTSIRT(newf, obj.irts{n_layers}.approx, ref_opt, 'var', deb, ...
                                'defensive', obj.dirt_opt.defensive);
                        end
                    else
                        ip = length(obj.prev_approx);
                        ip = min(ip, n_layers+1);
                        % debugging and initialization samples
                        deb = get_inputdata(obj, obj.prev_approx{ip}.base, n_layers, samples, density);
                        irt = TTSIRT(newf, obj.prev_approx{ip}, ref_opt, 'var', deb, ...
                            'defensive', obj.dirt_opt.defensive);
                    end
                    %{
                    %keep_ind = 1:obj.d;
                    if n_layers > 1
                        tmp_base = obj.irts{n_layers}.approx;
                        tmp_base.opt = ref_opt;
                    else
                        tmp_base = bases{n_layers+1};
                    end
                    % approximate
                    if isempty(obj.prev_approx)
                        irt = TTSIRT(newf, tmp_base, ref_opt, 'var', deb, ...
                            'defensive', obj.dirt_opt.defensive);
                    else
                        ip = length(obj.prev_approx);
                        ip = min(ip, n_layers+1);
                        irt = TTSIRT(newf, tmp_base, ref_opt, 'var', deb, ...
                            'defensive', obj.dirt_opt.defensive, ...
                            'previous', obj.prev_approx{ip});
                    end
                    %}
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
                    deb = get_inputdata(obj, tmp_base, n_layers, samples, density);
                    irt = TTSIRT(newf, tmp_base, ref_opt, 'var', deb, ...
                        'defensive', obj.dirt_opt.defensive);
            end
        end
    end
end
