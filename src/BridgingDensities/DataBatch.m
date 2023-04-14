classdef DataBatch < SingleBeta

    properties
        masks
        masks_inc
        %
        data_order
        %
        data_layers
        %
        beta_p
    end

    methods

        function obj = DataBatch(varargin)
            %
            obj@SingleBeta(varargin{:});
            %
            obj.is_adapt = true;
            obj.data_layers = 0;
            obj.n_layers = 0;
            obj.beta_p = 0;
            obj.betas = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [str,log_weights] = print_str(obj, mllkds_vec, mlps, mlogfx)
            if obj.data_layers > 1
                ind1 = obj.data_order(obj.masks{obj.data_layers-1});
                ind2 = obj.data_order(obj.masks_inc{obj.data_layers});
                mllkds1 = sum(mllkds_vec(ind1,:),1);
                mllkds2 = sum(mllkds_vec(ind2,:),1);
                [~,dh2,~] = f_divergence(-mlogfx, -mllkds1-obj.betas(obj.n_layers,1)*mllkds2-mlps);
                log_weights = -mllkds1-obj.betas(obj.n_layers,2)*mllkds2-mlps+mlogfx;
            else
                ind = obj.data_order(obj.masks{obj.data_layers});
                mllkds = sum(mllkds_vec(ind,:),1);
                [~,dh2,~] = f_divergence(-mlogfx, -obj.betas(obj.n_layers,1)*mllkds-mlps);
                log_weights = -obj.betas(obj.n_layers,2)*mllkds-mlps+mlogfx;
            end
            log_weights = log_weights - max(log_weights);
            ess = ess_ratio(log_weights);
            %
            str = sprintf('\n\nIter=%2d, Ndata=%2d, DHell=%3.3e, beta_p=%3.3e, beta=%3.3e, ess=%3.3e', ...
                obj.n_layers, length(obj.masks_inc{obj.data_layers}), sqrt(dh2), ...
                obj.betas(obj.n_layers,1), obj.betas(obj.n_layers,2), ess);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = adapt_temperature(obj, method, mllkds_vec, mlps, mlogfx)
            beta = max(obj.beta_p,obj.min_beta);
            switch method
                case {'Aratio'}
                    ind = obj.data_order(obj.masks_inc{obj.data_layers});
                    mllkds = sum(mllkds_vec(ind,:),1);
                    % compute ess over sample size
                    ess = ess_ratio((obj.beta_p-beta)*mllkds);
                    while ess > obj.ess_tol
                        beta = beta*obj.beta_factor;
                        ess = ess_ratio((obj.beta_p-beta)*mllkds);
                    end
                case {'Eratio'}
                    if obj.data_layers > 1
                        ind1 = obj.data_order(obj.masks{obj.data_layers-1});
                        ind2 = obj.data_order(obj.masks_inc{obj.data_layers});
                        mllkds1 = sum(mllkds_vec(ind1,:),1);
                        mllkds2 = sum(mllkds_vec(ind2,:),1);
                        ess = ess_ratio(-mllkds1-beta*mllkds2-mlps+mlogfx);
                        while ess > obj.ess_tol
                            beta = beta*obj.beta_factor;
                            ess = ess_ratio(-mllkds1-beta*mllkds2-mlps+mlogfx);
                        end
                    else
                        ind = obj.data_order(obj.masks{obj.data_layers});
                        mllkds = sum(mllkds_vec(ind,:),1);
                        % compute ess over sample size
                        ess = ess_ratio(-beta*mllkds-mlps+mlogfx);
                        while ess > obj.ess_tol
                            beta = beta*obj.beta_factor;
                            ess = ess_ratio(-beta*mllkds-mlps+mlogfx);
                        end
                    end
            end
            %beta = min(1, beta/obj.beta_factor);
            beta = min(1, beta);
            if beta > 0.8 && (beta-obj.beta_p) > (1-beta)
                beta = 1;
            end
            obj.betas = [obj.betas; obj.beta_p, beta];
            obj.n_layers = obj.n_layers + 1;
            obj.beta_p = beta;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = get_ratio_fun(obj, dirt, z, mllkds_vec, mlps, ft) 
            beta_p = obj.betas(obj.n_layers,1);
            beta = obj.betas(obj.n_layers,2);
            if obj.n_layers > 1 % not in the first layer
                % compute the reference density at z
                logfz = log_joint_pdf(dirt.ref, z);
                %
                % get current mask, previous mask, and the incremental mask
                switch dirt.dirt_opt.method
                    case {'Aratio'}
                        ind = obj.data_order(obj.masks_inc{obj.data_layers});
                        mllkds = sum(mllkds_vec(ind,:),1);
                        f =  mllkds*(beta-beta_p) - logfz;
                    case {'Eratio'}
                        if obj.data_layers > 1
                            ind1 = obj.data_order(obj.masks{obj.data_layers-1});
                            ind2 = obj.data_order(obj.masks_inc{obj.data_layers});
                            mllkds1 = sum(mllkds_vec(ind1,:),1);
                            mllkds2 = sum(mllkds_vec(ind2,:),1);
                            f =  mllkds1 + mllkds2*beta + mlps - ft - logfz;
                        else
                            ind2 = obj.data_order(obj.masks{obj.data_layers});
                            mllkds2 = sum(mllkds_vec(ind2,:),1);
                            f =  mllkds2*beta + mlps - ft - logfz;
                        end
                end
            else
                ind = obj.data_order(obj.masks{1});
                mllkds = sum(mllkds_vec(ind,:),1);
                f = beta*mllkds + mlps;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function f = ratio_fun(obj, func, dirt, z, V)
            % Evaluate the ratio function during the DIRT construction
            %   f = RATIO_FUN(dirt, func, k, z, V)
            %
            %   func - function handle of the target function, returns minus log
            %          likelihood and minus log prior
            %   beta - the next beta value
            %   z    - reference random variables, d x n
            %   V    - dimension of the reduced basis or reparametrization basis
            %   f    - minus log pdf of the ratio function
            
            if nargin == 5 && ~isempty(V)
                [x, ft] = eval_irt(dirt, V*z);
            else
                [x, ft] = eval_irt(dirt, z);
            end
            [mllkds_vec, mlps] = func(x, -1);
            f = get_ratio_fun(obj, dirt, z, mllkds_vec, mlps, ft);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = islast(obj)
            islastdata = false;
            islastbeta = false;
            if obj.data_layers > 0
                islastdata = obj.masks{obj.data_layers}(end) == numel(obj.data_order);
            end
            if obj.n_layers > 0
                islastbeta = abs(obj.betas(obj.n_layers,2)-1) < 1E-6;
            end
            flag = islastbeta && islastdata;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkds, mlps] = eval_wo_lis(obj, func, x)
            [mllkds, mlps] = func(x, -1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [mllkds, mlps, mlogfx, gz_func, gz_ref, gz_irt] = ...
                eval_with_lis(obj, dirt, func, ratio_method, samples)
            mllkds = [];
            mlps = [];
            mlogfx = [];
            gz_func = [];
            gz_ref = [];
            gz_irt = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [V,s,g] = build_lis(obj, ratio_method, gz_func, gz_ref, gz_irt, log_weights)
            V = [];
            s = [];
            g = [];
        end
    end
end