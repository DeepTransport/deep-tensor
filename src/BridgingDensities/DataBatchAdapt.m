classdef DataBatchAdapt < DataBatch
    
    properties
        init_mask
    end
    
    methods
        function obj = DataBatchAdapt(data_order,varargin)
            %
            p = inputParser;
            addRequired(p,'data_order');
            p.KeepUnmatched = true;
            parse(p,data_order,varargin{:});
            %
            obj@DataBatch(varargin{:});
            %
            obj.data_order = data_order;
            obj.init_mask = 1;
            %
            obj.masks = {};
            obj.masks_inc = {};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = set_init(obj, mllkds_vec, etol)
            if nargin < 3
                etol = 0.8;
            end
            etol = min(etol,obj.ess_tol_init);
            %
            end_i = 1;
            ind = obj.data_order(end_i);
            ess = ess_ratio(-sum(mllkds_vec(ind,:),1));
            while ess > etol && end_i < numel(obj.data_order)
                end_i = end_i+1;
                ind = obj.data_order(1:end_i);
                ess = ess_ratio(-sum(mllkds_vec(ind,:),1));
            end
            end_i = max(end_i-1, 1);
            obj.init_mask = 1:end_i; % update
            %
            ind = obj.data_order(obj.init_mask);
            mllkds = sum(mllkds_vec(ind,:),1);
            % define beta_p
            beta_p = 0;
            beta = max(beta_p,obj.min_beta);
            % compute ess over sample size
            ess = ess_ratio((beta_p-beta)*mllkds);
            while ess > etol
                beta = beta*obj.beta_factor;
                ess = ess_ratio((beta_p-beta)*mllkds);
            end
            %out.beta = min(1, beta/obj.beta_factor);
            beta = min(1, beta);
            %
            log_weights = -beta*mllkds;
            log_weights = log_weights - max(log_weights);
            ess = ess_ratio(log_weights);
            %
            fprintf('\nInitial: beta=%3.3e, ess=%3.3e\n', beta, ess);
            %
            obj.init_beta = beta;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = adapt_data(obj, method, mllkds_vec, mlps, mlogfx)
            % define previous mask
            if obj.data_layers > 0
                start_i = obj.masks{obj.data_layers}(end) + 1;
            else
                start_i = 1;
            end
            %
            if start_i > numel(obj.data_order)
                return;
            end
            end_i = start_i;
            switch method
                case {'Aratio'}  
                    ind = obj.data_order(end_i);
                    ess = ess_ratio(-sum(mllkds_vec(ind,:),1));
                    while ess > obj.ess_tol && end_i < numel(obj.data_order)
                        end_i = end_i+1;
                        ind = obj.data_order(start_i:end_i);
                        ess = ess_ratio(-sum(mllkds_vec(ind,:),1));
                    end
                case {'Eratio'}
                    ind = obj.data_order(1:end_i);
                    ess = ess_ratio(-sum(mllkds_vec(ind,:),1)-mlps+mlogfx);
                    while ess > obj.ess_tol && end_i < numel(obj.data_order)
                        end_i = end_i+1;
                        ind = obj.data_order(1:end_i);
                        ess = ess_ratio(-sum(mllkds_vec(ind,:),1)-mlps+mlogfx);
                    end
            end
            %
            if ess < obj.ess_tol
                end_i = max(end_i-1, start_i);
            end
            % add new data
            obj.data_layers = obj.data_layers +1;
            obj.masks{obj.data_layers} = 1:end_i;
            obj.masks_inc{obj.data_layers} = start_i:end_i;
            %
            obj.beta_p = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = adapt_density(obj, method, mllkds_vec, mlps, mlogfx)
            if obj.data_layers == 0 && obj.n_layers == 0
                obj.data_layers = 1;
                obj.masks{obj.data_layers} = obj.init_mask;
                obj.masks_inc{obj.data_layers} = obj.init_mask;
                %
                obj.betas = [0, obj.init_beta];
                obj.n_layers = obj.n_layers + 1;
                obj.beta_p = obj.init_beta;
            else
                if abs(obj.beta_p-1) < 1E-6
                    obj = adapt_data(obj, method, mllkds_vec, mlps, mlogfx);
                end
                obj = adapt_temperature(obj, method, mllkds_vec, mlps, mlogfx);
            end
        end
        
    end
    
end