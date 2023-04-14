classdef Bridging
    
    properties
        is_adapt
        n_layers
    end
    
    methods (Abstract)
        ratio_fun(obj, func, dirt, z, V)
        %
        get_ratio_fun(obj, dirt, z, mllkds, mlps, mlogfx)
        %
        eval_wo_lis(obj, func, x)
        %
        eval_with_lis(obj, dirt, func, ratio_method, samples)
        %
        set_init(obj, mllkds, etol)
        %
        adapt_density(obj, method, mllkds, mlps, mlogfx)
        %
        print_str(obj, mllkds_vec, mlps, mlogfx)
        %
        islast(obj)
    end
    
    methods
        function n = num_layers(obj)
            n = obj.n_layers;
        end
    end
end