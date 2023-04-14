classdef SparseSIRT < SIRT 
    % OrthogonalSIRT class
    %
    % OrthogonalSIRT Properties:
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also ONED, ONEDCDF, SPARSEOPTION, and SPARSEFUN
    
    % the current implementation only uses polynomial basis with bounded 
    % domain, all operations are mapped to the reference domain [-1,1]
    % z = (x - shifts)./scales; for Legendre, z in [-1,1]
    % x in domain
    
    properties
        bases_cdf_nodes
        cum_P
        batch_size
    end

    methods
        J = eval_rt_jac_reference_inner(obj, r, z)
    end

    methods        
        
        function obj = SparseSIRT(potential_func, base, varargin)
            defaultOption   = SparseOption();
            obj@SIRT(defaultOption, potential_func, base, varargin{:});
            obj.batch_size = 1E2;
            obj.int_dir = 0;
        end

        function approx = approximate(obj, func, base, opt, samples, deb, previous)
            if isempty(previous)
                approx = SparseFun(func, base, opt, 'debug', deb);
            elseif isa(previous, 'SparseFun')
                approx = SparseFun(func, base, opt, 'debug', deb, 'init_indices', previous.indices);
            else
                warning('previous data should be a SparseFun object')
                approx = SparseFun(func, base, opt, 'debug', deb);
            end
        end
    end
    
end