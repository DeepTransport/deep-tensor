classdef IRT < AbstractIRT
    % IRT class
    %
    % IRT Properties:
    %   int_dir - 
    %           The direction for marginalising the approximation
    %           >0: marginalise from x_d to x_1
    %           <0: marginalise from x_1 to x_d
    %   fun_z - Normalising constant of the function approximation part
    %   z     - Normalising constant
    %   tau   - denfensive term
    %   oned_cdfs - 
    %           One dimensional bases for building CDFs.
    %   approx  - 
    %           Approximation of the sqrt of the target density.
    %
    % IRT Methods:
    %   marginalise - 
    %           Marginalises the approximation.
    %   eval_pdf  - 
    %           Evaluates the normalised (marginal) pdf.
    %   eval_irt  - 
    %           X = R^{-1}(Z), where X is the target random variable, R is
    %           the Rosenblatt transport, and Z is the uniform random variable.
    %           * Can map marginal random variables for TT approximations.
    %           * Cannot map marginal random variables for sparse polynomial 
    %             approximations.
    %   eval_cirt - 
    %           Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly follow
    %           the target represented by SIRT, Z is uniform.
    %           * This function cannot handle marginal random variables.
    %   eval_rt   - 
    %           Z = R(X), where Z is uniform and X is target.
    %           * Can map marginal random variables for TT approximations.
    %           * Cannot map marginal random variables for sparse polynomial 
    %             approximations.
    %   eval_rt_jac - 
    %           Evaluates the Jacobian of Z = R(X).
    %           * This function cannot handle marginal random variables.
    %   random  - 
    %           Generates random variables
    %   sobol - Generates transformed sobol points
    %   set_defensive - 
    %           Resets the defensive term.
    %
    %%%%%%%%%%%%%%%%%
    %
    % see also TTIRT and PolyIRT
    
    methods (Static)
        function y = get_potential_to_density(base, y, z)
            [~,dxdz] = reference2domain(base, z);
            % log density of the reference measure
            mlogw = eval_measure_potential_reference(base, z);
            %
            % y is the potential function, use change of coordinates
            % y = exp( - 0.5*y + 0.5*sum(log(base.dxdz)) + 0.5*mlogw );
            y = exp( - (y-mlogw-sum(log(dxdz),1)) );
        end
        
        function y = potential_to_density(base, potential_fun, z)
            [x, dxdz] = reference2domain(base, z);
            y = feval(potential_fun, x);
            % log density of the reference measure
            mlogw = eval_measure_potential_reference(base, z);
            %
            % y is the potential function, use change of coordinates
            % y = exp( - 0.5*y + 0.5*sum(log(base.dxdz)) + 0.5*mlogw );
            y = exp( - (y-mlogw-sum(log(dxdz),1)) );
        end
    end
    
    methods
        function obj = IRT(potential, arg, varargin)
            [potential,arg,opt,var,tau,base] = parse_input(obj,potential,arg,varargin{:});
            %
            if all(arrayfun(@(x)isa(x{1},'Lagrange1'),base.oneds))
                func = @(z) IRT.potential_to_density(base, potential, z);
                obj.approx = approximate(obj, func, arg, opt, var);
                %
                obj.oned_cdfs = cell(size(obj.approx.base.oneds));
                for i = 1:length(obj.approx.base.oneds)
                    obj.oned_cdfs{i} = Lagrange1AbsCDF(obj.approx.base.oneds{i}, obj.approx.opt.cdf_tol);
                end
                obj.tau = tau;
                obj = marginalise(obj);
            else
                error('IRT only works for Lagrange1 basis, use SIRT')
            end
        end
    end
end