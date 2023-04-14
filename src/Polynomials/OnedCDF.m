classdef OnedCDF
    % onedCDF class - Superclass for all one dimensional basis for
    %                 evaluating CDF and inverse of the CDF
    %
    % onedCDF Properties:
    %   tol             - Tolerance for inverting the CDF
    %
    % onedCDF Methods:
    %   invert_cdf      - Inverts CDF. Newton's method is used. If Newton's
    %                     method does not converge in 10 iterations, the
    %                     regula falsi method is then applied
    %   eval_cdf        - Evaluates CDF.
    %   eval_cdf_deri   - Evaluates the derivative of the conditional CDF.
    %                     This function is used for computing the Jacobian
    %                     of the Rosenblatt transport.
    %
    % See also PIECEWISECDF and SPECTRALCDF
    
    properties
        tol
        num_Newton
    end
    
    methods (Abstract)
        invert_cdf(obj)
        eval_cdf(obj)
        eval_int_deri(obj)
    end
    
    methods
        function obj = OnedCDF(varargin)
            defaultErrTol = 1E-8;
            defaultNumNewton = 10;
            p = inputParser;
            %
            addOptional(p,'err_tol',defaultErrTol, @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1));
            addOptional(p,'num_Newton',defaultNumNewton, @(x) isnumeric(x) && isscalar(x) && (x>0));
            parse(p,varargin{:});
            obj.tol = p.Results.err_tol;
            obj.num_Newton = p.Results.num_Newton;
        end
    end
end