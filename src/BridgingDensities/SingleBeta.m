classdef SingleBeta < Bridging
    
    properties
        min_beta
        ess_tol
        ess_tol_init
        beta_factor
        %
        betas
        init_beta
    end
    
    methods
        function obj = SingleBeta(varargin)
            %
            defaultMinBeta  = 1E-4;
            defaultESSTol   = 0.5;
            defaultBetaFactor = 1.05;
            %
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && all(x >= 0);
            %
            addOptional(p,'null',[]);
            addParameter(p,'min_beta',defaultMinBeta,validScalarPosNum);
            addParameter(p,'ess_tol',defaultESSTol,validScalarPosNum);
            addParameter(p,'ess_tol_init',defaultESSTol,validScalarPosNum);
            addParameter(p,'beta_factor',defaultBetaFactor,validScalarPosNum);
            %
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            %
            obj.min_beta = p.Results.min_beta;
            obj.ess_tol = p.Results.ess_tol;
            obj.ess_tol_init = p.Results.ess_tol_init;
            obj.beta_factor = p.Results.beta_factor;
            obj.init_beta = obj.min_beta;
        end
    end
    
end