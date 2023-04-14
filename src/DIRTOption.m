classdef DIRTOption

    properties
        method
        max_layers
        %qmc_flag
        n_samples
        n_debugs
        defensive
        dhell_tol
    end

    properties (Access = private, Constant = true)
        defaultMethod   = 'Aratio';
        expectedMethod  = {'Eratio','Aratio'};
        defaultMLayers  = 50;
        %defaultQMCFlag  = false;
        defaultNSamples = 1E3;
        defaultNDebugs  = 1E3;
        defaultTau      = 1E-8;
    end

    methods
        function obj = DIRTOption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && all(x >= 0);
            %
            addParameter(p,'method',obj.defaultMethod, ...
                @(x) any(validatestring(x,obj.expectedMethod)));
            addParameter(p,'max_layers',obj.defaultMLayers, validScalarPosNum);
            %addParameter(p,'qmc_flag',  defaultQMCFlag, @(x) islogical(x) && isscalar(x));
            addParameter(p,'n_samples', obj.defaultNSamples,validScalarPosNum);
            addParameter(p,'n_debugs',  obj.defaultNDebugs, validScalarPosNum);
            %addParameter(p,'dhell_tol', obj.defaultDHellTol,validScalarPosNum);
            addParameter(p,'defensive', obj.defaultTau, validScalarPosNum);
            %
            p.KeepUnmatched = true;
            parse(p,varargin{:});
            %
            obj.method = p.Results.method;
            obj.max_layers = p.Results.max_layers;
            %obj.qmc_flag = p.Results.qmc_flag;
            obj.n_samples = p.Results.n_samples;
            obj.n_debugs = p.Results.n_debugs;
            %
            %obj.dhell_tol = p.Results.dhell_tol;
            obj.defensive = p.Results.defensive;
        end
    end
end