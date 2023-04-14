classdef SparseOption
    % Set the options for constructing sparse polynomial
    %

    properties
        init_sample_size
        max_sample_size
        enrich_sample_size
        %
        bulk_parameter
        max_dim_basis
        init_total_degree
        enrich_degree
        %
        adaptation_rule
        weight_rule
        indexset
        fast
        %
        tol
        opt_tol
        stagnation_tol
        overfit_tol
        cdf_tol
        %
        display_iterations
    end

    properties (Access = private, Constant = true)
        defaultInitSampleSize = 3;
        defaultMaxSampleSize = 1E5;
        defaultEnrichSampleSize = 1;
        %
        defaultBulkParameter = 0.5;
        defaultMaxDimBasis = 1E3;
        defaultInitTotalDegree = 10;
        defaultEnrichDegree = 2;
        %
        defaultTol = 1E-2;
        defaultOptTol = 2;
        defaultStagnationTol = 5E-2;
        defaultOverfitTol = 1.1;
        defaultCDFTol   = 1E-8;
        %
        defaultAdaptationRule  = 'reducedmargin';
        expectedAdaptationRule = {'margin', 'reducedmargin'};
        %
        defaultWeightRule  = 'Christoffel';
        expectedWeightRule = {'none', 'Christoffel'};
        %
        defaultIndexSet  = 'total';
        expectedIndexSet = {'total', 'hyperbolic'};
        %
        defaultDisplayIterations = true;
        %
        defaultFastUpdate = false;
    end

    methods
        function obj = SparseOption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'init_sample_size',  obj.defaultInitSampleSize, ...
                validScalarPosNum);
            addParameter(p,'max_sample_size',  obj.defaultMaxSampleSize, ...
                validScalarPosNum);
            addParameter(p,'enrich_sample_size', obj.defaultEnrichSampleSize, ...
                validScalarPosNum);
            %
            addParameter(p,'bulk_parameter', obj.defaultBulkParameter);
            addParameter(p,'max_dim_basis', obj.defaultMaxDimBasis);
            addParameter(p,'init_total_degree', obj.defaultInitTotalDegree);
            addParameter(p,'enrich_degree', obj.defaultEnrichDegree);
            %
            addParameter(p,'tol', obj.defaultTol, validErrTol);
            addParameter(p,'opt_tol', obj.defaultOptTol, ...
                @(x) isnumeric(x) && isscalar(x) && (x>=0));
            addParameter(p,'stagnation_tol', obj.defaultStagnationTol, validErrTol);
            addParameter(p,'overfit_tol', obj.defaultOverfitTol, ...
                @(x) isnumeric(x) && isscalar(x) && (x>=0));
            addParameter(p,'cdf_tol', obj.defaultCDFTol, validErrTol);
            %
            addParameter(p,'display_iterations', obj.defaultDisplayIterations);
            %
            addParameter(p,'adaptation_rule', obj.defaultAdaptationRule, ...
                @(x) any(validatestring(x,obj.expectedAdaptationRule)));
            %
            addParameter(p,'fast_update', obj.defaultFastUpdate);
            %
            addParameter(p,'weight_rule', obj.defaultWeightRule, ...
                @(x) any(validatestring(x,obj.expectedWeightRule)));
            %
            addParameter(p,'indexset', obj.defaultIndexSet, ...
                @(x) any(validatestring(x,obj.expectedIndexSet)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            %
            obj.init_sample_size = p.Results.init_sample_size;
            obj.max_sample_size = p.Results.max_sample_size;
            obj.enrich_sample_size = p.Results.enrich_sample_size;
            %
            obj.bulk_parameter = p.Results.bulk_parameter;
            obj.max_dim_basis = p.Results.max_dim_basis;
            obj.init_total_degree = p.Results.init_total_degree;
            obj.enrich_degree = p.Results.enrich_degree;
            obj.adaptation_rule = p.Results.adaptation_rule;
            obj.weight_rule = p.Results.weight_rule;
            obj.indexset = p.Results.indexset;
            obj.fast = p.Results.fast_update;
            %
            obj.tol = p.Results.tol;
            obj.opt_tol = p.Results.opt_tol;
            obj.stagnation_tol = p.Results.stagnation_tol;
            obj.overfit_tol = p.Results.overfit_tol;
            obj.cdf_tol = p.Results.cdf_tol;
            obj.display_iterations = p.Results.display_iterations;
        end
    end
end