classdef TTOption 
    % Set the options for constructing functional tensor train
    %
    % TToption Properties:
    %   max_rank  - Max rank of each core, default is 20.
    %   init_rank - Rank of the initial tensor train, default is 10.
    %   kick_rank - Rank of the enrichment sample size, default is 2.
    %   max_als   - Max number of ALS iterations. Default is 4.
    %   als_tol   - Tolerance for terminating ALS. Default is 1E-2.
    %   local_tol - Truncation tolerance of local SVD, default is 1E-10.
    %               The SVD is truncated at singular values that is about
    %               1E-10 relative to the largest singular value.
    %   cdf_tol   - Tolerance for evaluating the inverse CDF function.
    %               Used by SIRT. Default is 1E-10.
    %   tt_method - Construction method. Default is option is 'amen'. 
    %               Options are 'amen' and 'random'.
    %   int_method  
    %             - Interpolation method for choosing cross indices. 
    %               Default is 'MaxVol'.
    %
    % TToption Methods:
    %   TToption  - Constructor. If no parameter is passed in, it returns 
    %               the default obj.
    %
    
    properties
        cdf_tol
        local_tol
        max_rank
        init_rank
        kick_rank
        max_als
        als_tol
        tt_method
        int_method
        %
        display_iterations
    end
    
    properties (Access = private, Constant = true)
        defaultMaxALS   = 4;
        defaultALSTol   = 1E-4;
        defaultInitRank = 20;
        defaultKickRank = 2;
        defaultMaxRank  = 30;
        defaultLocTol   = 1E-10;
        defaultCDFTol   = 1E-10;
        defaultTTMethod = 'amen';
        expectedTTMethod  = {'random','amen','fix_rank'};
        defaultIntM     = 'MaxVol';
        expectedIntM    = {'QDEIM','WDEIM','MaxVol'};
        defaultDisplayIterations = true;
    end
    
    methods
        function obj = TTOption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'max_als',  obj.defaultMaxALS,  validScalarPosNum);
            addParameter(p,'als_tol',  obj.defaultALSTol,  validErrTol);
            addParameter(p,'init_rank',obj.defaultInitRank,validScalarPosNum);
            addParameter(p,'kick_rank',obj.defaultKickRank,@(x) isnumeric(x) && isscalar(x) && (x >= 0));
            addParameter(p,'max_rank', obj.defaultMaxRank, validScalarPosNum);
            addParameter(p,'local_tol',obj.defaultLocTol,  validErrTol);
            addParameter(p,'cdf_tol',  obj.defaultCDFTol,  validErrTol);
            addParameter(p,'tt_method',obj.defaultTTMethod, ...   
                @(x) any(validatestring(x,obj.expectedTTMethod)));
            addParameter(p,'int_method',  obj.defaultIntM, ...   
                @(x) any(validatestring(x,obj.expectedIntM)));
            addParameter(p,'display_iterations', obj.defaultDisplayIterations);
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            obj.max_als   = tmp.max_als;
            obj.als_tol   = tmp.als_tol;
            obj.init_rank = tmp.init_rank;
            obj.kick_rank = tmp.kick_rank;
            obj.max_rank  = tmp.max_rank;
            obj.local_tol = tmp.local_tol;
            obj.cdf_tol   = tmp.cdf_tol;
            obj.tt_method = tmp.tt_method;
            obj.int_method  = tmp.int_method;
            obj.display_iterations = p.Results.display_iterations;
            
            if obj.kick_rank == 0
                obj.tt_method = 'fix_rank';
            end
        end
    end
end