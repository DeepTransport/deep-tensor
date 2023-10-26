classdef LISOption 
    % Set the options for constructing functional tensor train
    %
    % LISOption Properties:
    %   tol         - Truncation tolerance
    %   min_rank    - Min LIS rank 
    %   rank_frac   - Max LIS rank, as a fraction of the total dimension
    %   method      - One of 'reduction', 'unitary', and 'none'
    %
    
    properties
        tol
        min_rank
        max_rank
        rank_frac
        method
        ratio_method
    end
    
    properties (Access = private)
        defaultTol = 1E-2;
        defaultMinRank = 2;
        defaultRankFrac = 0.2;
        defaultMethod = 'none';
        expectedMethod  = {'none','reduction','unitary','permute'};
        defaultRaMethod   = 'Aratio';
        expectedRaMethod  = {'Eratio','Aratio'};
    end
    
    methods
        function obj = LISOption(varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'tol',  obj.defaultTol,  validErrTol);
            addParameter(p,'min_rank',obj.defaultMinRank,validScalarPosNum);
            addParameter(p,'rank_frac',obj.defaultRankFrac,validScalarPosNum);
            addParameter(p,'method',obj.defaultMethod, ...   
                @(x) any(validatestring(x,obj.expectedMethod)));
            addParameter(p,'ratio_method',obj.defaultRaMethod, ...
                @(x) any(validatestring(x,obj.expectedRaMethod)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            obj.tol   = tmp.tol;
            obj.min_rank = tmp.min_rank;
            obj.max_rank = tmp.min_rank; % unset initially 
            obj.rank_frac = tmp.rank_frac;
            obj.method = tmp.method;
            obj.ratio_method = tmp.ratio_method;
        end
        
        function obj = update(obj, varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validErrTol = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<1);
            %
            addParameter(p,'tol',  obj.defaultALSTol,  validErrTol);
            addParameter(p,'min_rank', obj.defaultMinRank,validScalarPosNum);
            addParameter(p,'max_rank');
            addParameter(p,'rank_frac',obj.defaultRankFrac,validScalarPosNum);
            addParameter(p,'method',obj.defaultMethod, ...   
                @(x) any(validatestring(x,obj.expectedMethod)));
            addParameter(p,'ratio_method',obj.defaultRaMethod, ...
                @(x) any(validatestring(x,obj.expectedRaMethod)));
            %
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            tmp = p.Results;
            %
            if ~ismember('tol',cellstr(p.UsingDefaults))
                obj.tol = tmp.tol;
            end
            if ~ismember('min_rank',cellstr(p.UsingDefaults))
                obj.min_rank = tmp.min_rank;
            end
            if ~ismember('rank_frac',cellstr(p.UsingDefaults))
                obj.rank_frac = tmp.rank_frac;
            end
            if ~ismember('max_rank',cellstr(p.UsingDefaults))
                obj.max_rank = tmp.max_rank;
            end
            if ~ismember('method',cellstr(p.UsingDefaults))
                obj.method = tmp.method;
            end
            if ~ismember('ratio_method',cellstr(p.UsingDefaults))
                obj.ratio_method = tmp.ratio_method;
            end
        end

        function obj = set_max_rank(obj, r)
            if r < obj.min_rank
                disp('the max LIS rank is less than the min LIS rank');
            end
            obj.max_rank = max(r, obj.min_rank);
        end
    end
end