classdef SemiUnboundedDomain < LinearDomain
    
    methods
        function obj = SemiUnboundedDomain(varargin)
            defaultLeftBound = 0;
            defaultDxDz = 1;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %
            p = inputParser;
            addOptional(p,'left_bound',defaultLeftBound);
            addOptional(p,'dxdz',defaultDxDz,validScalarPosNum);
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            left = p.Results.left_bound;
            dxdz = p.Results.dxdz;
            %
            obj.bound = [left, inf];
            obj.mean = left;
            obj.dxdz = dxdz;
        end
    end
    
end