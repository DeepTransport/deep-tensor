classdef BoundedDomain < LinearDomain
    
    methods
        function obj = BoundedDomain(varargin)
            defaultBound = [-1,1];
            p = inputParser;
            addOptional(p,'bound',defaultBound);
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            bound = p.Results.bound;
            if bound(1) < bound(2)
                obj.bound = bound;
                obj.mean = 0.5*(bound(:,2)+bound(:,1));
                obj.dxdz = 0.5*(bound(:,2)-bound(:,1));
            else
                error('left boundary should be less than the right boundary')
            end
        end
    end
    
end