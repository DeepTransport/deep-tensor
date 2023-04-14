classdef UnboundedDomain < LinearDomain
    
    methods
        function obj = UnboundedDomain(varargin)
            defaultMean = 0;
            defaultDxDz = 1;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %
            p = inputParser;
            addOptional(p,'mean',defaultMean,@(x) isnumeric(x)&&isscalar(x));
            addOptional(p,'dxdz',defaultDxDz,validScalarPosNum);
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            mean = p.Results.mean;
            dxdz = p.Results.dxdz;
            %
            obj.bound = [-inf, inf];
            obj.mean = mean;
            obj.dxdz = dxdz;
        end
    end
    
end