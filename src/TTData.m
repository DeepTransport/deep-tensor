classdef TTData
    % TTData class: data of a functional TT approximation
    %
    % TTData Properties:
    %   direction -
    %           ALS direction, >0: built from left to right
    %                          <0: built from right to left
    %   cores - Nodal values or coefficent tensors of 1D functions, the
    %           dimension of the current core is organised as: previous
    %           rank x oned{k}.num_nodes x current rank.
    %   interp_x  -   
    %           Interpolation coordinates for functional TT.
    %   res_x - Interpolation coordinates for the residual.
    %   res_w - Interpolation weights for the residual.
    %   nout  - Dimension of the codomain.
    %
    % TTData Methods:
    %   clean - Cleans intermediate data used for builing TT and keeps 
    %           cores and direction for evaluations.
    %   rank  - Ranks of TT cores
    % 
    % See also TTFUN

    properties
        direction
        %
        cores
        interp_x
        %
        res_x
        res_w
        %
        nout
    end
    
    
    methods
        function obj = clean(obj)
            obj.interp_x = [];
            obj.res_x = [];
            obj.res_w = [];
        end
        
        function rs = rank(obj)
            % Find the ranks of the TT cores and the degrees of freedom of
            % the approximation basis for each coordinate.
            %   rs = SIZE(tt)
            %
            %   rs  - ranks, r(0) = 1 is not included
            %
            nd = length(obj.cores);
            rs  = ones(1,nd);
            for k = 1:nd-1
                rs(k) = size(obj.cores{k+1}, 1);
            end
        end
        
        function obj = TTData()
            obj.direction = 0;
            %
            obj.cores = [];
            obj.interp_x = [];
            %
            obj.res_x = [];
            obj.res_w = [];
            %
            obj.nout = 0;
        end

    end
end