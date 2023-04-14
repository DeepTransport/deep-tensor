classdef Spectral < Oned
% Spectral class    - Superclass for Fourier, Chebyshev1st, Chebyshev2nd,
%                     and recurr subclasses. The recurr class implements
%                     Legendre, Jabobi11, Hermite, and Laguerre subclasses
%
% Constructors:
%   * poly(order, domain)
%
%   poly        - Choose one from Fourier, Chebyshev1st, Chebyshev2nd, 
%                 Legendre, Jabobi11, Hermite, and Laguerre. 
%   order       - (Required) order of the polynomial 
%   domain      - (Optional) default is [0,1]
%
% See also SPECTRALCDF
    
    properties
        weights(:,1) 
        omegas(:,1) 
        normalising(1,:) 
        basis2node(:,:) 
        node2basis(:,:) 
    end
    
    methods        
        function obj = post_construction(obj)
            obj.basis2node = eval_basis(obj, obj.nodes);
            obj.omegas = eval_measure(obj, obj.nodes);
            obj.omegas = reshape(obj.omegas, size(obj.nodes));
            obj.weights = reshape(obj.weights, size(obj.nodes));
            obj.node2basis  = obj.basis2node'*diag(obj.weights);
            %
            obj.mass_R  = eye(cardinal(obj));
            obj.int_W   = reshape(obj.weights(:), 1, [])*obj.basis2node;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function coeff = approximate(obj, func, nfout)
            coeff = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        function coeff = collocation(obj, func)
            f = func(obj.nodes);
            coeff = obj.node2basis*(f./obj.omegas);
        end
        %}
        

        function debug(obj, xs)
            df = eval_log_measure_deri(obj, xs);
            f = eval_log_measure(obj, xs);
            fd = (f(3:end)-f(1:end-2))./(xs(3:end)-xs(1:end-2));
            figure
            subplot(1,2,1)
            plot(xs(2:end-1), df(2:end-1))
            hold on
            plot(xs(2:end-1), fd)    
            title('log density')
            norm(df(2:end-1) - fd)
            
            df = eval_measure_deri(obj, xs);
            f = eval_measure(obj, xs);
            fd = (f(3:end)-f(1:end-2))./(xs(3:end)-xs(1:end-2));
            subplot(1,2,2)
            plot(xs(2:end-1), df(2:end-1))
            hold on
            plot(xs(2:end-1), fd)    
            title('density')
            norm(df(2:end-1) - fd)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ind,B,interp_atx] = point_selection(obj, int_method, H)
            % Build the cross indices
            %
            rold  = size(H,1)/cardinal(obj);
            nodes = reshape(obj.basis2node*reshape(H, cardinal(obj), []), length(obj.nodes), rold, []);
            nodes = reshape(nodes, rold*length(obj.nodes), []);
            switch int_method
                case {'QDEIM'}
                    [~,~,e] = qr(nodes', 'vector');
                    ind = e(1:size(nodes,2));
                    interp_atx = nodes(ind,:);
                    B =  H/interp_atx;
                case {'MaxVol'}
                    [ind,~] = maxvol(nodes);
                    interp_atx = nodes(ind,:);
                    B =  H/interp_atx;
            end
            if cond(interp_atx) > 1E5
                disp('warning: poor condition number in the interpolation')
            end
        end
    end
    
end
