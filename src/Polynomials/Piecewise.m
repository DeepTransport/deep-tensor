classdef Piecewise < Oned
% Piecewise class   - Superclass for Lagrange1 (piecewise linear basis) and 
%                     Larangep (piecewise high order Lagrange)
%
% Constructors:
%   * Lagrange1(n, domain)
%   * Lagrangep(order, n, domain)
%
%   order       - (Required) order of the polynomial 
%   n           - (Required) number of elements
%   domain      - (Optional) default is [0,1]
%
% See also PIECEWISECDF
    
    properties
        %{
        gs % ghostsize
        bc % boundary condition
        %}        
        grid(:,1) 
        num_elems 
        elem_size 
        jac 
        mass(:,:) 
        unweighed_mass(:,:) 
        %
        unweighed_mass_R(:,:) 
        unweighed_int_W(1,:) 
    end
    
    methods
        function obj = Piecewise(order, num_elems)
            %[obj.order,obj.num_elems,obj.domain,obj.bc,obj.gs] = ...
            %    Piecewise.process_input(order, num_elems, varargin{:});
            %
            obj.domain = [-1,1];
            obj.order = order;
            obj.num_elems = num_elems;
            obj.grid = linspace(obj.domain(1), obj.domain(2), obj.num_elems+1);
            %obj.grid = linspace(obj.domain(1)+obj.gs, obj.domain(2)-obj.gs, obj.num_elems+1);
            obj.elem_size = (obj.grid(end)-obj.grid(1))/double(obj.num_elems);
            obj.constant_weight = true;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function coeff = approximate(obj, func, nfout)
            coeff = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function x = sample_measure(obj, n)
            x = rand(1,n)*(obj.grid(end)-obj.grid(1))+obj.grid(1);
        end

        function x = sample_measure_skip(obj, n)
            x = sample_measure(obj, n);
        end

        function w = eval_measure(obj, x)
            w = ones(size(x)) * ( 1/(obj.grid(end)-obj.grid(1)) );
        end        

        function w = eval_log_measure(obj, x)
            w = ones(size(x)) *(-log(obj.grid(end)-obj.grid(1)));
        end   

        function w = eval_measure_deri(obj, x)
            w = zeros(size(x));
        end       

        function w = eval_log_measure_deri(obj, x)
            w = zeros(size(x));
        end     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ind,B,interp_atx] = point_selection(obj, int_method, H)
            % Build the cross indices
            rold = size(H,1)/cardinal(obj);
            %
            switch int_method
                case {'QDEIM'}
                    [~,~,e] = qr(H', 'vector');
                    ind = e(1:size(H,2));
                    interp_atx = H(ind,:);
                    B =  H/interp_atx;
                case {'DEIM'}
                    disp('Not implemented')
                case {'MaxVol'}
                    [ind,B] = maxvol(H);
                    interp_atx = H(ind,:);
            end
            if cond(interp_atx) > 1E5
                disp('warning: poor condition number in the interpolation')
            end
        end
        
    end
    
end
