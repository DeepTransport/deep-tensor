classdef Lagrange1CDF < Lagrange1 & PiecewiseCDF
    
    properties
        T
        iV
    end
    
    methods
        function obj = Lagrange1CDF(poly, varargin)
            obj@Lagrange1(poly.num_elems);
            obj@PiecewiseCDF(varargin{:});
            %
            obj.nodes = linspace(obj.domain(1), obj.domain(2), obj.num_elems*2+1);
            num_nodes = length(obj.nodes);
            %
            dhh = obj.elem_size/2;
            %
            ii = zeros(3,obj.num_elems);
            jj = zeros(3,obj.num_elems);
            %
            % L = [1, 0, 0; 1, 1, 0; 1, 2, 1];
            % U = [1, 0, 0; 0, dhh, dhh^2; 0, 0, dhh^2*2];
            % LU = Vandermore
            % Vandermore = [1, 0, 0; 1, dhh, dhh^2; 1, dhh*2, 4*dhh^2];
            obj.iV = [1, 0, 0; -3/(2*dhh), 2/dhh, -1/(2*dhh); 1/(2*dhh^2), -1/dhh^2, 1/(2*dhh^2)];
            for i = 1:obj.num_elems
                ind = (1:3)+(i-1)*3;
                ii(:,i) = ind;
                jj(:,i) = (i-1)*2 + (1:3);
            end
            obj.T = sparse(ii(:), jj(:), ones(obj.num_elems*3,1), obj.num_elems*3, num_nodes);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function data = pdf2cdf(obj, pdf)
            data.size = size(pdf,2);
            if data.size > 1
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                data.poly_coef = obj.iV*reshape(obj.T*pdf, 3, []);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                % data.poly_coef = permute(reshape(obj.TT*pdf, 3, obj.num_elems, []), [1,3,2]);
                %
                cdf_elems = reshape([obj.elem_size, obj.elem_size^2/2, obj.elem_size^3/3]*reshape(data.poly_coef, 3, []), obj.num_elems, []);
                %
                data.cdf_poly_grid  = zeros(obj.num_elems+1, data.size);
                data.cdf_poly_grid(2:end,:) = cumsum(cdf_elems,1);
                %
                data.poly_norm = data.cdf_poly_grid(end,:);
                %
            else
                data.poly_coef = obj.iV*reshape(obj.T*pdf, 3, []);
                %
                cdf_elems = [obj.elem_size, obj.elem_size^2/2, obj.elem_size^3/3]*reshape(data.poly_coef, 3, []);
                %
                data.cdf_poly_grid  = zeros(obj.num_elems+1, 1);
                data.cdf_poly_grid(2:end) = cumsum(cdf_elems);
                %
                data.poly_norm = data.cdf_poly_grid(end);
                %
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_int_lag_local(obj, data, ei, mask, r)
            r = r - reshape(obj.grid(ei),size(r));
            if data.size > 1
                coi = find(mask);
                ind = ei + (coi-1)*obj.num_elems;
                jnd = ei + (coi-1)*(obj.num_elems+1);
                %
                % A: data.A(ind), B: data.B(ind), C: data.C(ind),
                % cdf: data.cdf_grid(jnd), base: data.base(ind)
                %
                % z = A(:).*r(:) + B(:).*r(:).^2/2 + C(:).*r(:).^3/3;
                % z = z + cdf(:);
                %
                f = sum([r(:), r(:).^2/2, r(:).^3/3].*data.poly_coef(:,ind)', 2) + reshape(data.cdf_poly_grid(jnd),[],1);
            else
                f = sum([r(:), r(:).^2/2, r(:).^3/3].*data.poly_coef(:,ei)', 2) + reshape(data.cdf_poly_grid(ei),[],1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,df] = eval_int_lag_local_deri(obj, data, ei, mask, r)
            r = r - reshape(obj.grid(ei),size(r));
            if data.size > 1
                coi = find(mask);
                ind = ei + (coi-1)*obj.num_elems;
                jnd = ei + (coi-1)*(obj.num_elems+1);
                %
                % A: data.A(ind), B: data.B(ind), C: data.C(ind),
                % cdf: data.cdf_grid(jnd), base: data.base(ind)
                %
                % z = A(:).*r(:) + B(:).*r(:).^2/2 + C(:).*r(:).^3/3;
                % z = z + cdf(:);
                %
                f  = sum([r(:), r(:).^2/2, r(:).^3/3].*data.poly_coef(:,ind)', 2) + reshape(data.cdf_poly_grid(jnd),[],1);
                df = sum([ones(length(r),1), r(:), r(:).^2].*data.poly_coef(:,ind)', 2);
            else
                f  = sum([r(:), r(:).^2/2, r(:).^3/3].*data.poly_coef(:,ei)', 2) + reshape(data.cdf_poly_grid(ei),[],1);
                df = sum([ones(length(r),1), r(:), r(:).^2].*data.poly_coef(:,ei)', 2);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf_local(obj, data, ei, mask, rhs)
            %
            a = obj.grid(ei);
            b = obj.grid(ei+1);
            %
            r = newton(obj, data, ei, mask, rhs(:), a(:), b(:));
            %if sum( r>b(:) | r<a(:) ) ~=0
            %    warning('newton failed')
            %    r = regula_falsi(obj, data, ei, mask, rhs(:), a(:), b(:));
            %end
        end

    end
    
end
