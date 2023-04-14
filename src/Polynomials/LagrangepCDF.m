classdef LagrangepCDF < Lagrangep & PiecewiseCDF
    
    properties
        cheby Chebyshev2ndUnweighted
        cdf_basis2node(:,:)
    end
    
    methods (Static)
        function [cheby, cdf_basis2node] = lag2cheby(poly)
            %
            %Define data structure that maps lagrange polynomials to
            %Chebyshev polynomials with preserved boundary values. This
            %function is mainly useful for the squared mode--in which the
            %order of the polynomial is doubled.
            %
            %%%%%
            %Inputs:
            %
            %domain:
            %  domain of the input polynomial
            %
            %n_nodes:
            %  number of nodes used in the transformation
            %
            %
            %%%%%
            %Output:
            %
            %def:
            %  A data structure contains:
            %
            %  nodes:       nodes used for the transformation from Lagrange to Chebyshev
            %               we need to evaluate the pdf function on these nodes
            %  num_nodes:   number of nodes used in the transformation
            %  domain:      domain of the transformation
            %  basis2node:  inverse of the vandermonde matrix, transform nodal values
            %               to coefficient of
            %  node2basis:  Vandermonde matrix
            %
            %Tiangang Cui, August, 2019
            
            cheby = BoundedPolyCDF(poly);
            n_nodes = cardinal(cheby);
            if n_nodes < 3
                error('Must use more than three nodes')
            end
            
            tmp   = Chebyshev2ndUnweighted(n_nodes-3);
            %
            % need the jacobian
            % reference: [-1,1], target: [0,1], J = x2z = 0.5
            %
            ref_nodes = [cheby.domain(1); tmp.nodes(:); cheby.domain(2)];
            cheby.basis2node = eval_basis(cheby, ref_nodes);
            cheby.basis2node = cheby.basis2node;
            %
            [L,U] = lu(vpa(cheby.basis2node));
            cheby.node2basis  = double(U\(L\eye(n_nodes)));
            cheby.mass_R = []; % sqrt(cheby.weights./cheby.omegas).*cheby.basis2node;
            cheby.int_W  = []; % reshape(cheby.weights(:)./cheby.omegas(:), 1, [])*cheby.basis2node;
            %
            cheby.nodes = 0.5*(ref_nodes+1);
            %
            cdf_basis2node = eval_int_basis(cheby, ref_nodes(:));
            cdf_basis2node = cdf_basis2node*0.5;
        end
    end
    
    
    methods
        function obj = LagrangepCDF(poly, varargin)
            obj@Lagrangep(poly.order, poly.num_elems);
            obj@PiecewiseCDF(varargin{:});
            
            % local CDF poly
            [obj.cheby, obj.cdf_basis2node] = LagrangepCDF.lag2cheby(poly);
            %
            obj.mass    = [];
            obj.mass_R  = [];
            obj.int_W   = [];
            %
            % setup global nodes
            num_nodes   = obj.num_elems*(cardinal(obj.cheby)-1)+1;
            obj.nodes       = zeros(1, num_nodes);
            for i = 1:obj.num_elems
                ind = ( 1:cardinal(obj.cheby) ) + (cardinal(obj.cheby)-1)*(i-1);
                obj.nodes(ind) = obj.cheby.nodes*obj.elem_size + obj.grid(i);
            end
            
            % map the function value y to each local element
            if cardinal(obj.cheby) > 2
                j = cardinal(obj.cheby):(cardinal(obj.cheby)-1):num_nodes;
                obj.global2local = reshape([reshape(1:(num_nodes-1), cardinal(obj.cheby)-1, obj.num_elems); j], [], 1);
            else
                obj.global2local = reshape([1:(num_nodes-1); 2:num_nodes], [], 1);
            end
            % vectorise
            obj.global2local = reshape(obj.global2local, [], 1);
            obj.nodes = obj.nodes(:);
        end       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function data = pdf2cdf(obj, pdf)
            data.size = size(pdf,2);
            if data.size > 1
                % 1st coord: local, 2nd coord: elems, 3rd: pdfs
                local_pdf   = reshape(pdf(obj.global2local,:), cardinal(obj.cheby), obj.num_elems, data.size);
                % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
                local_pdf   = permute(local_pdf, [1,3,2]);
                data.poly_coef   = reshape(obj.cheby.node2basis*reshape(local_pdf, cardinal(obj.cheby), []), ...
                    cardinal(obj.cheby), data.size, obj.num_elems);
                data.cdf_poly_grid   = zeros(obj.num_elems+1, data.size);
                data.cdf_poly_nodes  = zeros(cardinal(obj), data.size);
                data.poly_base       = zeros(obj.num_elems, data.size);
                %
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*data.poly_coef(:,:,i))*obj.jac;
                    data.poly_base(i,:)  = tmp(1,:);
                    data.cdf_poly_nodes(ind(:,i),:)  = tmp - data.poly_base(i,:) + data.cdf_poly_grid(i,:);
                    data.cdf_poly_grid(i+1,:)        = data.cdf_poly_grid(i,:) + tmp(end,:) - data.poly_base(i,:);
                end
                %
                data.poly_norm = data.cdf_poly_grid(end,:);
            else
                % 1st coord: local, 2nd coord: elems
                local_pdf   = reshape(pdf(obj.global2local), cardinal(obj.cheby), obj.num_elems);
                data.poly_coef   = obj.cheby.node2basis*local_pdf;
                %
                data.cdf_poly_grid   = zeros(obj.num_elems+1, 1);
                data.cdf_poly_nodes  = zeros(cardinal(obj), 1);
                data.poly_base       = zeros(obj.num_elems, 1);
                %
                ind = reshape(obj.global2local, [], obj.num_elems);
                for i = 1:obj.num_elems
                    tmp             = (obj.cdf_basis2node*data.poly_coef(:,i))*obj.jac;
                    data.poly_base(i)    = tmp(1);
                    data.cdf_poly_nodes(ind(:,i)) = tmp - data.poly_base(i) + data.cdf_poly_grid(i);
                    data.cdf_poly_grid(i+1) = data.cdf_poly_grid(i) + tmp(end) - data.poly_base(i);
                end
                %
                data.poly_norm = data.cdf_poly_grid(end,:);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function F = eval_int_lag_local(obj, data, ei, mask, x)
            domains = [reshape(obj.grid(ei),[],1), reshape(obj.grid(ei+1),[],1)];
            x2z = 0.5*(domains(:,2)-domains(:,1));
            mid = 0.5*(domains(:,2)+domains(:,1));
            x   = (x(:)-mid)./x2z;
            if data.size == 1
                b   = eval_int_basis(obj.cheby, x); % rewrite
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.poly_coef(:,ei)',2), size(x));
                F   = tmp - reshape(data.poly_base(ei),size(tmp)) + reshape(data.cdf_poly_grid(ei),size(tmp));
            else
                pi  = find(mask);
                j1  = pi + (ei-1)*data.size;
                %
                b   = eval_int_basis(obj.cheby, x); % rewrite
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.poly_coef(:,j1)',2), size(x));
                %
                j2  = (pi-1)* obj.num_elems    + ei;
                j3  = (pi-1)*(obj.num_elems+1) + ei;
                F   = tmp - reshape(data.poly_base(j2),size(tmp)) + reshape(data.cdf_poly_grid(j3),size(tmp));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,df] = eval_int_lag_local_deri(obj, data, ei, mask, x)
            domains = [reshape(obj.grid(ei),[],1), reshape(obj.grid(ei+1),[],1)];
            x2z = 0.5*(domains(:,2)-domains(:,1));
            mid = 0.5*(domains(:,2)+domains(:,1));
            x   = (x(:)-mid)./x2z;
            if data.size == 1
                %b  = eval_ref_int_basis(obj.cheby, x); % rewrite
                [b,db] = eval_int_basis_newton(obj.cheby, x);
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.poly_coef(:,ei)',2), size(x));
                f   = tmp - reshape(data.poly_base(ei),size(tmp)) + reshape(data.cdf_poly_grid(ei),size(tmp));
                df  = reshape(sum(db.*data.poly_coef(:,ei)',2), size(x));
            else
                pi  = find(mask);
                j1  = pi + (ei-1)*data.size;
                %
                %b  = eval_ref_int_basis(obj.cheby, x); % rewrite
                [b,db] = eval_int_basis_newton(obj.cheby, x);
                b   = b.*x2z;
                tmp = reshape(sum(b.*data.poly_coef(:,j1)',2), size(x));
                %
                j2  = (pi-1)* obj.num_elems    + ei;
                j3  = (pi-1)*(obj.num_elems+1) + ei;
                f   = tmp - reshape(data.poly_base(j2),size(tmp)) + reshape(data.cdf_poly_grid(j3),size(tmp));
                df  = reshape(sum(db.*data.poly_coef(:,j1)',2), size(x));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        
        function r = invert_cdf_local(obj, data, ei, mask, rhs)
            if data.size == 1
                ind = sum(reshape(data.cdf_poly_nodes,1,[]) <= reshape(rhs,[],1), 2)';
            else
                ind = sum(data.cdf_poly_nodes(:,mask) <= reshape(rhs,1,[]), 1);
            end
            ind = max(ind, 1);
            ind = min(ind, numel(obj.nodes)-1);
            a = obj.nodes(ind);
            b = obj.nodes(ind+1);
            %
            %r = regula_falsi(obj, data, ei, mask, rhs(:), a(:), b(:));
            r = newton(obj, data, ei, mask, rhs(:), a(:), b(:));
            %
            %r = (r+1)/2;
        end
    end
    
end
