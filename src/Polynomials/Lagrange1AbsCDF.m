classdef Lagrange1AbsCDF < Lagrange1 & OnedCDF
    
    properties
        T
        iV
    end
    
    methods
        function obj = Lagrange1AbsCDF(poly, varargin)
            obj@Lagrange1(poly.num_elems);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function p = eval_pdf_radon(obj, pk, r)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            r = reshape(r,[],1);
            data = pdf2cdf(obj, pk);
            p = zeros(size(r));
            %
            ei = sum(reshape(obj.nodes,1,[]) < r,2)';
            mask1 = ei==0;
            mask2 = ei==(obj.num_elems+1);
            mask3 = ~mask1 & ~mask2;
            %
            if sum(mask3) > 0
                ei = ei(mask3);
                %
                if data.size == 1
                    p(mask3) = r(mask3).*(reshape(data.a(ei),[],1)*2) + reshape(data.b(ei),[],1);
                    p(mask3) = p(mask3)./data.norm;
                else
                    coi = find(mask3);
                    ind = ei + (coi-1)*obj.num_elems; 
                    p(mask3) = r(mask3).*(reshape(data.a(ind),[],1)*2) + reshape(data.b(ind),[],1);
                    p(mask3) = p(mask3)./reshape(data.norm(mask3),[],1);
                end
            end
            w = eval_measure(obj, r(:));
            p = p./w;
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function data = pdf2cdf(obj, pk)
            % p(x) = 2ax + b
            % F(x) = ax^2 + bx + c
            %
            h = obj.nodes(2:end) - obj.nodes(1:end-1);
            a = (pk(2:end,:) - pk(1:end-1,:))./h;
            data.a = a/2;
            data.b = (pk(1:end-1,:).*obj.nodes(2:end)-pk(2:end,:).*obj.nodes(1:end-1))./h;
            %
            left  = obj.nodes(1:end-1).^2.*data.a + obj.nodes(1:end-1).*data.b;
            right = obj.nodes(2:end  ).^2.*data.a + obj.nodes(2:end  ).*data.b;
            tmp   = right-left;
            %
            data.cdf_grid_unnormalised_right = cumsum(tmp, 1);
            data.cdf_grid_unnormalised_left  = data.cdf_grid_unnormalised_right - tmp;
            data.norm = data.cdf_grid_unnormalised_right(end,:);
            %
            data.c = data.cdf_grid_unnormalised_left - left;
            data.size = size(pk,2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_int_lag1(obj, data, r)
            if data.size > 1 && data.size ~= length(r)
                error('Error: dimenion mismatch')
            end
            [mxi, nxi] = size(r);
            r = reshape(r,[],1);
            z = zeros(size(r));
            % index to locate the element
            ei = sum(reshape(obj.nodes,1,[]) < r,2)';
            mask1 = ei==0;
            %
            mask2 = ei==(obj.num_elems+1);
            if sum(mask2) > 0
                if data.size ==1
                    z(mask2) = data.norm;
                else
                    z(mask2) = data.norm(mask2);
                end
            end
            mask3 = ~mask1 & ~mask2;
            if sum(mask3) > 0
                ei = ei(mask3);
                if data.size > 1
                    coi = find(mask3);
                    ind = ei + (coi-1)*obj.num_elems; 
                    z(mask3) = r(mask3).^2.*reshape(data.a(ind),[],1) + r(mask3).*reshape(data.b(ind),[],1) + reshape(data.c(ind),[],1);
                else
                    z(mask3) = r(mask3).^2.*reshape(data.a(ei),[],1) + r(mask3).*reshape(data.b(ei),[],1) + reshape(data.c(ei),[],1);
                end
            end
            z = reshape(z, mxi, nxi);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            r = reshape(r,[],1);
            data = pdf2cdf(obj, pk);
            z = eval_int_lag1(obj, data, r);
            %
            %
            if numel(data.norm) > 1
                z = z./reshape(data.norm,size(z));
            else
                z = z./data.norm;
            end
            z = reshape(z, size(r));
            %
            z(isnan(z)) = eps;
            z(isinf(z)) = 1-eps;
            z(z>(1-eps)) = 1-eps;
            z(z<eps) = eps;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_int_deri(obj, pk, r)
            data = pdf2cdf(obj, pk);
            z = eval_int_lag1(obj, data, r);
            z = reshape(z, size(r));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [r,p] = invert_cdf(obj, pk, xi)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
                %pk = abs(pk);
            end
            data = pdf2cdf(obj, pk);
            %
            if data.size > 1 && data.size ~= length(xi)
                error('Error: dimenion mismatch')
            end
            %
            [mxi, nxi] = size(xi);
            xi = reshape(xi,[],1); % vertical
            %
            r = zeros(size(xi)); % vertical
            p = zeros(size(r));
            % index to locate the element
            if data.size == 1
                rhs = xi.*data.norm;
                ei  = sum(reshape(data.cdf_grid_unnormalised_right,1,[]) < rhs,2)';
            else
                rhs = xi(:).*data.norm(:);
                ei  = sum(data.cdf_grid_unnormalised_right <= reshape(rhs,1,[]), 1);
            end
            ei = ei + 1; % count the left cdf point (=0)
            %
            mask1 = ei==0; % left cell
            if sum(mask1) > 0
                r(mask1) = obj.domain(1);
            end
            mask2 = ei==(obj.num_elems+1); % right cell, need to implement
            if sum(mask2) > 0
                r(mask2) = obj.domain(end);
            end
            %
            mask3 = ~mask1 & ~mask2; % middle cells
            if sum(mask3) > 0
                ei = ei(mask3);
                if data.size == 1
                    a = reshape(data.a(ei),[],1);
                    b = reshape(data.b(ei),[],1);
                    c = reshape(data.c(ei),[],1) - reshape(rhs(mask3),[],1);
                else
                    coi = find(mask3);
                    ind = ei + (coi-1)*obj.num_elems; 
                    a = reshape(data.a(ind),[],1);
                    b = reshape(data.b(ind),[],1);
                    c = reshape(data.c(ind),[],1) - reshape(rhs(mask3),[],1);
                end
                r(mask3) = quadratic_root(a,b,c,obj.nodes(ei),obj.nodes(ei+1));
                p(mask3) = r(mask3).*(a*2) + b;
                %
                if data.size == 1
                    p(mask3) = p(mask3)./data.norm;
                else
                    p(mask3) = p(mask3)./reshape(data.norm(mask3),[],1);
                end
            end
            %
            r = reshape(r, mxi, nxi);
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            %
            w = eval_measure(obj, r(:));
            p = p./w;
        end
    end
    
    
end
