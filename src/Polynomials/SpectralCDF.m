classdef SpectralCDF < OnedCDF
    % SpectralCDF class
    %
    % For Fourier basis, FourierCDF is used. For other spectral polynomials
    % in bounded domains, we first transform the polynomial to the 2nd
    % Chebyshev basis, and then apply the inversion. See Chebyshev2ndCDF.
    %
    % Before applying root findings, a grid search based on sampling_nodes
    % is applied to locate the left and right boundary of root finding.
    %
    % See also ChebyshevCDF and FourierCDF.
    
    properties
        sampling_nodes(1,:)
        cdf_basis2node(:,:)
    end
    

    methods (Abstract)
        grid_measure(obj)
        eval_int_basis(obj)
        eval_int_basis_newton(obj)
    end
    
    methods
        function obj = SpectralCDF(varargin)
            obj@OnedCDF(varargin{:});
            %
            obj.sampling_nodes = grid_measure(obj, max(cardinal(obj)*2, 200));
            %[x,J] = domain2reference(obj,obj.sampling_nodes(:));
            %b = eval_int_basis(obj, x);
            %obj.cdf_basis2node      = b.*J;
            obj.cdf_basis2node = eval_int_basis(obj, obj.sampling_nodes(:));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = update_sampling_nodes(obj, sampling_nodes)
            obj.sampling_nodes = sampling_nodes;
            obj.cdf_basis2node = eval_int_basis(obj, obj.sampling_nodes(:));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_int(obj, coef, x)
            %[x,J] = domain2reference(obj,x(:));
            %b = eval_int_basis(obj, x);
            %b = b.*J;
            b = eval_int_basis(obj, x);
            %
            if size(coef,2) > 1
                if size(coef,2) == length(x)
                    f = reshape(sum(b.*coef',2), [], 1);
                else
                    error('Error: dimenion mismatch')
                end
            else
                f = reshape((b*coef), [], 1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_int_search(obj, coef, cdf_poly_base, rhs, x)
            f = eval_int(obj, coef, x);
            f = f - cdf_poly_base - rhs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [f,df] = eval_int_newton(obj, coef, cdf_poly_base, rhs, x)
            %
            [b,db] = eval_int_basis_newton(obj, x(:));
            %
            if size(coef,2) > 1
                if size(coef,2) == length(x)
                    f  = reshape(sum(b .*coef',2), [], 1);
                    df = reshape(sum(db.*coef',2), [], 1);
                else
                    error('Error: dimenion mismatch')
                end
            else
                f  = reshape((b *coef), [], 1);
                df = reshape((db*coef), [], 1);
            end
            f  = f - cdf_poly_base - rhs;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = eval_cdf(obj, pk, r)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            if size(pk,2) > 1 && size(pk,2) ~= length(r)
                error('Error: dimenion mismatch')
            end
            r = reshape(r,[],1);
            %
            coef = obj.node2basis*pk;
            poly_base = obj.cdf_basis2node(1,:)*coef;
            poly_norm = (obj.cdf_basis2node(end,:)-obj.cdf_basis2node(1,:))*coef;
            %
            mask1 = r<=obj.sampling_nodes(1);
            mask3 = r>=obj.sampling_nodes(end);
            mask2 = ~(mask1|mask3);
            z = zeros(size(r));
            if sum(mask3) > 0
                if size(pk,2) == 1
                    z(mask3) = poly_norm;
                else
                    z(mask3) = poly_norm(mask3);
                end
            end
            if sum(mask2) > 0
                if size(pk,2) == 1
                    z(mask2) = eval_int(obj, coef, r(mask2)) - poly_base;
                else
                    tmp = eval_int(obj, coef(:,mask2), r(mask2));
                    z(mask2) = reshape(tmp,[],1) - reshape(poly_base(mask2),[],1);
                end
            end
            %
            if numel(poly_norm) > 1
                z = z./reshape(poly_norm, size(z));
            else
                z = z/poly_norm;
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
            r = reshape(r,[],1);
            coef = obj.node2basis*pk;
            base = obj.cdf_basis2node(1,:)*coef;
            if size(pk,2) == 1
                z = eval_int(obj, coef, r) - base;
                z = reshape(z, size(r));
            else
                tmp = eval_int(obj, coef, r);
                z = (reshape(tmp, size(r)) - reshape(base, size(r)));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function r = invert_cdf(obj, pk, xi)
            if (sum(pk(:)<0)>0)
                disp(['negative pdf ' num2str(sum(pk(:)<0))])
            end
            data_size = size(pk,2);
            coef = obj.node2basis*pk;
            cdf_poly_nodes = obj.cdf_basis2node*coef;
            %
            cdf_poly_base  = cdf_poly_nodes(1,:);
            cdf_poly_nodes = cdf_poly_nodes - cdf_poly_base;
            cdf_poly_norm  = cdf_poly_nodes(end,:);
            if data_size > 1 && data_size ~= length(xi)
                error('Error: dimenion mismatch')
            end
            %
            xi = reshape(xi,[],1);
            r = zeros(size(xi));
            
            if data_size == 1
                rhs = xi.*cdf_poly_norm; % vertical
                ind = sum(reshape(cdf_poly_nodes,1,[]) < rhs(:),2)';
            else
                rhs = xi(:).*cdf_poly_norm(:); % vertical
                ind = sum(cdf_poly_nodes < reshape(rhs,1,[]), 1);
            end
            mask1 = ind==0 | reshape(xi,1,[])<=eps;
            mask3 = ind==length(obj.sampling_nodes) | reshape(xi,1,[])>=1-eps;
            mask2 = ~(mask1|mask3);
            %
            % left and right tails
            if sum(mask1) > 0
                r(mask1) = obj.sampling_nodes(1);
            end
            if sum(mask3) > 0
                r(mask3) = obj.sampling_nodes(end);
            end
            %
            if sum(mask2) > 0
                a = obj.sampling_nodes(ind(mask2));
                b = obj.sampling_nodes(ind(mask2)+1);
                %
                if data_size == 1
                    r(mask2) = newton(obj, coef, cdf_poly_base, rhs(mask2), a(:), b(:));
                else
                    r(mask2) = newton(obj, coef(:,mask2), reshape(cdf_poly_base(mask2),[],1), rhs(mask2), a(:), b(:));
                end
            end
            %
            r = reshape(r, size(xi));
            r(isnan(r)) = 0.5*(obj.domain(1) + obj.domain(2));
            r(isinf(r)) = 0.5*(obj.domain(1) + obj.domain(2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c = regula_falsi(obj, coef, cdf_poly_base, rhs, a, b)
            fa = eval_int_search(obj,coef,cdf_poly_base,rhs,a);
            fb = eval_int_search(obj,coef,cdf_poly_base,rhs,b);
            if any((fb.*fa) > 0)
                disp(['Root finding: initial guesses on one side, # violation: ' num2str(sum((fb.*fa) > 0))])
                %disp([a(fb.*fa > 0), b(fb.*fa > 0)])
            end
            c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
            cold = inf;
            %i = 2;
            while ( norm(c-cold, Inf) > obj.tol )
                cold = c;
                fc  = eval_int_search(obj,coef,cdf_poly_base,rhs,c);
                if norm(fc, Inf) < obj.tol
                    break;
                end
                I1  = (fc < 0);
                I2  = (fc > 0);
                I3  = ~I1 & ~I2;
                a   = I1.*c + I2.*a + I3.*c;
                b   = I1.*b + I2.*c + I3.*c;
                fa  = I1.*fc + I2.*fa + I3.*fc;
                fb  = I1.*fb + I2.*fc + I3.*fc;
                step    = -fb.*(b - a)./(fb - fa);
                step(isnan(step)) = 0;
                c = b + step;
                %norm(fc, inf)
                %i = i+1;
                %disp(i)
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c = newton(obj, coef, cdf_poly_base, rhs, a, b)
            %i = 0;
            fa = eval_int_search(obj,coef,cdf_poly_base,rhs,a);
            fb = eval_int_search(obj,coef,cdf_poly_base,rhs,b);
            if any(fb.*fa > 0)
                disp(['Root finding: initial guesses on one side, # violation: ' num2str(sum(fb.*fa > 0))])
                %disp([a(fb.*fa > 0), b(fb.*fa > 0)])
            end
            c = b - fb.*(b - a)./(fb - fa);  % Regula Falsi
            rf_flag = true;
            for iter = 1:obj.num_Newton
                cold = c;
                [f,df] = eval_int_newton(obj, coef, cdf_poly_base, rhs, cold);
                step = f./df;
                step(isnan(step)) = 0;
                c = cold - step;
                I1 = c<a;
                I2 = c>b;
                I3 = ~I1 & ~I2;
                c  = a.*I1 + b.*I2 + c.*I3;
                if ( norm(f, Inf) < obj.tol ) || ( norm(step, Inf) < obj.tol )
                    rf_flag = false;
                    break;
                end
            end
            %norm(f, Inf)
            if rf_flag
                disp('newton does not converge, continue with regula falsi')
                fc = eval_int_search(obj,coef,cdf_poly_base,rhs,c);
                I1 = (fc < 0);
                I2 = (fc > 0);
                I3 = ~I1 & ~I2;
                a  = I1.*c + I2.*a + I3.*a;
                b  = I1.*b + I2.*c + I3.*b;
                c = regula_falsi(obj, coef, cdf_poly_base, rhs, a, b);
            end
        end
    end
    
end