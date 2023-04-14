classdef LagrangeRef
    % LagrangeRef class
    %
    % Define the reference Lagrange basis, in the reference domain [0,1].
    % This function should not be explicitly used.
    %
    % Constructor input:
    %   n           - Number of interpolation points, should be greater
    %                 than or equal to 2
    %
    % LagrangeRef Properties:
    %
    %   order       - Order of the interpolation polynomial n-1
    %   domain      - The reference domain
    %   nodes       - All the interpolation points, nx1
    %   num_nodes   - Number of nodes
    %   omega       - Barycentric weights of the Lagrange polynomial
    %   mass        - Reference mass matrix
    %   weights     - Reference weighting factors (integration of each basis
    %                 function), nx1
    
    properties
        domain(1,2) = [0, 1]
        order
        nodes(:,1)
        omega(:,1)
        mass(:,:)
        weights(:,1)
    end
    
    methods
        function obj = LagrangeRef(n)
            if n < 2
                error('We need more than two points to define Lagrange interpolation')
            end
            obj.nodes       = zeros(n,1);
            obj.nodes(1)    = 0;
            obj.nodes(end)  = 1;
            if n > 2
                order = n-3;
                Jacob = Jacobi11(order); % in [-1,1]
                % to the interval [0, 1]
                obj.nodes(2:n-1) = 0.5*(Jacob.nodes+1);
            end
            % compute the local omega coefficients
            obj.omega = zeros(n,1);
            for j = 1:n
                ind     = true(n,1);
                ind(j)  = false;
                obj.omega(j) = 1./prod( obj.nodes(j) - obj.nodes(ind) );
            end
            % define the mass matrix
            I = eye(n);
            obj.mass = zeros(n);
            for i = 1:n
                for j = 1:n
                    fij = @(x) eval(obj, I(:,i), x).*eval(obj, I(:,j), x);
                    obj.mass(i,j) = integral(fij, obj.domain(1), obj.domain(2));
                end
            end
            % setup the intergration of each basis
            obj.weights = zeros(n,1);
            for i = 1:n
                fi = @(x) eval(obj, I(:,i), x);
                obj.weights(i) = integral(fi, obj.domain(1), obj.domain(2));
            end
        end
        
        function n = cardinal(obj)
            n = length(obj.nodes);
        end
        
        function f = eval(obj, f_at_x, x)
            tau = eps;
            m   = length(x);
            n   = cardinal(obj);
            f   = zeros(m,1);
            %
            outside = obj.nodes(1)-tau >= x | obj.nodes(n)+tau <= x;
            inside  = ~outside;
            if sum(outside)
                disp('warning: points outside of the domain')
                f(outside) = 0;
            end
            %
            tmp_x   = x(inside);
            diff    = repmat(tmp_x(:), 1, n) - repmat(obj.nodes(:)', length(tmp_x), 1);
            %
            % stablise
            diff(abs(diff)<tau) = tau;
            tmp_m   = repmat(obj.omega(:)', length(tmp_x), 1) ./ diff;
            %
            % evaluation of the internal interpolation
            f(inside)   = sum(repmat(f_at_x(:)',length(tmp_x),1).*tmp_m, 2)./sum(tmp_m, 2);
            %
            f = reshape(f, size(x));
        end
        
    end
end