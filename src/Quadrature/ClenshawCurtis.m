classdef ClenshawCurtis
    % ClenshawCurtis class
    %
    % ClenshawCurtis Properties:
    %   max_log_order - Max logarithm or the order with base 2
    %   min_log_order - Min logarithm or the order with base 2
    %   tol           - Quadrature error tolerance
    %   ref_pts       - Reference quadrature points in [-1, 1], nx1, 
    %                   defined at the maximum order
    %   weights       - A matrix containing quadrature weights in [-1, 1],
    %                   precomputed for all orders
    %   all_ind       - Boolean indices of active quadrature points at each
    %                   order
    %   new_ind       - Boolean indices of new quadrature points at each
    %                   order
    %
    % ClenshawCurtis Methods:
    %   quad_rule     - Generates quadrature weights and points
    %   int           - Performs adaptive integration
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example: 
    %   fun = @(x,c) 1./(x.^3-2*x-c); % define function
    %   cc = ClenshawCurtis(); % define the quadrature rule
    %   [y, n] = cc.int(@(x) fun(x,5), 0, 2) % integrate
    %   
    %   % We can also return the points where the function is evaluated at
    %   % and the corresponsing weights
    %   [y, n, x, w] = cc.int(@(x) fun(x,5), 0, 2);
    %
    %   % Vector valued boundaries
    %   [y, n, x, w] = cc.int(@(x) fun(x,5), [0,-20], [2,2]);
    %   % the default order may not be sufficient, refine the quadrature rule
    %   cc = ClenshawCurtis('max_order', 1E4); 
    %   y = cc.int(@(x) fun(x,5), [0,-20], [2,2]);
    %
    % See also QUADCC (without using this class)
    
    properties
        max_log_order
        min_log_order
        tol
        %
        ref_pts
        weights
        all_ind
        new_ind
    end
    
    methods (Static)
        function [x, w] = quad_rule(log_order)
            % Generating Clenshaw-Curtis quadrature rule using the method
            % of Jorg Waldvogel
            %
            N = 2^log_order;
            x = cos((0:N)'*pi/N);
            
            N2 = mod(N,2); 
            u0 = 1/(N^2-1+N2); % Boundary weights of CC
            %
            % Clenshaw-Curtis nodes: k = 0,1,...,N; 
            % vector of weights w_0 = w_n = u0
            %
            % auxiliary vectors 
            L = (0:N-1)'; 
            m = min(L,N-L); 
            r = 2./(1-4*m.^2); 
            %
            w = [ifft(r-u0); u0]; % Clenshaw-Curtis weights
        end
    end
    
    methods
        function obj = ClenshawCurtis(varargin)
            defaultMaxOrder = 256;
            defaultMinOrder = 8;
            defaultTol      = 1E-12;
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && all(x > 0);
            %
            addParameter(p, 'max_order', defaultMaxOrder, validScalarPosNum);
            addParameter(p, 'min_order', defaultMinOrder, validScalarPosNum);
            addParameter(p, 'tol',       defaultTol,      validScalarPosNum);
            
            %
            p.KeepUnmatched = false;
            parse(p, varargin{:});
            %
            obj.min_log_order = ceil(log2(p.Results.min_order));
            obj.max_log_order = max(ceil(log2(p.Results.max_order)), obj.min_log_order+1);
            obj.tol = p.Results.tol;
            
            max_n = 2^obj.max_log_order + 1;
            num_k = obj.max_log_order - obj.min_log_order + 1;
            obj.weights = zeros(max_n, num_k);
            for k = num_k:-1:1
                lo = k+obj.min_log_order-1;
                [x, w] = ClenshawCurtis.quad_rule(lo);
                if k == num_k
                    obj.ref_pts = reshape(x, [], 1);
                end
                ind = 1:(2^(num_k-k)):max_n;
                obj.weights(ind,k) = reshape(w, [], 1);
            end
            %
            obj.all_ind = false(max_n, num_k);
            obj.new_ind = false(max_n, num_k);
            for k = 1:num_k
                %
                ind = 1:(2^(num_k-k)):max_n;
                obj.all_ind(ind,k) = true;
                if k == 1
                    obj.new_ind(:,k) = obj.all_ind(:,k);
                else
                    obj.new_ind(:,k) = xor(obj.all_ind(:,k), obj.all_ind(:,k-1));
                end
            end
            
        end
        
        function [y,fcnt,x,w] = int(obj,fun,a,b)
            % Adaptive numerical integration using Clenshaw-Curtis,
            % integrates the function from a to b.
            %   [y,fcn,x,w] = INT(cc,fun,a,b)
            %
            %   fun - given as either a string or an inline function
            %   a   - left boundary, 1 x m vector
            %   b   - right boundary, 1 x m vector
            %   y   - the result of the integration
            %   fcn - the number of function evaluations
            %   x   - quadrature points used
            %   w   - quadrature weights
            
            m = length(a);
            if m ~= length(b)
                error('Boundary points have mismatch dimensions')
            end
            
            f_cache = zeros(2^obj.max_log_order+1, m);
            
            % change of variable: x in [a, b] and z in [-1, 1]
            %   z = 2 (x-a)/(b-a) - 1
            %   x = z (b-a)/2 + (a+b)/2
            jac = reshape((b-a)/2, 1, []); % dx/dz
            %
            tmp = reshape((b+a)/2, 1, []); % shift
            
            y = 0;
            conv_flag = false;
            for lo = obj.min_log_order:obj.max_log_order
                yp = y;
                %
                k = lo-obj.min_log_order+1;
                %
                x = obj.ref_pts(obj.new_ind(:,k)).*jac + tmp;
                f_cache(obj.new_ind(:,k),:) = feval(fun,x);
                y = sum(f_cache(obj.all_ind(:,k),:).*obj.weights(obj.all_ind(:,k),k),1).*jac;
                
                if k > 1 && norm(y - yp, Inf) < obj.tol
                    conv_flag = true;
                    break;
                end
            end
            
            fcnt = 2^lo+1;
            
            if ~conv_flag
                warning('CC quad does not converge')
            end
            
            if nargout > 2
                x = obj.ref_pts(obj.all_ind(:,k)).*jac + tmp;
                w = obj.weights(obj.all_ind(:,k),k).*jac;
            end
        end
    end
end