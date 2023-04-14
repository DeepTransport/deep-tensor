classdef Recurr < Spectral
    
    properties
        a(:,1) 
        b(:,1) 
        c(:,1) 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = Recurr(order, a, b, c, normalising)
            %
            obj.a = a;
            obj.b = b;
            obj.c = c;
            %
            a = vpa(a);
            b = vpa(b);
            c = vpa(c);
            % assemble J matrix
            % diagonal term
            T0  = -b./a;
            J1  = sqrt( c(2:end)./(a(1:end-1).*a(2:end)));
            J   = diag(T0,0) + diag(J1,1) + diag(J1,-1);
            %
            % off diagonal
            [V,L] = eig(vpa(J));
            [q_pts,ind] = sort(diag(L));
            %
            obj.order = order;
            obj.nodes = double(q_pts(:));
            obj.weights = double(V(1,ind)'.^2);
            obj.normalising = normalising;
            %
            obj = post_construction(obj);
        end
        
        function f = eval_basis(obj, x)
            %
            % alpha = 1, beta = 1
            % change of variable, x in [0, 1], then y = 2*x - 1
            % all the normalising terms and recurrence terms can be pre-computed
            %
            if obj.order == 0
                f   = ones(length(x), 1).*obj.normalising;
                return;
            end
            %
            f = zeros(length(x), obj.order+1);
            f(:,1) = 1;
            f(:,2) = (obj.a(1)*x(:)+obj.b(1)).*f(:,1);
            for j = 2:obj.order
                f(:,j+1) = (obj.a(j)*x(:) + obj.b(j)).*f(:,j) - obj.c(j)*f(:,j-1);
            end
            f = f.*obj.normalising;
            %w = eval_measure(obj, x);
            %f = (f.*w(:)).*obj.normalising;
        end
        
        function df = eval_basis_deri(obj, x)
            %
            if obj.order == 0
                df = zeros(length(x), 1).*obj.normalising;
                return;
            end
            %
            f = zeros(length(x), obj.order);
            f(:,1) = 1;
            f(:,2) = (obj.a(1)*x(:)+obj.b(1)).*f(:,1);
            for j = 2:obj.order
                f(:,j+1) = (obj.a(j)*x(:) + obj.b(j)).*f(:,j) - obj.c(j)*f(:,j-1);
            end
            
            df = zeros(length(x), obj.order+1);
            %df(:,1) = 0;
            %
            %df(:,2) = obj.a(1)*f(:,1) + (obj.a(1)*x(:)+obj.b(1)).*df(:,1);
            df(:,2) = obj.a(1)*f(:,1);
            for j = 2:obj.order
                %df(:,j+1) = j*f(:,j) + x(:).*df(:,j);
                df(:,j+1) = obj.a(j)*f(:,j) + (obj.a(j)*x(:) + obj.b(j)).*df(:,j) - obj.c(j)*df(:,j-1);
            end
            df = df.*obj.normalising;
            %w = eval_measure(obj, x);
            %df = (df.*w(:)).*obj.normalising;
        end
    end
end