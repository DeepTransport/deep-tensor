classdef OU

    properties
        B
        Q
        C
        a
        d
        norm
    end

    methods
        function data = OU(d, a)
            A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
            D = diag([1, a*ones(1,d-1)]);
            %
            data.B = D\A;
            data.Q = data.B'*data.B;
            data.C = inv(data.Q);
            %
            data.norm = sqrt((2*pi)^d / det(data.Q));
            data.d = d;
            data.a = a;
        end

        function f = eval_potential_marginal(data, ind, x)
            C = data.C(ind, ind);
            f = 0.5*sum((C\x).*x,1);
            z = 0.5*log(det(C)) + 0.5*length(ind)*log(2*pi);
            f = z + f;
        end

        function f = eval_potential(data, x)
            global neval__
            f = 0.5*sum((data.B*x).^2,1) + log(data.norm);
            neval__ = neval__ + numel(f);
        end
        
        function [f,g] = eval_potential_nuts(data,x)
            f = 0.5*sum((data.B*x).^2,1) + log(data.norm);
            g = data.Q*x;
        end

        function [f1,f2] = eval_potential_pcn(data,x)
            f1 = 0.5*sum((data.B*x).^2,1);
            f2 = 0.5*sum(x.^2,1);
            f1 = f1 - f2;
        end

        function f = eval_potential_conditional(data, x1, x2, dir)
            %order of the samples (x1, x2)
            x = [x1; x2];
            f = 0.5*sum((data.B*x).^2,1) + log(data.norm);

            if dir > 0
                %from the last dimension, marginalise to the first
                %conditional (x2 | x1)
                ind = 1:size(x1,1);
                fm  = eval_potential_marginal(data, ind, x1);
            else
                %from the first dimension, marginalise to the last
                %conditional (x1 | x2)
                ind = size(x1, 1) + ( 1:size(x2,1) );
                fm  = eval_potential_marginal(data, ind, x2);
            end
            f = f - fm;
        end

    end

end