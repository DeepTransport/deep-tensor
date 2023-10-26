classdef Banana

    properties
        sigma
        a
    end

    methods
        function obj = Banana(sigma, a)
            obj.sigma = sigma;
            obj.a = a;
        end

        function f = eval_potential(obj, X)
            [f1, f2] = eval_potential_dirt(obj, X);
            f = f1+f2;
        end

        function [f1, f2] = eval_potential_dirt(obj, X)
            x = X(1,:);
            y = X(2,:);
            %
            f1 = (y-obj.a*x.^2+2).^2*10;
            f1 = f1/obj.sigma;
            %
            f2 = (y.^2+x.^2)/2;
        end
    end
end