classdef DoubleBanana

    properties
        sigma
        data
    end

    methods
        function obj = DoubleBanana(sigma, data)
            obj.sigma = sigma;
            obj.data = data;
        end

        function [mllkd, mlp, gmllkd, gmlp] = eval_potential_dirt(obj, u)
            F   = log((1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2);
            mllkd = sum((F-obj.data).^2,1)/(2*obj.sigma^2);
            mlp = 0.5*sum(u.^2,1);
            %lpt = -sum((F-data).^2,1)*beta/(2*sigma^2) - sum(u.^2,1)*beta/2;%
            %p   = exp(lpt);
            tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
                ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
            gmllkd = (F*numel(obj.data) - sum(obj.data(:))).*tmp/obj.sigma^2;
            gmlp = u;
        end
    end
end