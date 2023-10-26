classdef ConditionalBanana < ConditionalDensity

    properties
        sigma
    end

    methods
        function obj = ConditionalBanana(sigma)
            obj.sigma = sigma;
        end

        function [mllkd, mlp, gmllkd, gmlp] = eval_potential_joint(obj, z)
            y = z(1,:);
            u = z(2:3,:);
            F = log( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
            G = [F-3; F-5];
            %
            mllkd_post = sum((G-y).^2,1)/(2*obj.sigma^2);
            mlp_post = 0.5*sum(u.^2,1);
            %
            mlp = 0.5*sum(z.^2,1);
            mllkd = mllkd_post + mlp_post - mlp;
            %
            if nargout > 2
                tmp = [2*(u(1,:)-1)+400*(u(1,:).^2-u(2,:)).*u(1,:); 200*(u(2,:)-u(1,:).^2)]...
                    ./( (1-u(1,:)).^2 + 100*(u(2,:)-u(1,:).^2).^2 );
                gmllkd = [-(2*F-2*y-8); (2*F-2*y-8).*tmp]/obj.sigma^2;
                gmlp = z;
            end
        end
    end
end