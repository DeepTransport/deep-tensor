classdef ConditionalDensity

    methods (Abstract)
        eval_potential_joint(obj, x)
    end

    methods
        function [mlf,gmlf] = eval_potential_conditional(obj, y, x)
            [mllkd,mlp,gmllkd,gmlp] = eval_potential_joint(obj, [repmat(y,1,size(x,2));x]);
            %
            mlf = mllkd + mlp;
            ind = length(y) + (1:size(x,1));
            gmlf = gmllkd(ind,:)+gmlp(ind,:);
        end
        
        function [mlf,gmlf] = pullback_potential(obj, irt, ry, z)
            if nargout == 1
                mlf = pullback(irt, @(x)eval_potential_joint(obj,x), [repmat(ry,1,size(z,2));z]);
            else
                [mlf,gmlf] = pullback(irt, @(x)eval_potential_joint(obj,x), [repmat(ry,1,size(z,2));z]);
            end
        end
        
        function [mlf,gmlf] = pullback_potential_nuts(obj, irt, ry, z)
            [mlf,gmlf] = pullback_potential(obj, irt, ry, z);
            %
            ind = length(ry) + (1:size(z,1));
            gmlf = gmlf(ind,:);
        end

        function [mllkd,mlp] = pullback_potential_pcn(obj, irt, ry, z)
            mlf = pullback_potential(obj, irt, ry, z);
            %
            mlp = 0.5*sum(z.^2,1);
            mllkd = mlf - mlp;
        end
    end
end