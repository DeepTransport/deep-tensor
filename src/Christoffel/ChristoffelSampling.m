classdef ChristoffelSampling 
    
    properties
        poly
    end
    
    methods (Abstract)
        z = sampling(obj, j, n)
        plot(obj, j)
    end
    
    methods
        function w = eval_weight(obj, I, z)
            nb = cardinal(I);
            z = z';
            [m,d] = size(z);
            w = zeros(m,1);
            for i = 1:nb
                ind = I.array(i,:);
                wi = ones(m,1);
                for j = 1:d
                    A = eval_basis(obj.poly, z(:,j));
                    wi = wi.*(A(:,ind(j)+1).^2);
                end
                w = w+wi;
            end
            w = nb./w;
        end
        
        function debug(obj, j, n)
            r = sampling(obj, j, n);
            figure
            plot(obj, j)
            hold on
            histogram(r, ceil((j-1)/10+1)*50, 'Normalization', 'pdf')
        end
    end
end