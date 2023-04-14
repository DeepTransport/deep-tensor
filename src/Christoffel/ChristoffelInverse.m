classdef ChristoffelInverse < ChristoffelSampling 
    
    properties
        poly_cdf
        B
    end
    
    methods
        function obj = ChristoffelInverse(poly)
            %
            class_ind = 0;
            poly_names = {'Legendre', 'Chebyshev1st', 'Chebyshev2nd'};
            cdf_names = {'BoundedPolyCDF', 'Chebyshev1stTrigoCDF', 'Chebyshev2ndTrigoCDF'};
            for i = 1:length(poly_names)
                if isa(poly,poly_names{i})
                    class_ind = i;
                    break;
                end
            end
            if class_ind == 0
                error('polynomial type not supported')
            end
            cdf_builder = str2func(cdf_names{class_ind});
            obj.poly_cdf = cdf_builder(poly, 'err_tol', 1E-10, 'num_Newton', 20);
            %
            B = eval_basis(poly, obj.poly_cdf.nodes);
            obj.B = B.^2 + 1E-6; % collocation points
            obj.poly = poly;
        end
        
        function r = sampling(obj, j, n)
            r = invert_cdf(obj.poly_cdf,obj.B(:,j),rand(1,n));
        end
        
        function plot(obj, j)
            xs = linspace(obj.poly.domain(1), obj.poly.domain(2), 1E3);
            A = eval_basis(obj.poly, xs);
            fi = A(:,j).^2.*eval_measure(obj.poly, xs(:));
            plot(xs, fi)
        end
    end
end