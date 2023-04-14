classdef ChristoffelReject < ChristoffelSampling 
    
    properties
        bound
        maxpdf
        factor
    end
    
    methods
        function obj = ChristoffelReject(poly)
            %
            % clipping cdfs
            tol = 1E-8;
            obj.bound = zeros(cardinal(poly), 2);
            a = max(poly.nodes(1)*2, poly.domain(1));
            b = min(poly.nodes(end)*2, poly.domain(2));
            xs = linspace(a, b, 1E4);
            A = eval_basis(poly, xs);
            fi = A.^2.*eval_measure(poly, xs(:));
            max_fi = max(fi);
            obj.maxpdf = max_fi*1.1;
            obj.factor = zeros(size(obj.maxpdf));
            for i = 1:cardinal(poly)
                ileft = find(fi(:,i)>tol, 1, 'first');
                iright = find(fi(:,i)>tol, 1, 'last');
                obj.bound(i,:) = [xs(ileft), xs(iright)];
                % for rejection sampling
                obj.factor(i) = (xs(iright)-xs(ileft))*obj.maxpdf(i);
            end
            %
            obj.poly = poly;
        end
        
        function z = rejection_sampling_inner(obj, j, n)
            u1 = rand(1,n)*(obj.bound(j,2)-obj.bound(j,1)) + obj.bound(j,1);
            u2 = rand(1,n)*obj.maxpdf(j);
            A = eval_basis(obj.poly, u1(:));
            f = A(:,j).^2.*eval_measure(obj.poly, u1(:));
            z = u1(f'>u2);
        end
        
        function z = sampling(obj, j, n)
            z = rejection_sampling_inner(obj, j, ceil(n*obj.factor(j)));
            while numel(z) < n
                znew = rejection_sampling_inner(obj, j, ceil((n-numel(z))*obj.factor(j)));
                z = [z, znew];
            end
            z = z(1:n);
        end
        
        function plot(obj, j)
            xs = linspace(obj.bound(j,1), obj.bound(j,2), 1E3);
            A = eval_basis(obj.poly, xs);
            fi = A(:,j).^2.*eval_measure(obj.poly, xs(:));
            plot(xs, fi)
        end
    end
end