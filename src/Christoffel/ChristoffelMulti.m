classdef ChristoffelMulti
    
    properties
        poly
        %
        bound
        maxpdf
        factor
    end
    
    methods (Static)
        function [bound, maxpdf, factor] = set_rejection(poly)
            %
            tol = 1E-8;
            % clipping cdfs
            a = max(poly.nodes(1)*2, poly.domain(1));
            b = min(poly.nodes(end)*2, poly.domain(2));
            xs = linspace(a, b, 1E4);
            %
            A = eval_basis(poly, xs);
            fi = A.^2.*eval_measure(poly, xs(:));
            %
            max_fi = max(fi);
            maxpdf = max_fi*1.1;
            factor = zeros(size(maxpdf));
            bound  = zeros(cardinal(poly),2);
            for i = 1:cardinal(poly)
                ileft = find(fi(:,i)>tol, 1, 'first');
                iright = find(fi(:,i)>tol, 1, 'last');
                bound(i,:) = [xs(ileft), xs(iright)];
                % for rejection sampling
                factor(i) = (xs(iright)-xs(ileft))*maxpdf(i);
            end
        end
        
        function z = rejection_sampling_inner(poly, bound, maxpdf, factor, order, n)
            n = ceil(n*factor(order+1));
            u1 = rand(1,n)*(bound(order+1,2)-bound(order+1,1)) + bound(order+1,1);
            u2 = rand(1,n)*maxpdf(order+1);
            A = eval_basis(poly, u1(:));
            f = A(:,order+1).^2.*eval_measure(poly, u1(:));
            z = u1(f'>u2);
        end
    end
    
    methods
        function obj = Christoffel(arg,varargin)
            p = inputParser;
            %
            addRequired (p,'arg');
            parse(p,arg,varargin{:});
            %
            if isa(arg, 'ApproxBases')
                kk = 1;
                names = cell(1,ndims(arg));
                nc = 0;
                for i = 1:ndims(arg)
                    if nc < cardinal(arg.oneds{i})
                        nc = cardinal(arg.oneds{i});
                        kk = i;
                    end
                    names{i} = class(arg.oneds{i});
                end
                if length(unique(names)) > 1
                    obj.maxpdf = cell(1,ndims(arg));
                    obj.factor = cell(1,ndims(arg));
                    obj.bound = cell(1,ndims(arg));
                    obj.poly = arg.oneds;
                    for i = 1:ndims(arg)
                        [obj.bound{i},obj.maxpdf{i},obj.factor{i}] = Christoffel.set_rejection(obj.poly{i});
                    end
                else
                    obj.poly = cell(1);
                    obj.maxpdf = cell(1);
                    obj.factor = cell(1);
                    obj.bound = cell(1);
                    obj.poly{1} = arg.oneds{kk};
                    [obj.bound{1},obj.maxpdf{1},obj.factor{1}] = Christoffel.set_rejection(obj.poly{1});
                end
            elseif isa(arg, 'Spectral')
                obj.poly = cell(1);
                obj.maxpdf = cell(1);
                obj.factor = cell(1);
                obj.bound = cell(1);
                obj.poly{1} = arg;
                [obj.bound{1},obj.maxpdf{1},obj.factor{1}] = Christoffel.set_rejection(obj.poly{1});
            else
                error('inputs need to be a spectral polynomial or ApproxBases of spectral polynomials')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = rejection_sampling(obj, j, order, n)
            if length(obj.poly) == 1
                pi = 1;
            else
                pi = j;
            end
            z = Christoffel.rejection_sampling_inner(obj.poly{pi}, obj.bound{pi}, obj.maxpdf{pi}, obj.factor{pi}, order, n);
            while numel(z) < n
                zadd = Christoffel.rejection_sampling_inner(obj.poly{pi}, obj.bound{pi}, obj.maxpdf{pi}, obj.factor{pi}, order, n);
                z = [z, zadd];
            end
            z = z(1:n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function z = random(obj, I, tau)
            n = cardinal(I);
            d = ndims(I);
            z = zeros(d,tau,n);
            for i = 1:cardinal(I)
                ind = I.array(i,:);
                for j = 1:d
                    z(j,:,i) = rejection_sampling(obj, j, ind(j), tau);
                end
            end
            z = reshape(z,d,[]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function w = eval_weight(obj, multi, z)
            nb = cardinal(multi);
            z = z';
            [m,d] = size(z);
            w = zeros(m,1);
            if length(obj.poly) == 1
                for i = 1:nb
                    ind = multi.array(i,:);
                    wi = ones(m,1);
                    for j = 1:d
                        A = eval_basis(obj.poly{1}, z(:,j));
                        wi = wi.*(A(:,ind(j)+1).^2);
                    end
                    w = w+wi;
                end
                w = nb./w;
            else
                for i = 1:nb
                    ind = multi.array(i,:);
                    wi = ones(m,1);
                    for j = 1:d
                        A = eval_basis(obj.poly{j}, z(:,j));
                        wi = wi.*(A(:,ind(j)+1).^2);
                    end
                    w = w+wi;
                end
                w = nb./w;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot(obj, j, order)
            if length(obj.poly) == 1
                pi = 1;
            else
                pi = j;
            end
            xs = linspace(obj.bound{pi}(order+1,1), obj.bound{pi}(order+1,2), 1E3);
            A = eval_basis(obj.poly{pi}, xs);
            fi = A(:,order+1).^2.*eval_measure(obj.poly{pi}, xs(:));
            plot(xs, fi)
        end
        
        function debug(obj, j, order, n)
            r = rejection_sampling(obj, j, order, n);
            figure
            plot(obj, j, order)
            hold on
            histogram(r, ceil(order/10+1)*50, 'Normalization', 'pdf')
        end
    end
    
end