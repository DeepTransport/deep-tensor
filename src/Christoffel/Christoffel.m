classdef Christoffel
    
    properties
        sampler
    end
    
    properties (Constant)
        cache = ChristoffelCache
    end
   
    methods
        function obj = Christoffel(arg,varargin)
            defaultCacheSize = 1E4;
            p = inputParser;
            %
            addRequired (p,'arg');
            addOptional(p,'cache_size', defaultCacheSize, @(x) isnumeric(x) && isscalar(x) && (x>0));
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
                    error('only support homogeneous basis in this class, check ChristoffelMulti for inhomogeneous basis')
                else
                    poly = arg.oneds{kk};
                end
            elseif isa(arg, 'Spectral')
                poly = arg;
            else
                error('inputs need to be a spectral polynomial or ApproxBases of spectral polynomials')
            end
            %

            if isa(poly, 'Chebyshev1st') 
                obj.sampler = ChristoffelInverse(poly);
            else
                obj.sampler = ChristoffelReject(poly);
            end
            
            %
            obj.cache.size = [cardinal(poly), p.Results.cache_size];
            obj.cache.data = zeros(obj.cache.size);
            obj.cache.counters = ones(cardinal(poly),1)*obj.cache.size(2);
            %
            fill_cache(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fill_cache_j(obj, j)
            ind = 1:obj.cache.counters(j);
            obj.cache.data(j,ind) = sampling(obj.sampler, j, obj.cache.counters(j));
            obj.cache.counters(j) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function fill_cache(obj)
            for j = 1:obj.cache.size(1)
                fill_cache_j(obj, j);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function z = sample_cache_j(obj, j, tau)
            tau = ceil(tau);
            if (obj.cache.counters(j) + tau) > obj.cache.size(2)
                fill_cache_j(obj, j);
            end
            z = obj.cache.data(j,obj.cache.counters(j)+(1:tau));
            obj.cache.counters(j) = obj.cache.counters(j) + tau;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function z = random(obj, I, tau)
            n = cardinal(I);
            d = ndims(I);
            z = zeros(d,tau,n);
            for i = 1:n
                ind = I.array(i,:);
                for j = 1:d
                    z(j,:,i) = sample_cache_j(obj, ind(j)+1, tau);
                end
            end
            z = reshape(z, d, []);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function w = eval_weight(obj, I, z)
            w = eval_weight(obj.sampler, I, z);
        end
    end
    
end

%{
classdef Christoffel
   
    properties
        method
        %
        poly
        poly_cdf
        %
        max_pdf
        domain_bound
        rej_factor
        %
        tau
        ref
        B
    end
    
    properties (Constant)
        cache = ChristoffelCache
    end
   
    methods
        function obj = Christoffel(polyin,varargin)
            defaultMethod = 'rejection';
            expectedMethod  = {'rejection','inverse'};
            defaultCDFTol = 1E-8;
            defaultNumNewton = 20;
            defaultCacheSize = 1E3;
            p = inputParser;
            %
            addRequired (p,'polyin');
            addParameter(p,'method', defaultMethod, @(x) any(validatestring(x,expectedMethod)));
            addParameter(p,'cache_size', defaultCacheSize, @(x) isnumeric(x) && isscalar(x) && (x>0));
            addParameter(p,'cdf_tol', defaultCDFTol, @(x) isnumeric(x) && isscalar(x) && (x>0) && (x<1));
            addParameter(p,'num_Newton', defaultNumNewton, @(x) isnumeric(x) && isscalar(x) && (x>0));
            parse(p,polyin,varargin{:});
            %
            class_ind = 0;
            poly_names = {'Legendre', 'Hermite', 'Laguerre'};
            cdf_names = {'BoundedPolyCDF', 'HermiteCDF', 'LaguerreCDF'};
            refs = {UniformReference(), ...
                GaussReference(0,1,UnboundedDomain()), ...
                ExpReference(0,1,SemiUnboundedDomain())};
            for i = 1:length(poly_names)
                if isa(polyin,poly_names{i})
                    class_ind = i;
                    break;
                end
            end
            if class_ind == 0
                error('polynomial type not supported')
            end
            %
            cdf_builder = str2func(cdf_names{class_ind});
            %
            obj.poly = polyin;
            obj.poly_cdf = cell(1,cardinal(obj.poly));
            obj.domain_bound = zeros(cardinal(obj.poly), 2);
            %
            obj.method = p.Results.method;
            %
            poly_cdf = cdf_builder(polyin, 'err_tol', p.Results.cdf_tol, 'num_Newton', p.Results.num_Newton);
            % clipping cdfs
            a = poly_cdf.nodes(1);
            b = poly_cdf.nodes(end);
            xs = linspace(a, b, 1E4);
            A = eval_basis(obj.poly, xs);
            fi = A.^2.*eval_measure(obj.poly, xs(:));
            max_fi = max(fi);
            obj.max_pdf = max_fi*1.1;
            obj.rej_factor = zeros(size(obj.max_pdf));
            for i = 1:cardinal(obj.poly)
                ileft = find(fi(:,i)>poly_cdf.tol, 1, 'first');
                iright = find(fi(:,i)>poly_cdf.tol, 1, 'last');
                sampling_nodes = linspace(xs(ileft), xs(iright), 200*ceil(i/10));
                obj.poly_cdf{i} = poly_cdf;
                obj.poly_cdf{i} = update_sampling_nodes(poly_cdf, sampling_nodes);
                obj.domain_bound(i,:) = [xs(ileft), xs(iright)];
                % for rejection sampling
                obj.rej_factor(i) = (xs(iright)-xs(ileft))*obj.max_pdf(i);
            end
            %
            obj.tau = 1E-6;
            obj.ref = refs{class_ind};
            obj.cache.size = [cardinal(obj.poly), p.Results.cache_size];
            obj.cache.data = zeros(obj.cache.size);
            obj.cache.counters = ones(cardinal(obj.poly),1)*obj.cache.size(2);
            %
            B = eval_basis(obj.poly, poly_cdf.nodes);
            obj.B = B.^2; % collocation points
            %
            fill_cache(obj);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function z = rejection_sampling_inner(obj, j, n)
            u1 = rand(1,n)*(obj.domain_bound(j,2)-obj.domain_bound(j,1)) + obj.domain_bound(j,1);
            u2 = rand(1,n)*obj.max_pdf(j);
            A = eval_basis(obj.poly, u1(:));
            f = A(:,j).^2.*eval_measure(obj.poly, u1(:));
            z = u1(f'>u2);
        end
        
        function z = rejection_sampling(obj, j, n)
            z = rejection_sampling_inner(obj, j, ceil(n*obj.rej_factor(j)));
            while numel(z) < n
                znew = rejection_sampling_inner(obj, j, ceil((n-numel(z))*obj.rej_factor(j)));
                z = [z, znew];
            end
            z = z(1:n);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fill_cache_j(obj, j)
            ind = 1:obj.cache.counters(j);
            switch obj.method
                case{'inverse'}
                    z = rand(1,obj.cache.counters(j));
                    obj.cache.data(j,ind) = invert_cdf(obj.poly_cdf,obj.B(:,j),obj.tau,obj.ref,z);
                case{'rejection'}
                    obj.cache.data(j,ind) = rejection_sampling(obj, j, obj.cache.counters(j));
            end
            obj.cache.counters(j) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function fill_cache(obj)
            for j = 1:cardinal(obj.poly)
                fill_cache_j(obj, j);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function z = sample_cache_j(obj, j, tau)
            tau = ceil(tau);
            if (obj.cache.counters(j) + tau) > obj.cache.size(2)
                fill_cache_j(obj, j);
            end
            z = obj.cache.data(j,obj.cache.counters(j)+(1:tau));
            obj.cache.counters(j) = obj.cache.counters(j) + tau;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function z = random(obj, multi_array, tau)
            [n,d] = size(multi_array);
            z = zeros(d,tau,n);
            for i = 1:n
                ind = multi_array(i,:);
                for j = 1:d
                    z(j,:,i) = sample_cache_j(obj, ind(j)+1, tau);
                end
            end
            z = reshape(z, size(z,1), []);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function w = eval_weight(obj, multi, z)
            nb = cardinal(multi);
            z = z';
            [m,d] = size(z);
            w = zeros(m,1);
            for i = 1:nb
                ind = multi.array(i,:);
                wi = ones(m,1);
                for j = 1:d
                    A = eval_basis(obj.poly, z(:,j));
                    wi = wi.*(A(:,ind(j)+1).^2);
                end
                w = w+wi;
            end
            w = nb./w;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot(obj, j)
            A = eval_basis(obj.poly, obj.poly_cdf{j}.sampling_nodes);
            fi = A(:,j).^2.*eval_measure(obj.poly, obj.poly_cdf{j}.sampling_nodes(:));
            plot(obj.poly_cdf{j}.sampling_nodes, fi)
        end
        
        function debug(obj, j, n)
            switch obj.method
                case{'inverse'}
                    z = rand(1,n);
                    r = invert_cdf(obj.poly_cdf{j},obj.B(:,j),obj.tau,obj.ref,z);
                case {'rejection'}
                    r = rejection_sampling(obj, j, n);
            end
            figure
            plot(obj, j)
            hold on
            histogram(r, ceil(j/10)*50, 'Normalization', 'pdf')
        end
    end
    
end
%}