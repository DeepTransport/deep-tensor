classdef EmpiricalMap
    
    properties
        d
        n
        nbins
        domains
        edges
        counts
        hs
        pdfs
        cdfs
        map_bound
        map_tail
    end
    
    methods
        function [u,f] = eval_cdf(obj, x)
            % Map to uniform
            if size(x,1) ~= obj.d
                error('dimension of the input parameter does not match the dimension of the map')
            end
            u = zeros(size(x));
            f = ones(obj.d,size(x,2));
            for i = 1:obj.d
                ind = sum(x(i,:)' > obj.edges(i,:), 2)';
                u(i,:) = (obj.pdfs(i,ind)).*(x(i,:)-obj.edges(i,ind)) + obj.cdfs(i,ind);
                f(i,:) = obj.pdfs(i,ind);
%                 f = f.*obj.pdfs(i,ind);
            end
        end
        
        function [z,dzdx] = eval_map(obj, x)
            [u,dudx] = eval_cdf(obj, x);
            %
            scale = 1 - obj.map_tail*2;
            z = erfinv((2*u-1)*scale)*sqrt(2);
            dudz = exp( -0.5*z.^2-0.5*log(2*pi) )/scale;
            dzdx = dudx./dudz;
        end
                
        function [x,dxdz] = eval_imap(obj, z)
            scale = 1 - obj.map_tail*2;
            u = 0.5 + erf(z/sqrt(2))/(2*scale); 
            dudz = exp( -0.5*z.^2-0.5*log(2*pi) )/scale;
            [x,dudx] = eval_icdf(obj, u);
            dxdz = dudz./dudx;
        end
        
        function [x,f] = eval_icdf(obj, u)
            % Map from uniform
            if size(u,1) ~= obj.d
                error('dimension of the input parameter does not match the dimension of the map')
            end
            x = zeros(size(u));
            f = ones(1,size(x,2));
            for i = 1:obj.d
                ind = sum(u(i,:)' > obj.cdfs(i,:), 2)';
                x(i,:) = (u(i,:) - obj.cdfs(i,ind))./obj.pdfs(i,ind) + obj.edges(i,ind);
                f = f.*obj.pdfs(i,ind);
            end
        end
        
        function f = eval_pdf(obj, x)
            if size(x,1) ~= obj.d
                error('dimension of the input parameter does not match the dimension of the map')
            end
            f = ones(1,size(x,2));
            for i = 1:obj.d
                ind = sum(x(i,:)' > obj.edges(i,:), 2)';
                f = f.*obj.pdfs(i,ind);
            end
        end
        
        function h = plot_cdf(obj, i)
            x = linspace(0, obj.hs(i), 5)';
            xpts = reshape(x + obj.edges(i,1:end-1),[],1);
            ypts = reshape((obj.pdfs(i,:)).*x + obj.cdfs(i,1:end-1),[],1);
            h = plot(xpts, ypts);
        end
        
        function h = plot_pdf(obj, i)
            xs = reshape([obj.edges(i,1:end-1); obj.edges(i,2:end)], [], 1);
            ys = reshape([obj.pdfs(i,:); obj.pdfs(i,:)], [], 1);
            h = plot([obj.edges(i,1);xs;obj.edges(i,end)], [0;ys;0]);
        end
        
        function obj = EmpiricalMap(data,varargin)
            % data: d x n
            
            defaultNbins = -1;
            defaultGBound = 4;
            p = inputParser;
            addRequired(p, 'data',  @(x) isnumeric(x));
            addOptional(p, 'nbins', defaultNbins, @(x) isnumeric(x) && isscalar(x));
            addParameter(p,'gauss_bound', defaultGBound, @(x) isnumeric(x) && isscalar(x) && all(x > 0));
            %
            p.KeepUnmatched = false;
            parse(p,data,varargin{:});
            %
            [obj.d, obj.n] = size(data);
            obj.nbins = p.Results.nbins;
            obj.map_bound = p.Results.gauss_bound;
            obj.map_tail = ( 1+erf(-obj.map_bound/sqrt(2)) )/2;
            %
            ub = max(data, [], 2);
            lb = min(data, [], 2);
            obj.domains = [lb-(ub-lb)*1E-3,  ub+(ub-lb)*1E-3];
            
            if obj.nbins == -1
                N = zeros(obj.d,1);
                for i = 1:obj.d
                    tmp = histcounts(data(i,:));
                    N(i) = length(tmp);
                end
                obj.nbins = max(N);
            end
            obj.edges = zeros(obj.d,obj.nbins+1);
            obj.counts = zeros(obj.d,obj.nbins);
            obj.hs = zeros(obj.d,1);
            obj.pdfs = zeros(obj.d,obj.nbins);
            obj.cdfs = zeros(obj.d,obj.nbins+1);
            for i = 1:obj.d
                obj.edges(i,:) = linspace(obj.domains(i,1), obj.domains(i,2), obj.nbins+1);
                obj.counts(i,:) = histcounts(data(i,:),obj.edges(i,:));
                obj.hs(i) = (obj.domains(i,2)-obj.domains(i,1))/obj.nbins;
                obj.pdfs(i,:) = (obj.counts(i,:)/obj.n/obj.hs(i));
                obj.cdfs(i,:) = cumsum([0,obj.counts(i,:)])/obj.n;
            end
        end
    end
    
    %[N,EDGES] = histcounts(X,M)
    
end