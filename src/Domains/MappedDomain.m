classdef MappedDomain < Domain
    % forward mapping: from infinite (-\infty, \infty) to finite [-1,1]
    
    properties
        scale
    end
    
    methods 
        function obj = MappedDomain(varargin)
            defaultScale = 1;
            p = inputParser;
            addOptional(p,'scale',defaultScale);
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            obj.bound = [-inf, inf];
            obj.scale = p.Results.scale;
        end
        
        function debug_deri(obj, d, n)
            xs = linspace(-d,d,n);
            [~,dzdx] = domain2reference(obj, xs);
            [logdzdx,logdzdx2] = domain2reference_log_density(obj, xs);
            figure
            subplot(2,2,1)
            plot(xs, log(dzdx), 'linewidth', 2)
            hold on
            plot(xs, logdzdx)
            title('logdzdx')
            subplot(2,2,2)
            fdx = (logdzdx(2:end) - logdzdx(1:end-1))/(2*d/n);
            plot(xs, logdzdx2, 'linewidth', 2)
            hold on
            plot((xs(1:end-1)+xs(2:end))/2, fdx)
            title('logdzdx2')
            norm(log(dzdx) - logdzdx)
            
            zs = linspace(-1+1E-2,1-1E-2,n);
            [~,dxdz] = reference2domain(obj, zs);
            [logdxdz,logdxdz2] = reference2domain_log_density(obj, zs);
            subplot(2,2,3)
            plot(zs, log(dxdz), 'linewidth', 2)
            hold on
            plot(zs, logdxdz)
            title('logdxdz')
            subplot(2,2,4)
            fdz = (logdxdz(2:end) - logdxdz(1:end-1))/(2/n);
            plot(zs, logdxdz2, 'linewidth', 2)
            hold on
            plot((zs(1:end-1)+zs(2:end))/2, fdz)
            title('logdxdz2')
            norm(log(dxdz) - logdxdz)
        end
        
        function debug(obj, d, n)
            xs = linspace(-d,d,n);
            [z, dzdx] = domain2reference(obj, xs);
            figure
            subplot(2,2,1)
            plot(xs, z)
            title('z(x)')
            subplot(2,2,2)
            fdzfdx = (z(2:end) - z(1:end-1))/(2*d/n);
            plot(xs, dzdx, 'linewidth', 2)
            hold on
            plot((xs(1:end-1)+xs(2:end))/2, fdzfdx)
            title('dxdz')
            norm(xs - reference2domain(obj,z))
            
            zs = linspace(-1+1E-2,1-1E-2,n);
            [x, dxdz] = reference2domain(obj, zs);
            subplot(2,2,3)
            plot(zs, x)
            title('x(z)')
            subplot(2,2,4)
            fdxfdz = (x(2:end) - x(1:end-1))/(2/n);
            plot(zs, dxdz, 'linewidth', 2)
            hold on
            plot((zs(1:end-1)+zs(2:end))/2, fdxfdz)
            title('dxdz')
            norm(zs - domain2reference(obj, x))
        end
    end
end