classdef ApproxBases
    
    properties
        oneds
        oned_domains
    end
    
    methods
        [f,g] = eval_reference(obj, z)
        g = grad_reference(obj, z)
    end
    
    methods
        function obj = ApproxBases(polys, in_domains, d)
            switch nargin
                case {0,1}
                    error('at least need to specify bases and domain')
                case {2}
                    nps = length(polys);
                    nbs = length(in_domains);
                    d = max(nps,nbs);
                    if nps > 1 && nbs > 1
                        if nps ~= nbs
                            error('dimension of the domain and the number of bases mismatch')
                        end % passing this means d == nps == nbs
                    else
                        if d > 1
                            if nbs == 1
                                in_domains = repmat({in_domains}, d, 1);
                            end
                            if nps == 1
                                polys = repmat({polys}, d, 1);
                            end
                        end
                    end
                case {3}
                    nps = length(polys);
                    nbs = length(in_domains);
                    if nps > 1 && nbs > 1
                        if nps ~= nbs
                            error('dimension of the domain and the number of bases mismatch')
                        end % passing this means nps == nbs
                        if nps ~= d
                            warning('dimesion does not match the number of polynomial bases');
                            d = nps;
                        end % now d == nps == nbs
                    else
                        % determining d
                        if max(nps,nbs) > 1 && max(nps,nbs) ~= d
                            warning('dimesion does not match the number of polynomial bases');
                            d = max(nps,nbs);
                        end
                        %
                        if nbs == 1
                            in_domains = repmat({in_domains}, d, 1);
                        end
                        if nps == 1
                            polys = repmat({polys}, d, 1);
                        end
                    end
            end
            obj.oneds = polys(:);
            obj.oned_domains = in_domains;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [z, mlogf] = sample_measure_reference(obj,n)
            z = zeros(ndims(obj),n);
            mlogf = zeros(1,n);
            for k = 1:ndims(obj)
                z(k,:) = sample_measure(obj.oneds{k},n);
                mlogf = mlogf - eval_log_measure(obj.oneds{k},z(k,:));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x, mlogfx] = sample_measure(obj,n)
            [z, mlogfz] = sample_measure_reference(obj,n);
            [x, dxdz] = reference2domain(obj,z);
            mlogfx = mlogfz + sum(log(dxdz),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function n = cardinals(obj,k)
            if nargin==1
                n = cellfun(@cardinal,obj.oneds);
            else
                n = cellfun(@cardinal,obj.oneds(k));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function d = ndims(obj)
            d = length(obj.oneds);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x, dxdz] = reference2domain(obj, z, ind)
            if nargin == 3
                d = numel(ind);
                z = reshape(z, d, []);
            else
                d = ndims(obj);
                ind = 1:d;
            end
            [x,dxdz] = cellfun(@reference2domain, reshape(obj.oned_domains(ind),[],1), ...
                mat2cell(z,ones(1,d)), 'UniformOutput', false);
            %
            x = reshape(cell2mat(x),size(z));
            dxdz = reshape(cell2mat(dxdz),size(z));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [z, dzdx] = domain2reference(obj, x, ind)
            if nargin == 3
                d = numel(ind);
                x = reshape(x, d, []);
            else
                d = ndims(obj);
                ind = 1:d;
            end
            [z,dzdx] = cellfun(@domain2reference, reshape(obj.oned_domains(ind),[],1), ...
                mat2cell(x,ones(1,d)), 'UniformOutput', false);
            %
            z = reshape(cell2mat(z),size(x));
            dzdx = reshape(cell2mat(dzdx),size(x));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [dlogxdz,d2logxdz2] = reference2domain_log_density(obj, z, ind)
            if nargin == 3
                d = numel(ind);
                z = reshape(z, d, []);
            else
                d = ndims(obj);
                ind = 1:d;
            end
            [dlogxdz,d2logxdz2] = cellfun(@reference2domain_log_density, reshape(obj.oned_domains(ind),[],1), ...
                mat2cell(z,ones(1,d)), 'UniformOutput', false);
            %
            dlogxdz = reshape(cell2mat(dlogxdz),size(z));
            d2logxdz2 = reshape(cell2mat(d2logxdz2),size(z));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [dlogzdx,d2logzdx2] = domain2reference_log_density(obj, x, ind)
            if nargin == 3
                d = numel(ind);
                x = reshape(x, d, []);
            else
                d = ndims(obj);
                ind = 1:d;
            end
            [dlogzdx,d2logzdx2] = cellfun(@domain2reference_log_density, reshape(obj.oned_domains(ind),[],1), ...
                mat2cell(x,ones(1,d)), 'UniformOutput', false);
            %
            dlogzdx = reshape(cell2mat(dlogzdx),size(x));
            d2logzdx2 = reshape(cell2mat(d2logzdx2),size(x));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function arg = duplicate_bases(obj)
            arg = ApproxBases(obj.oneds, obj.oned_domains);
            arg.sampling_measure = obj.sampling_measure;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = remove_bases(obj, ind)
            obj.oneds(ind) = [];
            obj.oned_domains(ind,:) = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function y = eval_measure_potential_reference(obj, z, ind)
            if nargin < 3
                d = ndims(obj);
                ind = 1:d;
            else
                d = length(ind);
            end
            w = cellfun(@eval_log_measure, reshape(obj.oneds(ind),[],1), mat2cell(z,ones(1,d)), 'UniformOutput', false);
            % radon nikodym derivative w.r.t. the base measure
            y = - sum(cell2mat(w),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function g = eval_measure_potential_reference_grad(obj, z, ind)
            if nargin < 3
                d = ndims(obj);
                ind = 1:d;
            else
                d = length(ind);
            end
            g = cellfun(@eval_log_measure_deri, reshape(obj.oneds(ind),[],1), mat2cell(z,ones(1,d)), 'UniformOutput', false);
            g = -cell2mat(g);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [fx, gx] = eval_measure_potential(obj, x, varargin)
            % permute the input
            if length(varargin) > 1
                error('too many optional arguments')
            elseif length(varargin) == 1
                ind = varargin{1};
            else
                ind = 1:size(x,1);
            end
            [z, dzdx] = domain2reference(obj, x, ind);
            fz = eval_measure_potential_reference(obj, z, ind);
            fx = fz + sum(log(dzdx),1);
            if nargout == 2
                gz = eval_measure_potential_reference_grad(obj, z, ind);
                gx(pt,:) = gz.*dzdx(:);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval_measure(obj, x)
            fx = eval_measure_potential(obj, x);
            fx = exp(-fx);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function fx = eval(obj, x)
            % Evaluate the approximated function. f = EVAL(approx, x)
            %
            %   x   - input variables, d x n.
            %   f   - function values at x, 1 x n
            %
            
            z = domain2reference(obj, x);
            fx = eval_reference(obj, z);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [gx,fx] = grad(obj, x)
            % Evaluate the gradient of the approximation
            %   g = GRAD(approx, x)
            %
            %   x   - input variables, d x n, in the reference domain
            %   g   - gradient at x, d x n
            
            [z,dzdx] = domain2reference(obj, x);
            [gz, fx] = grad_reference(obj, z);
            gx = gz.*dzdx;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end