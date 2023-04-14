classdef SparseFun < ApproxBases
    
    properties
        opt
        %
        indices MultiIndices
        data
        ls
        x
        A
        w
        y
        %
        n_evals
        err
        %
        importance
    end
    
    methods (Static)
        [x,err,Q,R] = ls_solve(A, y, w)
        [x,err,Q,R] = ls_solve_update(Q,R,A,y,w)
    end
    
    methods
        
        obj = re_fit(obj,y)
        obj = least_square(obj,func)
        obj = call_ApproximationToolbox(obj, func)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function n = cardinal(obj)
            n = cardinal(obj.indices);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = clean(obj)
            obj.ls = [];
            obj.x = [];
            obj.A = [];
            obj.w = [];
            obj.y = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function basis_at_z = eval_basis(obj, I, z)
            % basis_at_z: num(x) x cardinal
            if isa(I, 'MultiIndices')
                I = I.array;
            end
            basis_at_z = ones(size(z,2), size(I,1));
            for k = 1:size(I,2)
                ind = I(:,k);
                bk = eval_basis(obj.oneds{k}, reshape(z(k,:),[],1));
                basis_at_z = basis_at_z.*bk(:,ind+1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function f = eval_reference(obj, z)
            basis_at_z = eval_basis(obj, obj.indices, z);
            f = reshape(basis_at_z*obj.data,[],size(z,2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [g,f] = grad_reference(obj, z)
            basis_at_z = eval_basis(obj, obj.indices, z);
            g = zeros(ndims(obj), size(z,2));
            for k = 1:ndims(obj)
                ind = obj.indices.array(:,k);
                bk = eval_basis(obj.oneds{k}, reshape(z(k,:),[],1));
                dbk = eval_basis_deri(obj.oneds{k}, reshape(z(k,:),[],1));
                dbk = dbk./bk;
                dbasis_at_z = basis_at_z.*dbk(:,ind+1);
                g(k,:) = reshape(dbasis_at_z*obj.data,[],size(z,2));
            end
            f = reshape(basis_at_z*obj.data,[],size(z,2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = SparseFun(func, arg, varargin)
            
            % parsing inputs
            defaultPoly     = Legendre(20);
            defaultDomain   = BoundedDomain([-1,1]);
            %
            defaultOption   = SparseOption();
            defaultDebugger = Debugger();
            
            p = inputParser;
            %
            addRequired(p,'func');
            addRequired(p,'arg');
            addOptional(p,'option',defaultOption);
            addParameter(p,'debug', defaultDebugger);
            addParameter(p,'init_indices',[]);
            %
            p.KeepUnmatched = true;
            parse(p,func,arg,varargin{:});
            %
            opt = p.Results.option;
            deb = p.Results.debug;
            init_indices = p.Results.init_indices;
            %
            if isa(arg, 'ApproxBases')
                arg1 = arg.oneds;
                arg2 = arg.oned_domains;
                d = ndims(arg);
            else
                if isnumeric(arg) && isscalar(arg) && (arg > 0)
                    arg1 = defaultPoly;
                    arg2 = defaultDomain;
                    d = arg;
                else
                    error('dimension should be a positive scalar')
                end
            end
            %
            obj@ApproxBases(arg1,arg2,d);
            %
            if isa(arg, 'SparseFun')
                obj.opt = arg.opt;
                obj.data = arg.data;
                obj.indices = arg.indices;
                obj.err = arg.err;
                obj.A = arg.A;
                obj.w = arg.w;
                obj.x = arg.x;
                obj.y = arg.y;
                obj.importance = arg.importance;
            else
                obj.opt = opt;
                obj.data = [];
                obj.indices = [];
                obj.err = Inf;
                obj.A = [];
                obj.w = [];
                obj.x = [];
                obj.y = [];
                switch obj.opt.weight_rule
                    case {'none'}
                        obj.importance = [];
                    case {'Christoffel'}
                        obj.importance = Christoffel(obj);
                end
            end
            %
            if ~isempty(init_indices)
                obj.data = [];
                obj.indices = init_indices;
                obj.err = [];
                obj.A = [];
                obj.w = [];
                obj.x = [];
                obj.y = [];
            end
            %
            if ~isempty(deb) && ~isevaluated(deb) % first evaluate
                deb = set_functions(deb, func(deb.samples));
            end
            %
            switch obj.opt.weight_rule
                case {'none'}
                    obj = least_square(obj, func, deb);
                case {'Christoffel'}
                    if obj.opt.fast
                        disp('Use low rank update of QR');
                        obj = least_square_fast(obj, func, deb);
                    else
                        obj = least_square_weighted(obj, func, deb);
                    end
            end
        end
        
    end
    
end