classdef Lagrange1 < Piecewise
    
    properties
        local_mass(2,2) = [2, 1; 1, 2]/6
        local_weights(1,2) = [1, 1]/2
        local_domain(1,2) = [0, 1]
    end
    
    methods
        function obj = Lagrange1(num_elems)
            obj@Piecewise(1, num_elems);
            %
            obj.nodes = obj.grid;
            %
            obj.jac     = obj.elem_size;
            obj.unweighed_mass    = zeros(cardinal(obj));
            unweighed_weights = zeros(cardinal(obj),1);
            for i = 1:obj.num_elems
                ind = [i, i+1];
                obj.unweighed_mass(ind,ind) = obj.unweighed_mass(ind,ind) + obj.local_mass*obj.jac;
                unweighed_weights(ind)  = unweighed_weights(ind)  + obj.local_weights(:)*obj.jac;
            end
            obj.unweighed_mass    = sparse(0.5*(obj.unweighed_mass+obj.unweighed_mass'));
            %
            obj.unweighed_mass_R  = chol(obj.unweighed_mass);
            obj.unweighed_int_W   = reshape(unweighed_weights, 1, []);
            %
            obj.mass = obj.unweighed_mass/(obj.domain(2)-obj.domain(1));
            obj.mass_R  = obj.unweighed_mass_R/sqrt(obj.domain(2)-obj.domain(1));
            obj.int_W   = obj.unweighed_int_W/(obj.domain(2)-obj.domain(1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function bas = eval_basis(obj, x)
            n   = length(x);
            mask_left   = obj.domain(1) > x(:);
            mask_right  = obj.domain(2) < x(:);
            mask_inside = ~(mask_left | mask_right);
            %{
            if sum(mask_left | mask_right)
                disp('warning: points outside of the domain')
            end
            %}
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x-obj.domain(1))./obj.elem_size);
                ind(ind==0) = 1;
                %
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                %left:  ind: 1-local_x
                %right: ind+1: local_x
                %
                coi = [ind(:); ind(:)+1];
                roi = repmat(reshape(find(mask_inside),[],1), 2, 1);
                val = [1-local_x(:); local_x(:)];
                bas = sparse(roi, coi, val, n, cardinal(obj));
            else
                bas = spalloc(n, cardinal(obj),0);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function bas = eval_basis_deri(obj, x)
            n   = length(x);
            mask_left   = obj.domain(1) > x(:);
            mask_right  = obj.domain(2) < x(:);
            mask_inside = ~(mask_left | mask_right);
            if sum(mask_left | mask_right)
                disp('warning: points outside of the domain')
            end
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x-obj.domain(1))./obj.elem_size);
                ind(ind==0) = 1;
                %
                coi = [ind(:); ind(:)+1];
                roi = repmat(reshape(find(mask_inside),[],1), 2, 1);
                val = [-ones(length(ind), 1)./obj.elem_size; ones(length(ind),1)./obj.elem_size];
                bas = sparse(roi, coi, val, n, cardinal(obj));
            else
                bas = spalloc(n, cardinal(obj),0); 
            end
        end
        
    end
end