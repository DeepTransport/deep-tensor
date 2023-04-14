classdef Lagrangep < Piecewise
    
    properties
        local LagrangeRef
        global2local(:,:) 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = Lagrangep(order, num_elems)
            obj@Piecewise(order, num_elems);
            %
            if order == 1
                disp('should use Lagrange1')
            end
            %
            obj.local = LagrangeRef(obj.order+1);
            % setup global nodes
            num_nodes = obj.num_elems*(cardinal(obj.local)-1)+1;
            obj.nodes = zeros(num_nodes, 1);
            for i = 1:obj.num_elems
                ind = ( 1:cardinal(obj.local) ) + (cardinal(obj.local)-1)*(i-1);
                obj.nodes(ind) = obj.local.nodes*obj.elem_size + obj.grid(i);
            end
            % map the function value y to each local element
            if cardinal(obj.local) > 2
                j = cardinal(obj.local):(cardinal(obj.local)-1):num_nodes;
                obj.global2local = [reshape(1:(num_nodes-1), cardinal(obj.local)-1, obj.num_elems); j]';
            else
                obj.global2local = [1:(num_nodes-1); 2:num_nodes]';
            end
            % setup the weights, mass matrix and its inverse
            obj.jac     = obj.elem_size/(obj.local.domain(2) - obj.local.domain(1));
            obj.unweighed_mass    = zeros(num_nodes);
            unweighed_weights = zeros(num_nodes,1);
            for i = 1:obj.num_elems
                ind = ( 1:cardinal(obj.local) ) + (cardinal(obj.local)-1)*(i-1);
                obj.unweighed_mass(ind,ind) = obj.unweighed_mass(ind,ind) + obj.local.mass*obj.jac;
                unweighed_weights(ind)  = unweighed_weights(ind) + obj.local.weights(:)*obj.jac;
            end
            %
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
        
        function [bas, w] = eval_basis(obj, x)
            tau = eps; % safe guard thershold
            n   = length(x);
            bas = zeros(n,cardinal(obj));
            mask_left   = obj.domain(1) > x(:);
            mask_right  = obj.domain(2) < x(:);
            mask_inside = ~(mask_left | mask_right);
            %{
            if sum(mask_left | mask_right)
                disp('warning: points outside of the domain')
            end
            %}
            %
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x(:)-obj.domain(1))./obj.elem_size);
                ind(ind==0) = 1;
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                % evaluate the barycentric formula
                diff    = repmat(local_x(:), 1, cardinal(obj.local)) - repmat(obj.local.nodes(:)', length(tmp_x), 1);
                % stablise
                diff(abs(diff)<tau) = tau;
                tmp_m   = repmat(obj.local.omega(:)', length(tmp_x), 1) ./ diff;
                lbs     = tmp_m./sum(tmp_m, 2);
                % embed lbs into the global grid
                % obj.global2local(ind,:) are the col indices
                % repmat(find(mask_inside), 1, cardinal(obj.local))  are the row indices
                %
                coi = obj.global2local(ind,:);
                roi = repmat(find(mask_inside), 1, cardinal(obj.local));
                ii  = (coi-1)*n + roi;
                % evaluation of the internal interpolation
                bas(ii(:))  = lbs(:);
            end
            w = ones(size(x(:)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [bas, w] = eval_basis_deri(obj, x)
            n   = length(x);
            bas = zeros(n,cardinal(obj));
            mask_left   = obj.domain(1) > x(:);
            mask_right  = obj.domain(2) < x(:);
            mask_inside = ~(mask_left | mask_right);
            if sum(mask_left | mask_right)
                disp('warning: points outside of the domain')
            end
            if sum(mask_inside) > 0
                tmp_x   = x(mask_inside);
                % find the element indices for each x
                ind = ceil((tmp_x(:)-obj.domain(1))./obj.elem_size);
                ind(ind==0) = 1;
                % map each x into local coordinate
                local_x = (reshape(tmp_x, 1, []) - reshape(obj.grid(ind), 1, []))./obj.elem_size;
                % evaluate the barycentric formula
                diff    = repmat(local_x(:), 1, cardinal(obj.local)) - repmat(obj.local.nodes(:)', length(tmp_x), 1);
                % stablise
                diff(abs(diff)<eps) = eps;
                tmp_m1 = repmat(obj.local.omega(:)', length(tmp_x), 1) ./ diff;
                tmp_m2 = repmat(obj.local.omega(:)', length(tmp_x), 1) ./ (diff.^2);
                %original function
                %lbs     = tmp_m./sum(tmp_m, 2);
                %
                a = 1./sum(tmp_m1, 2);
                b = sum(tmp_m2, 2).*(a.^2);
                lbs = (tmp_m1.*b - tmp_m2.*a)./obj.jac; 
                % embed lbs into the global grid
                % obj.global2local(ind,:) are the col indices
                % repmat(find(mask_inside), 1, cardinal(obj.local))  are the row indices
                %
                coi = obj.global2local(ind,:);
                roi = repmat(find(mask_inside), 1, cardinal(obj.local));
                ii  = (coi-1)*n + roi;
                % evaluation of the internal interpolation
                bas(ii(:))  = lbs(:);
            end
            w = ones(size(x(:)));
        end
        
    end
end