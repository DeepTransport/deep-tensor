classdef AdaptLagrangeP
    
    properties
        domain
        order
        %
        local LagrangeRef
        ref_pts
        ref_basis
        ref_eval_pts
        ref_left_ind_old
        ref_left_ind_new
        ref_right_ind_old
        ref_right_ind_new
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function obj = AdaptLagrangeP(order)
            if order == 1
                disp('should use AdaptLagrange1')
            end
            %
            obj.domain = [-1,1];
            obj.order = order;
            obj.local = LagrangeRef(order+1);
            % 
            obj.ref_pts = [obj.local.nodes(:); obj.local.nodes(:)+1]/2;
            obj.ref_basis = [];
            if rem(obj.order,2) == 0
                ref_eval_ind = [2:obj.order, (obj.order+3):(2*obj.order+1)];
                obj.ref_eval_pts = obj.ref_pts(ref_eval_ind);
                %
                obj.ref_left_ind_old = [1, obj.order+1, obj.order+2, 2*obj.order+2];
                obj.ref_right_ind_old = [1, obj.order/2+1, obj.order/2+1, obj.order+1];
                obj.ref_left_ind_new = ref_eval_ind;
                obj.ref_right_ind_new = 1:(obj.order-1)*2;
            else
                ref_eval_ind = [2:(obj.order+1), (2:obj.order)+obj.order+1]; % delete the repeated point 
                obj.ref_eval_pts = obj.ref_pts(ref_eval_ind);
                %
                obj.ref_left_ind_old = [1, 2*obj.order+2];
                obj.ref_right_ind_old = [1, obj.order+1];
                obj.ref_left_ind_new = 2:(2*obj.order+1);
                obj.ref_right_ind_new = [1:obj.order, obj.order:(2*obj.order-1)];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [grids,f_at_z,l2_err] = approximate(obj, func, tol, max_iter, num_elems)
            nnodes = obj.order + 1; 
            grids = inspace(obj.domain(1),obj.domain(2),num_elems+1);
            ref_ind = 1:num_elems;
            sizes = (2/num_elems)*ones(1,num_elems);
            % initialize
            f_at_z = zeros(nnodes,num_elems);
            l2_err = zeros(1,num_elems);
            %
            pts = reshape(obj.local.nodes(1:obj.order),[],1).*reshape(sizes,1,[]) + reshape(grids(1:num_elems),1,[]);
            pts = [pts(:); obj.domain(2)];
            tmp_f_at_z = feval(func, pts(:));
            f_at_z(1:obj.order,:) = reshape(tmp_f_at_z(1:end-1), obj.order, []);
            f_at_z(end,:) = tmp_f_at_z(nnodes:obj.order:end);


            for iter = 2:max_iter
                ref_l_at_z = reshape(obj.ref_basis*f_at_z(:,ref_ind), nnodes, []);
                ref_sizes = reshape(repmat(reshape(sizes(ref_ind)*0.5,1,[]),2,1), 1,[]);
                %
                ref_f_at_z = zeros(nnodes*2, numel(ref_ind));
                pts = obj.ref_eval_pts(:).*reshape(sizes(ref_ind),1,[]) + reshape(grids(ref_ind),1,[]);
                tmp_f_at_z = reshape(feval(func, pts(:)), size(pts));
                % convert this into ref_f_at_z 
                ref_f_at_z(1,:) = f_at_z(1,ref_ind);
                ref_f_at_z(ind_ref,:) = tmp_f_at_z(ind_eval,:);
                ref_f_at_z(end,:) = f_at_z(end,ref_ind);
                ref_f_at_z = reshape(ref_f_at_z, nnodes, []);
                % local error
                diff = ref_l_at_z - ref_f_at_z;
                ref_l2_err = sqrt(sum(diff.*(obj.local.M*diff),1).*ref_sizes);
                err_flag = ref_l2_err < tol;
                
            end
        end
        
    end
end