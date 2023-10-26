function obj = round(obj, thres)
% Round the TT cores. tt = ROUND(tt, tol)
% 
%   tol - the truncation threshold of each SVD relative to the largest 
%         singular value)

rs = rank(obj);
if nargin == 1
    thres = obj.opt.local_tol;
end
% Apply double rounding to get to the starting direction
for ii = 1:2
    obj.data.direction = -obj.data.direction;
    if obj.data.direction > 0
        ind = 1:(ndims(obj.base)-1);
    else
        ind = ndims(obj.base):-1:2;
    end
    %
    % start
    for k = ind
        if obj.data.direction > 0
            if k == 1
                Jx_left = [];
            else
                Jx_left = obj.data.interp_x{k-1};
            end
            [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k+1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                Jx_left, obj.data.cores{k+1}, obj.data.cores{k}, ...
                obj.data.direction, obj.opt.int_method, thres, obj.opt.max_rank);
            rs(k) = size(obj.data.cores{k}, 3);
        else
            if k == ndims(obj.base)
                Jx_right = [];
            else
                Jx_right = obj.data.interp_x{k+1};
            end
            [obj.data.cores{k}, obj.data.interp_x{k}, obj.data.cores{k-1}] = TTFun.build_basis_svd(obj.base.oneds{k},...
                Jx_right, obj.data.cores{k-1}, obj.data.cores{k}, ...
                obj.data.direction, obj.opt.int_method, thres, obj.opt.max_rank);
            rs(k-1) = size(obj.data.cores{k}, 1);
        end
    end
end
if use_amen(obj) 
    obj.data.res_x = [];
    obj.data.res_w = [];
end
disp(' >> rounded TT, ranks:')
disp(['  >> ' num2str(rs)])
end