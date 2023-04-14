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
    obj.direction = -obj.direction;
    if obj.direction > 0
        ind = 1:(ndims(obj)-1);
    else
        ind = ndims(obj):-1:2;
    end
    %
    % start
    for k = ind
        if obj.direction > 0
            if k == 1
                Jx_left = [];
            else
                Jx_left = obj.interp_x{k-1};
            end
            [obj.cores{k}, obj.interp_x{k}, obj.cores{k+1}] = FTT.build_basis_svd(obj.oneds{k},...
                Jx_left, obj.cores{k+1}, obj.cores{k}, ...
                obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
            rs(k) = size(obj.cores{k}, 3);
        else
            if k == ndims(obj)
                Jx_right = [];
            else
                Jx_right = obj.interp_x{k+1};
            end
            [obj.cores{k}, obj.interp_x{k}, obj.cores{k-1}] = FTT.build_basis_svd(obj.oneds{k},...
                Jx_right, obj.cores{k-1}, obj.cores{k}, ...
                obj.direction, obj.opt.int_method, thres, obj.opt.max_rank);
            rs(k-1) = size(obj.cores{k}, 1);
        end
    end
end
if strcmp (obj.opt.tt_method, 'amen')
    obj.res_w = [];
end

disp(' >> rounded TT, ranks:')
disp(['  >> ' num2str(rs)])
end