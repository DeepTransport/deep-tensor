function fx = eval_block(obj, x, dir)
% Evaluate the FTT for either the first or last k variables.
%   f = EVAL_REFERENCE(tt, x, dir)
%
%   x   - input variables, k x n, in the reference domain
%   dir - direction evlauation, >0: from left to right
%                               <0: from right to left
%   f   - function values at x, m x n

k   = size(x, 1);
nx  = size(x, 2);
%
if dir > 0
    % start from the 1st dimension
    fx  = ones(nx,1);
    % all the intermediate dimensions, except the last dimension
    for j = 1:min(k,ndims(obj.base))
        nj  = size(obj.data.cores{j}, 2);
        rjm = size(obj.data.cores{j}, 1);
        %
        if j < ndims(obj.base) || (size(obj.data.cores{j}, 3) > 1 && size(obj.data.cores{j}, 4) == 1)
            tmp = reshape(permute(obj.data.cores{j}, [2,1,3]), nj, []);
        else
            % collapse the third dimension = 1
            tmp = reshape(permute(reshape(obj.data.cores{j}, rjm, nj, []), [2,1,3]), nj, []);
        end
        T   = reshape(permute(reshape( eval_radon(obj.base.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
        % how to speed up this part?
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        B   = sparse(ii(:), jj(:), fx(:), nx, rjm*nx);
        %
        fx  = B*T;
    end
    %fx  = fx';
else
    % start from the last dimension
    xind = k:-1:1;
    tind = ndims(obj.base):-1:1;
    fx   = ones(1,nx);
    % all the intermediate dimensions, except the first dimension
    % need to walk through d-k+1 dimensions
    for i = 1:min(k,ndims(obj.base))
        j = tind(i);
        nj  = size(obj.data.cores{j}, 2);
        rj  = size(obj.data.cores{j}, 3);
        %
        if j > 1 || (size(obj.data.cores{j}, 1) > 1 && size(obj.data.cores{j}, 4) == 1)
            tmp = reshape(permute(obj.data.cores{j}, [2,3,1]), nj, []);
        else
            % collapse the first dimension = 1
            tmp = reshape(obj.data.cores{j}, nj, []);
        end
        T   = reshape(permute(reshape( eval_radon(obj.base.oneds{j},tmp,x(xind(i),:)), nx, rj, []), [2,1,3]), rj*nx, []);
        % how to speed up this part?
        ii  = reshape(1:rj*nx, [], 1);
        jj  = reshape(repmat(1:nx, rj, 1), [], 1);
        B   = sparse(ii, jj, fx(:), rj*nx, nx);
        %
        fx  = T'*B;
    end
end
end