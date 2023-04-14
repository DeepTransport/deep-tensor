function [gx,fx] = grad_reference(obj, x)
% Evaluate the gradient of the TT
%   [g,f] = GRAD(tt, x)
%
%   x   - input variables, k x n, in the reference domain
%   g   - gradient at x, k x n
%   f   - function values at x, 1 x n

nx  = size(x, 2);
fxl = cell(ndims(obj),1);
fxr = cell(ndims(obj),1);
    
for j = 1:ndims(obj)
    nj  = size(obj.cores{j}, 2);
    rjm = size(obj.cores{j}, 1);
    %
    tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
    % rjm nx rj
    if j == 1
        fxl{j} = reshape(permute(reshape( eval_radon(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
    else
        T   = reshape(permute(reshape( eval_radon(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
        % how to speed up this part?
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        Bf  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        %
        fxl{j}  = Bf*T; % nx by rj
    end
end
%
for j = ndims(obj):-1:1
    nj  = size(obj.cores{j}, 2);
    rj  = size(obj.cores{j}, 3);
    %
    tmp = reshape(permute(obj.cores{j}, [2,3,1]), nj, []);
    % rj nx rjm
    if j == ndims(obj)
        fxr{j} = reshape(permute(reshape( eval_radon(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [3,1,2]), rjm, []);
    else
        T = reshape(permute(reshape( eval_radon(obj.oneds{j},tmp,x(j,:)), nx, rj, []), [2,1,3]), rj*nx, []);
        % how to speed up this part?
        ii  = reshape(1:rj*nx, [], 1);
        jj  = reshape(repmat(1:nx, rj, 1), [], 1);
        Bf  = sparse(ii, jj, fxr{j+1}(:), rj*nx, nx);
        %
        fxr{j} = T'*Bf; % rjm by nx
    end
end
%
gx = zeros(ndims(obj),nx);
for j = 1:ndims(obj)
    nj  = size(obj.cores{j}, 2);
    rjm = size(obj.cores{j}, 1);
    %
    tmp = reshape(permute(obj.cores{j}, [2,1,3]), nj, []);
    % (rjm nx) by rj
    D = reshape(permute(reshape( eval_radon_deri(obj.oneds{j},tmp,x(j,:)), nx, rjm, []), [2,1,3]), rjm*nx, []);
    if j == 1
        gx(j,:) = sum(D'.*fxr{j+1}, 1);
    elseif j == ndims(obj)
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        % nx by (rjm nx)
        Bl  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        % fxr{j+1}: rj by nx,   Bl*D: nx by rj
        gx(j,:) =(Bl*D)';
    else
        jj  = reshape(reshape(1:rjm*nx, rjm, nx)', [], 1);
        ii  = repmat((1:nx)', 1, rjm);
        % nx by (rjm nx)
        Bl  = sparse(ii(:), jj(:), fxl{j-1}(:), nx, rjm*nx);
        % fxr{j+1}: rj by nx,   Bl*D: nx by rj
        gx(j,:) = sum((Bl*D)'.*fxr{j+1}, 1);
    end
end
fx = fxr{1};

end
