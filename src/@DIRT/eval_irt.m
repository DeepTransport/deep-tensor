function [x,mlogf,gz,Juz,Jux] = eval_irt(obj, z, k)
% Evaluate the deep Rosenblatt transport X = T(Z), where Z is the reference
% random variable and X is the target random variable.
%   [x,f,gz,Juz,Jux] = EVAL_IRT(dirt, z, k)
%
%   z     - reference random variables, d x n
%   k     - number of layers in the evaluations
%   x     - random variable drawn form the pdf defined by DIRT
%   mlogf - negative log of the DIRT density
%   gz    - gradient of negative log of the DIRT density
%   Juz   - cell array of the Jacobian of the diagonal map, U = CDF(Z), at
%            each layer, dimension of each: d x n
%   Jux   - cell array of the Jacobian of the SRT (CDF map), U = R(X), at
%            each layer, dimension of each: d x (dn)

if nargin <= 2
    k = num_layers(obj);
else
    k = min(k, num_layers(obj));
end
x = z;

if nargout <= 2
    mlogf = log_joint_pdf(obj.ref, x);
    mlogf = -mlogf;
    for l = k:-1:1
        %
        if obj.basis_r(l) > 0 && obj.basis_r(l) < obj.d
            xr = obj.basis{l}'*x;
            xn = x - obj.basis{l}*xr;
            % evaluate the diagonal transform
            logd = log_joint_pdf(obj.ref, xr);
            ur = eval_cdf(obj.ref, xr);
            % evaluate sirt, ft is the potential (-log)
            [xr, mlogt] = eval_irt(obj.irts{l}, ur);
            % add complementary sample back for the identity part
            x = obj.basis{l}*xr + xn;
        elseif obj.basis_r(l) == obj.d % reparametrize
            %x = obj.basis{l}'*x;
            % evaluate the diagonal transform
            logd = log_joint_pdf(obj.ref, x);
            u = eval_cdf(obj.ref, x);
            % evaluate sirt, ft is the potential (-log)
            [x, mlogt] = eval_irt(obj.irts{l}, u);
            % double check this
            x = obj.basis{l}*x;
        else
            % evaluate the diagonal transform
            logd = log_joint_pdf(obj.ref, x);
            u = eval_cdf(obj.ref, x);
            % evaluate sirt, ft is the potential (-log)
            [x, mlogt] = eval_irt(obj.irts{l}, u);
        end
        % update density
        mlogf = mlogf + mlogt + logd;
    end
else
    [mlogf,gmlogf] = log_joint_pdf(obj.ref, x);
    mlogf = -mlogf;
    gmlogf = -gmlogf;
    %
    Juz = cell(k,1);
    Jux = cell(k,1);
    glogd = cell(k,1);
    gmlogt = cell(k,1);
    %
    for l = k:-1:1
        %
        if obj.basis_r(l) > 0 && obj.basis_r(l) < obj.d
            xr = obj.basis{l}'*x;
            xn = x - obj.basis{l}*xr;
            % evaluate the diagonal transform
            [logd,glogd{l}] = log_joint_pdf(obj.ref, xr);
            [ur,Juz{l}] = eval_cdf(obj.ref, xr);
            % evaluate sirt, ft is the potential (-log)
            [xr,mlogt,gmlogt{l}] = eval_irt(obj.irts{l}, ur);
            % compuate Jacobian
            Jux{l} = eval_rt_jac(obj.irts{l}, xr, ur);
            % add complementary sample back for the identity part
            x = obj.basis{l}*xr + xn;
        elseif obj.basis_r(l) == obj.d % reparametrize
            %x = obj.basis{l}'*x;
            % evaluate the diagonal transform
            [logd,glogd{l}] = log_joint_pdf(obj.ref, x);
            [u,Juz{l}] = eval_cdf(obj.ref, x);
            % evaluate sirt, ft is the potential (-log)
            [x,mlogt,gmlogt{l}] = eval_irt(obj.irts{l}, u);
            % compuate Jacobian
            Jux{l} = eval_rt_jac(obj.irts{l}, x, u);
            % double check this
            x = obj.basis{l}*x;
            gmlogt{l} = obj.basis{l}*gmlogt{l};
        else
            % evaluate the diagonal transform
            [logd,glogd{l}] = log_joint_pdf(obj.ref, x);
            [u,Juz{l}] = eval_cdf(obj.ref, x);
            % evaluate sirt, ft is the potential (-log)
            [x,mlogt,gmlogt{l}] = eval_irt(obj.irts{l}, u);
            % compuate Jacobian
            Jux{l} = eval_rt_jac(obj.irts{l}, x, u);
        end
        % update density
        mlogf = mlogf + mlogt + logd;
    end
    %
    [d,n] = size(z);
    gz = zeros(d,n);
    for l = 1:k
        if obj.basis_r(l) > 0 && obj.basis_r(l) < obj.d
            if l > 1
                gz = gz + obj.basis{l}*gmlogt{l} + obj.basis{l-1}*glogd{l-1};
            else
                gz = gz + obj.basis{l}*gmlogt{l};
            end
            gr = obj.basis{l}'*gz;
            gn = gz - obj.basis{l}*gr;
            for i = 1:n
                ind = (i-1)*obj.basis_r(l) + (1:obj.basis_r(l));
                gr(:,i) = Juz{l}(:,i) .* ( tril(Jux{l}(:,ind))' \ gr(:,i) );
            end
            gz = obj.basis{l}*gr + gn;
            if l == k
                gz = gz + obj.basis{l}*glogd{l} + gmlogf;
            end
        elseif obj.basis_r(l) == obj.d 
            if l > 1
                gz = gz + gmlogt{l} + glogd{l-1};
                %gz = gz + obj.basis{l}*gmlogt{l} + obj.basis{l-1}*glogd{l-1};
            else
                gz = gz + gmlogt{l};
                %gz = gz + obj.basis{l}*gmlogt{l};
            end
            % reparametrize, double check
            gz = obj.basis{l}'*gz;
            for i = 1:n
                ind = (i-1)*d + (1:d);
                gz(:,i) = Juz{l}(:,i) .* ( tril(Jux{l}(:,ind))' \ gz(:,i) );
            end
            %gz = obj.basis{l}*gz;
        else
            if l > 1
                gz = gz + gmlogt{l} + glogd{l-1};
            else
                gz = gz + gmlogt{l};
            end
            for i = 1:n
                ind = (i-1)*d + (1:d);
                gz(:,i) = Juz{l}(:,i) .* ( tril(Jux{l}(:,ind))' \ gz(:,i) );
            end
        end
    end

end

end