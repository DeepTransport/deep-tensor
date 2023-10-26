% Constructs Q1 FEM discretization of the fruitfly problem and generates
% cosine KLE
function [p,pm,bound,W1g,W1m,spind,Pua,Phi,lambda,Mass, W2] = build_kle_eig(meshlevel, bc_type, nu, corr_length, tol_kle, m0)
%% Initialize Q1 FEM
h = 2^(-4-meshlevel);
n = 2^(4+meshlevel)+1;
% Gridpoints
p = (0:h:1)';
x = repmat(p,1,n);
y = repmat(p',n,1);
p = [x(:)'; y(:)'];

if (strcmpi(bc_type, 'dn'))
    % Dirichlet-Neumann
    % Indices of Dirichlet boundary nodes. Make it separately left & right
    bound = [find(p(1,:)==0)'; find(p(1,:)==1)'];
else
    % Full Dirichlet
    bound = [find(p(1,:)==0)'; find(p(1,:)==1)'; find(p(2,:)==0)'; find(p(2,:)==1)'];
    bound = unique(bound, 'stable');
end

% u is computed at integer points, a is computed at midpoints 0.5:n-1.5
% Prepare a sparse 3D tensor for creating a matrix out of the coeff.
%%%%%%%%%%%%%%%% 1D gradient
W1g = sparse(n^2, n);
i0 = (1:n)';
% i,i-1,i-1
i = i0(2:n);
j = i-1;
k = i-1;
wijk = (-1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i
i = i0(1:n-1);
j = i+1;
k = i;
wijk = (-1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i-1
i = i0(2:n-1);
j = i;
k = i-1;
wijk = (1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i
i = i0(2:n-1);
j = i;
k = i;
wijk = (1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% 1,1,1
i = i0(1);
j = i;
k = i;
wijk = (1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% n,n,n-1
i = i0(n);
j = i;
k = i-1;
wijk = (1/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);

%%%%%%%%%%%%%%% 1D Mass
W1m = sparse(n^2, n);
% i,i-1,i-1
i = i0(2:n);
j = i-1;
k = i-1;
wijk = (h/6)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i
i = i0(1:n-1);
j = i+1;
k = i;
wijk = (h/6)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i-1
i = i0(2:n-1);
j = i;
k = i-1;
wijk = (h/3)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i
i = i0(2:n-1);
j = i;
k = i;
wijk = (h/3)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% 1,1,1
i = i0(1);
j = i;
k = i;
wijk = (h/3)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% n,n,n-1
i = i0(n);
j = i;
k = i-1;
wijk = (h/3)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);

% Unfortunately, we have to return W1g and W1m separately, since Matlab
% is so stupid that performs sparse*sparse somehow with the full memory
% We will assemble W2g = W1g*C*W1m' + W1m*C*W1g';
% However, we need to perform [1324] permute afterwards. Prepare
% indices.
A = W1g*sparse(ones(n,n))*W1g'; % We only need indices, any matrix is ok
[ij1,ij2]=find(A);
j1 = floor((ij1-1)/n)+1;
i1 = ij1-(j1-1)*n;
j2 = floor((ij2-1)/n)+1;
i2 = ij2-(j2-1)*n;
i = i1+(i2-1)*n;
i = int64(i);
j = j1+(j2-1)*n;
j = int64(j);
ij1 = int64(ij1); % Otherwise we will overflow double
ij2 = int64(ij2);
% Now the permute is computed as follows:
% A_real = sparse(i,j,A(ij1+(ij2-1)*n^2),n^2,n^2);
spind = [i, j, ij1+(ij2-int64(1))*int64(n)^2]; % keep them in the same storage

% Projection from u to inner u and A
Pua = speye(n^2);
Pua(bound,:)=[];


%% Try to assemble 2D constructing tensor directly

[ijg,kg,wg] = find(W1g);
[ijm,km,wm] = find(W1m);

assert(all(ijg==ijm), 'nonzero row positions in W1g and W1m differ')
assert(all(kg==km), 'nonzero column positions in W1g and W1m differ')

[ig, jg] = ind2sub([n n], ijg);
% ijg = tt_ind2sub([n n], ijg);
% ig2 = ijg(:,1);
% jg2 = ijg(:,2);
% assert(all(ig==ig2) && all(jg==jg2), 'broken ind2sub')

i2 = ig + (ig'-1)*n;  i2 = int64(i2(:));
j2 = jg + (jg'-1)*n;  j2 = int64(j2(:));
k2 = kg + (kg'-1)*n;  k2 = int64(k2(:));

w2 = wg * wm' + wm * wg';  w2 = w2(:);

ij2 = i2 + (j2 - int64(1))*int64(n^2);
W2 = sparse(ij2, k2, w2, int64(n)^int64(4), n^2);


%% KLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the KLE C0*(k-startpos)^(-nu)*cos(2*pi*f1(k)*x)*cos(2*pi*f2(k)*y),
% where f1(k) = k - s(k)*(s(k)+1)/2, f2(k) = s(k)-f1(k),
% s(k) = floor(-1/2 + sqrt(1/4+2*k)),
% and C0 is chosen such that sum(lambda)=sigma


fprintf('Computing eigenvalue KLE\n');
x = (0.5:n-1.5)'/(n-1);  % Coeff is defined at midpoints
w = ones(n-1,1)/(n-1);
w = kron(w,w);
pmx = repmat(x,1,n-1);
pmx = pmx(:);
pmy = repmat(x',n-1,1);
pmy = pmy(:);
dist2 = ((pmx-pmx').^2 + (pmy-pmy').^2)/(corr_length^2);
if (isinf(nu))
    % Gaussian covariance
    C = exp(-dist2/2); 
else
    % Matern covariance
    C = (2^(1-nu)/gamma(nu)).*((2*nu*dist2).^(nu/2)).*besselk(nu, sqrt(2*nu*dist2));
    C(dist2==0) = 1;
end
WC = diag(sqrt(w))*C*diag(sqrt(w));
WC = (WC+WC')*0.5;
[Psi,lambda] = eig(WC);
lambda = diag(lambda);
[lambda,prm] = sort(lambda, 'descend');
Psi = Psi(:,prm);
Phi = sqrt(1./w).*Psi;

L = find(1-cumsum(lambda)<tol_kle, 1)
lambda = lambda(1:L);

Phi = Phi(:,1:L);
Phi = reshape(Phi, n-1, n-1, L);
Phi(n, :, :) = 0;
Phi(:, n, :) = 0;
Phi = reshape(Phi, n^2, L);

pmx = reshape(pmx, n-1, n-1);
pmx(n,:) = 1+h/2;
pmx(:,n) = (0.5:n-0.5)'/(n-1);
pmy = reshape(pmy, n-1, n-1);
pmy(:,n) = 1+h/2;
pmy(n,:) = (0.5:n-0.5)'/(n-1);
pm = [pmx(:), pmy(:)]';

%% Windowed mass matrices for pressure observations
% In 1D first, since the domain is symmetric
Mass1 = cell(m0, 1);
for i=1:m0
    % our interval is [(i-1)H, (i+1)H], where H=1/(m0+1)
    ind = double((x>(i-1)/(m0+1)) & (x<(i+1)/(m0+1)))  / (2/(m0+1));
    ind(n) = 0;
%     ind = double(((0:(p(1,2)-p(1,1)):1)>=(i-1)/(m0+1))&((0:(p(1,2)-p(1,1)):1)<=(i+1)/(m0+1)))/(0.5/(m0+1));
    Mass1{i,1} = W1m*sparse(ind(:));
    Mass1{i,1} = reshape(Mass1{i,1}, size(W1m,2), []);
    
    % Point evaluations
    ind = double(abs((0:h:1)' - i/(m0+1))<1e-10);
    Mass1{i,1} = spdiags(ind, 0, n, n);
end
% 2D Masses
Mass = cell(m0, m0);
for j=1:m0
    for i=1:m0
        Mass{i,j} = kron(Mass1{j}, Mass1{i}) ; %  / (1-i/(m0+1));  % rescale observations to O(1)
    end
end
end
