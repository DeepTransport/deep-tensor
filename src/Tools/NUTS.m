function out = NUTS(fun, x, N, M, sigma)

%   fun - Given as either a string or an inline function, returns the minus
%         log target and its gradient
%   x   - Initial state
%   N   - Number of steps
%   M   - Number of steps for running adaptation

[mlp,gmlp] = feval(fun, x);

% initial state
state.x = x; % position
state.mlp = mlp; % minus log target at x
state.gmlp = gmlp; % gradient of minus log target at x
state.v = zeros(size(x)); % momentum

%%%% adapt sigma
if nargin < 4
    M = N;
end
% target acceptance probability
delta = 0.8;
% step size
if (nargin<5)||(isempty(sigma))
    sigma = find_sigma(fun,state);
end
%
mu = log(10) + sigma;
sigmabar = 1;
Hbar = 0;
t0 = 10;
gamma = 0.05;
kappa = 0.75;
%%%%

np = length(x);
out.mlp = zeros(1,N);
out.samples = zeros(np,N);
out.sigma = zeros(1,N);
out.tree = zeros(2,N);

out.samples(:,1) = state.x;
out.mlp(1) = state.mlp;
out.sigma(1) = sigma;
% start MCMC
for i = 2:N
    %%%% one iteration of NUTS
    % new momentum
    state.v = randn(size(state.x));
    state.E = state.mlp + 0.5*norm(state.v(:))^2;
    slice = exp(-state.E)*rand;
    %
    lstate = state;
    rstate = state;
    %
    j = 0;
    n = 1;
    s = true;
    while s
        dir = 2*randi([0,1]) - 1;
        if dir == -1
            [lstate,~,statep,np,sp,mh,nmh] = build_tree(fun,lstate,slice,dir,j,exp(sigma));
        else
            [~,rstate,statep,np,sp,mh,nmh] = build_tree(fun,rstate,slice,dir,j,exp(sigma));
        end
        if sp && rand < (np/n)
            state = statep;
        end
        first  = isturning(rstate.x,lstate.x,lstate.v);
        second = isturning(rstate.x,lstate.x,rstate.v);
        s = first & second & sp;
        j = j + 1;
        n = n + np;
    end
    %%%% end of NUTS iteration
    
    %%%% adapt
    if i <= M
        Hbar = ( (i+t0-1)*Hbar + delta - mh./nmh )/(i+t0);
        sigma = mu - sqrt(i)*Hbar/gamma;
        sigmabar = (i^(-kappa))*sigma + (1-i^(-kappa))*sigmabar;
    else
        sigma = sigmabar;
    end
    %%%% end adapt
    
    if mod(i,10)==0, fprintf('NUTS\t i=%d\t nmh=%d\n', i, nmh), end
    
    %
    out.samples(:,i) = state.x;
    out.mlp(i) = state.mlp;
    out.sigma(i) = sigma;
    out.tree(:,i) = [j,nmh];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lstate,rstate,staten,n,s,mh,nmh] = build_tree(fun,state,slice,dir,j,epsilon)

% state - position and momentum variables
% slice - slice variable
% dir   - direction variable
% j     - iteration variable
% epsilon - step size
%
% lstate,rstate - leftmost and rightmost position and momentum states
% n   - integer in Alg. 3
% s   - indicator variable of stopping criterion
% mh  - acceptance probability in MH
% nmh - number of evaluations of mh

deltaMax = 1000; % allowable deviation from Hamiltonian; Recommended on page 9

if j == 0
    staten = leapfrog(fun, state, dir*epsilon);
    %
    lstate = staten;
    rstate = staten;
    %
    n = log(slice) <= -staten.E;
    s = log(slice) < (deltaMax-staten.E);
    mh  = exp(min(state.E-staten.E, 0));
    nmh = 1;
else
    [lstate,rstate,staten,n,s,mh,nmh] = build_tree(fun,state,slice,dir,j-1,epsilon);
    if s
        if dir == -1
            [lstate,~,statep,np,sp,mhp,nmhp] = build_tree(fun,lstate,slice,dir,j-1,epsilon);
        else
            [~,rstate,statep,np,sp,mhp,nmhp] = build_tree(fun,rstate,slice,dir,j-1,epsilon);
        end
        if rand < (np / (n + np))
            staten = statep;
        end
        mh = mh + mhp;
        nmh = nmh + nmhp;
        n = n + np;
        first  = isturning(rstate.x,lstate.x,lstate.v);
        second = isturning(rstate.x,lstate.x,rstate.v);
        s = first & second & sp;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = isturning(rx, lx, v)

b = (rx(:)-lx(:))'*v(:) >= 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state = leapfrog(fun, state, epsilon)

state.v = state.v - (0.5*epsilon)*state.gmlp;
%
state.x = state.x + epsilon*state.v;
%
[state.mlp,state.gmlp] = feval(fun, state.x);
%
state.v = state.v - (0.5*epsilon)*state.gmlp;
%
state.E = state.mlp + 0.5*norm(state.v(:))^2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sigma = find_sigma(fun, state)

sigma = 0;

% momentum
state.v = randn(size(state.x));
state.E = state.mlp + 0.5*norm(state.v(:))^2;
%
statep = leapfrog(fun,state,exp(sigma));
logratio = state.E - statep.E;

a = 2*(logratio>log(0.5))- 1;

while a*logratio > -a*log(2)
    %epsilon = epsilon * 2^a;
    sigma = sigma + a*log(2);
    statep = leapfrog(fun,state,exp(sigma));
    logratio = state.E - statep.E;
end

end
