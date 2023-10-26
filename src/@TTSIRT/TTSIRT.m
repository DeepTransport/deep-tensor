classdef TTSIRT < SIRT
    % TTSIRT class
    %
    % TTSIRT Properties:
    %   ys,ms - Cell arrays for computing the marginal density.
    %
    % TTSIRT Methods:
    %   marginalise - 
    %           Marginalises the approximation.
    %   eval_pdf  - 
    %           Evaluates the normalised (marginal) pdf.
    %   eval_irt  - 
    %           X = R^{-1}(Z), where X is the target random variable, R is
    %           the Rosenblatt transport, and Z is the uniform random variable.
    %           * Can map marginal random variables.
    %   eval_cirt - 
    %           Y|X = R^{-1}(Z, X), where X is given, (X,Y) jointly follow
    %           the target represented by SIRT, Z is uniform.
    %           * This function cannot handle marginal random variables.
    %   eval_rt   - 
    %           Z = R(X), where Z is uniform and X is target.
    %           * Can map marginal random variables.
    %   eval_rt_jac - 
    %           Evaluates the Jacobian of Z = R(X).
    %           * This function cannot handle marginal random variables.
    %   random  - 
    %           Generates random variables
    %   sobol - Generates transformed sobol points
    %   set_defensive - 
    %           Resets the defensive term.
    %
    % See also SIRT, ONED, ONEDCDF, TTOPTION, and TTFun
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example:
    %
    % % Setup a d-dimensional Gausssian density with correlation controlled
    % % by a.
    %   d = 20; a = 0.5;
    %   A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
    %   D = diag([1, a*ones(1,d-1)]);
    %   B = D\A;
    %   Q = B'*B;   % precision matrix
    %   C = inv(Q); % covariance matrix
    %   z = sqrt((2*pi)^d/det(Q)); % normalising constant
    % % potential function, which is the negative logarithm of the joint 
    % % distribution, unnormalised
    %   potential = @(x) 0.5*sum((B*x).^2,1);
    %
    % % Set samples used for initilizing TT and debug
    %   debug_x  = B\randn(d, 1E4);
    %   sample_x = B\randn(d, 1E3);
    %   deb = InputData(sample_x, debug_x);
    %
    % % Build the SIRT using Fourier basis, with an approximation domain 
    % % [-5,5]^d.
    %   pol  = Fourier(20);
    %   dom  = BoundedDomain([-5,5]);
    %   base = ApproxBases(Legendre(40), dom, d);
    %
    % % Cross Options
    %   opt = TTOption('als_tol',1E-4,'max_rank',20);
    %
    %   irt = TTSIRT(potential, base, opt, 'var', deb);
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 1: inverse Rosenblatt and Rosenblatt transports
    %   Z   = rand(d, 1E4);     % uniform random seeds
    %   [X,f] = eval_irt(irt,Z);% run IRT, f is the potential function
    %   U   = eval_rt (irt,X);  % U should be the same as Z
    %   p   = eval_potential(irt,X);  % should be the same as f
    %   fe  = potential(X);     % exact potential
    %   norm(Z - U, 'fro')
    %   norm(f - p, 'fro')
    %   figure; plot(fe, f/z, '.')
    %   figure; plot(C - cov(X')); % error in the estimated covariance
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 2: inverse Rosenblatt and Rosenblatt transports for marginal
    % % random variables.
    % % The marginal distribution, unnormalised
    %   marginal_potential = @(x, ind) 0.5*sum((C(ind,ind)\x).*x,1);
    % % normalising constant of the marginal
    %   marginal_z = @(ind) sqrt(det(C(ind,ind))*(2*pi)^length(ind));
    %
    % % Case 2.1: Marginal for X_1, ..., X_k
    % %     - need to run the marginal function with int_dir=1 (default)
    %   if irt.int_dir ~= 1, irt = marginalise(irt, 1); end
    %   ind = 1:8; % the leftover coordinates after marginalisation
    %   [Xm,f] = eval_irt(irt, Z(ind,:));
    %   p = eval_potential(irt, Xm);
    %   Um = eval_rt(irt, Xm);
    %   norm(Z(ind,:) - Um, 'fro')
    %   norm(f - p, 'fro')
    %   fe = marginal_potential(Xm, ind) - log(marginal_z(ind));
    %   figure; plot(fe , f, '.');
    %   figure; plot(C(ind, ind) - cov(Xm'))
    %
    % % Case 2.2: Marginal for X_{k+1}, ..., X_d
    % %     - need to run the marginal function with int_dir=-1
    %   if irt.int_dir ~= -1, irt = marginalise(irt, -1); end
    %   ind = 13:d; % the leftover coordinates after marginalisation
    %   [Xm,f] = eval_irt(irt, Z(ind,:));
    %   p = eval_potential(irt, Xm);
    %   Um = eval_rt(irt, Xm);
    %   norm(Z(ind,:) - Um, 'fro')
    %   norm(f - p, 'fro')
    %   fe = marginal_potential(Xm, ind) - log(marginal_z(ind));
    %   figure; plot(fe , f, '.');
    %   figure; plot(C(ind, ind) - cov(Xm'))
    %
    %%%%%%%%%%%%%%%%%
    %
    % % Task 3: inverse Rosenblatt transport for conditional random variables.
    % % Case 3.1: X_{>=k} | x_{<k}
    % %     - need to run the marginal function with int_dir=1 (default)
    %   indx = 1:8; % index of variables that will be conditioned on
    %   indy = 9:d; % index of conditional variables
    %   X = B\randn(d,1);
    %   X = X(indx);
    %   % conditional mean
    %   my = C(indy,indx)*(C(indx,indx)\X);
    %   % conditional covariance
    %   Cy = C(indy,indy) - C(indy,indx)*(C(indx,indx)\C(indx,indy));
    %   Zy = Z(indy,:);
    %   if irt.int_dir ~= 1, irt = marginalise(irt, 1); end
    %   [Y,f] = eval_cirt(irt, X, Zy);
    %   fe = potential([repmat(X,1,size(Zy,2)); Y]) - marginal_potential(X,indx) + log(marginal_z(indx)) - log(z);
    %   figure; plot(fe , f, '.');
    %   figure; plot(Cy - cov(Y'))
    %   figure; plot(my - mean(Y,2))
    %
    % % Case 3.2: X_{<=k} | x_{>k}
    % %     - need to run the marginal function with int_dir=-1
    %   indx = 9:d; % index of variables that will be conditioned on
    %   indy = 1:8; % index of conditional variables
    %   X = B\randn(d,1);
    %   X = X(indx);
    %   % conditional mean
    %   my = C(indy,indx)*(C(indx,indx)\X);
    %   % conditional covariance
    %   Cy = C(indy,indy) - C(indy,indx)*(C(indx,indx)\C(indx,indy));
    %   Zy = Z(indy,:);
    %   if irt.int_dir ~= -1, irt = marginalise(irt, -1); end
    %   [Y,f] = eval_cirt(irt, X, Zy);
    %   fe = potential([Y; repmat(X,1,size(Zy,2))]) - marginal_potential(X,indx) + log(marginal_z(indx)) - log(z);
    %   figure; plot(fe , f, '.');
    %   figure; plot(Cy - cov(Y'))
    %   figure; plot(my - mean(Y,2))
    %
    % % Case 3.3: Alternative way for generating X_{>=k} | x_{<k}
    % %     - need to run the marginal function with int_dir=1 (default)
    %   indx = 1:8; % index of variables that will be conditioned on
    %   indy = 9:d; % index of conditional variables
    %   X = B\randn(d,1);
    %   X = X(indx);
    %   % conditional mean
    %   my = C(indy,indx)*(C(indx,indx)\X);
    %   % conditional covariance
    %   Cy = C(indy,indy) - C(indy,indx)*(C(indx,indx)\C(indx,indy));
    %   Zy = Z(indy,:);
    %   if irt.int_dir ~= 1, irt = marginalise(irt, 1); end
    %   Zx = eval_rt (irt, X);
    %   fx = eval_pdf(irt, X);
    %   [R,fxy] = eval_irt(irt, [repmat(Zx, 1, size(Zy,2)); Zy]);
    %   Y = R(indy,:);
    %   f = fxy./fx;
    %
    %   fe = potential([repmat(X,1,size(Zy,2)); Y]) - marginal_potential(X,indx) + log(marginal_z(indx)) - log(z);
    %   figure; plot(fe , f, '.');
    %   figure; plot(Cy - cov(Y'))
    %   figure; plot(my - mean(Y,2))
    %
    %%%%%%%%%%%%%%%%%
    %
    
    properties
        ys
        ms
    end
    
    properties (Constant)
        defaultVar = InputData()
        defaultTau = 1E-8
        defaultOpt = TTOption()
        defaultData = TTData()
        defaultPoly = Lagrangep(2,20)
        defaultDomain = BoundedDomain([-1,1])
    end
    
    methods (Static)
        function T = eval_oned_core_213(oned, core, x)
            rkm = size(core, 1);
            nn  = cardinal(oned);
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_radon(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_213_deri(oned, core, x)
            rkm = size(core, 1);
            nn  = cardinal(oned);
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_radon_deri(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_231(oned, core, x)
            rk  = size(core, 3);
            nn  = cardinal(oned);
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_radon(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        function T = eval_oned_core_231_deri(oned, core, x)
            rk  = size(core, 3);
            nn  = cardinal(oned);
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_radon_deri(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        % [z,ys,ms] = cumint(ftt, dir)
        % Marginalise the pdf represented by ftt dimension by dimension
        
    end
    
    methods

        function obj = TTSIRT(potential, arg, varargin)
            obj@SIRT(potential, arg, varargin{:});
            obj.order = [];
        end

        function approx = approximate(obj, func, arg, opt, var)
            approx = TTFun(func, arg, opt, 'var', var);
            if use_amen(approx)
                approx = round(approx);
            end
        end

        function obj = round(obj, tol)
            obj.approx = round(obj.approx, tol);
        end
    end
    
end