classdef LeastSquare

    properties
        A
        y
        w
        qA
        rA
        %
        x 
        err
        opt_err
    end

    methods
        function obj = LeastSquare()
            obj.A = [];
            obj.y = [];
            obj.w = [];
            obj.qA = [];
            obj.rA = [];
            %
            obj.x = [];
            obj.err = Inf;
            obj.opt_err = Inf;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function f = isempty(obj)
            f = isempty(obj.qA);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = initialize(obj,A,y,w)
            if nargin == 4
                obj.A = A.*sqrt(w);
                obj.y = y.*sqrt(w);
                obj.w = w(:);
            else
                obj.A = A;
                obj.y = y;
                obj.w = ones(size(A,1),1);
            end
            [obj.qA,obj.rA] = qr(obj.A,0);
            obj.x = obj.rA\(obj.qA'*obj.y);
            [obj.err, obj.opt_err] = cross_validate(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x = add_cols_solve(obj, Aadd)
            if ~isempty(obj.w)
                Aadd = Aadd.*sqrt(obj.w);
            end
            R1 = obj.qA'*Aadd;
            [Q2, R2] = qr(Aadd - obj.qA*R1, 0);
            %
            rA = [obj.rA, R1; zeros(size(Aadd,2),size(obj.qA,2)), R2];
            qA = [obj.qA, Q2];
            %
            x = rA\(qA'*obj.y);
            %A = [obj.A, Aadd];
            %
            %norm(qA*rA - A, 'fro')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = add_cols(obj, Aadd)
            Aadd = Aadd.*sqrt(obj.w);
            %
            R1 = obj.qA'*Aadd;
            [Q2, R2] = qr(Aadd - obj.qA*R1, 0);
            %
            obj.rA = [obj.rA, R1; zeros(size(Aadd,2),size(obj.qA,2)), R2];
            obj.qA = [obj.qA, Q2];
            %
            obj.x = obj.rA\(obj.qA'*obj.y);
            obj.A = [obj.A, Aadd];
            [obj.err, obj.opt_err] = cross_validate(obj);
            %
            %norm(obj.qA*obj.rA - obj.A, 'fro')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = add_rows(obj, Aadd, yadd, wadd)
            if isempty(obj)
                obj = initialize(obj, Aadd, yadd, wadd);
            else
                if nargin == 4
                    obj.y = [obj.y; yadd.*sqrt(wadd)];
                    obj.w = [obj.w(:); wadd(:)];
                    Aadd = Aadd.*sqrt(wadd);
                else
                    obj.y = [obj.y; yadd];
                    obj.w = [obj.w(:); ones(size(Aadd,1),1)];
                end
                n = size(obj.qA,2);
                [Q1,R1] = qr([obj.rA;Aadd],0);
                %
                %[m,n] = size(obj.qA);
                %r = size(Aadd,1);
                % obj.qA = [obj.qA, zeros(m,r); zeros(r,n), eye(r)]*Q1(:,1:n);
                obj.qA = [obj.qA*Q1(1:n,1:n); Q1((n+1):end,1:n)];
                obj.rA = R1(1:n,1:n);
                %
                obj.x = obj.rA\(obj.qA'*obj.y);
                obj.A = [obj.A; Aadd];
                [obj.err, obj.opt_err] = cross_validate(obj);
                %
                %norm(obj.qA*obj.rA - obj.A, 'fro')
                %[obj.qA, obj.rA] = qr(obj.A, 0);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [err, opt_err] = cross_validate(obj)
            [m,n] = size(obj.A);
            nout = size(obj.y,2);
            second_moment = sum(obj.y.^2,1)/size(obj.y,1);
            % Compute the relative leave-one-out cross-validation error the fast
            % leave-one-out cross-validation procedure [Cawlet & Talbot 2004] based on
            % the Sherman-Morrison-Woodbury formula
            %
            err = zeros(1,nout);
            if m-1 < n
                warning('Not enough samples for cross-validation');
                err(:) = Inf;
                return
            end
            % Compute the predicted residuals using the Bartlett matrix inversion formula
            % (special case of the Sherman-Morrison-Woodbury formula)
            T = sum(obj.qA.^2,2);
            del = (obj.y-obj.qA*(obj.rA*obj.x))./(1-T);
            err = sum(del.^2,1)/m;
            err = err./second_moment;
            %
            % Compute the relative cross-validation error
            %if rcond(obj.rA) > 1E-16
            dra = diag(obj.rA);
            if min(abs(dra))/max(abs(dra)) > 1E-16
                %invrA = inv(rA);
                %C = invrA*invrA';
                % (1-n/m)^(-1)*(1+trace(C))
                % Direct Eigenvalue Estimator [Chapelle, Vapnik & Bengio, 2002], [Blatman & Sudret, 2011]
                % -> accurate even when m is not >> n (corr ~ 1+2*n/m when m->Inf, as trace(C) ~ n/m when m->Inf)
                %corr2 = (m/(m-n))*(1+trace(C));
                %norm(sum(diag(rA).^(-1))^2 - trace(C))
                %norm(sum(diag(rA).^(-1)) - trace(invrA))
                %norm(corr-corr2)
                corr = (m/(m-n))*(1+sum(dra.^(-1))^2);
                if corr > 0
                    err = err*corr;
                end
            end
            err = sqrt(err);
            %
            opt_err = normest( obj.rA'*(obj.rA/size(obj.A,1)) - eye(size(obj.A,2)), 1e-2);
        end

    end

end