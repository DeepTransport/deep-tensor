function data = setup_ou_process(d, a)

A = diag(-sqrt(1-a^2)*ones(d-1,1), -1) + eye(d);
D = diag([1, a*ones(1,d-1)]);

data.B = D\A;
data.Q = data.B'*data.B;
data.C = inv(data.Q);

data.norm = sqrt((2*pi)^d / det(data.Q));
data.d = d;
data.a = a;

end