A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,-0.0380033,-1.0414,0;0,1,0,0];
B1 = [0.291586; 0.0708916; -0.0200573; 0];
B2 = [0, 10.4276 ;0.932639, 0.745253; -0.0242993, -3.1475; 0, 0];
omega_n = 2.45;
zeta = 0.6;
s_dutch = -(omega_n*zeta - j*omega_n*sqrt(1-zeta^2));
s_dutch1 = -(omega_n*zeta + j*omega_n*sqrt(1-zeta^2));
s_roll = -2.795;
s_spiral = -0.0138;
p = [s_spiral, s_dutch, s_dutch1, s_roll];
K = place(A, B2, p)
F = A-B2*K;
Q= eye(size(A)); P = ones(size(B2));
N = zeros(size(B2));
K = lq(A, B2, Q, P);