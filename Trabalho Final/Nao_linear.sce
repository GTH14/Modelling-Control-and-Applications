clc()
clear()
xdel(winsid())

Ix = 24.675886906e6;
Iz = 67.384152706e6;
Ixz = -2.1150760206e6;
g = 9.81;
c = 340;
M_inf = 0.8;
V = M_inf*c;
u0 = V;
rho = 1.3;
Q = rho*V^2/2;
m = 288773;
S = 541.2;
b = 64.4;
S = 541.2;
b = 64.4;
Cy_beta= -0.88;
Cl_beta = -0.277;
Cn_beta = 0.195;
Cl_p = -0.334;
Cn_p = -0.0415;
Cl_r = 0.300;
Cn_r = -0.327;
Cl_del_a = 0.0137;
Cn_del_a = 0.0002;
Cy_del_r = 0.1157;
Cl_del_r = 0.0070;
Cn_del_r = -0.1256;
theta = 4.6*%pi/180;
x0 = [0;0;0;0];
dt = 0.1
tf = 10;
t = 0:dt:tf;
function [u] = f1(t)
    u = 1/10*%pi/180*[0;1];
endfunction
//function [u] = f1(t)
//    u = [0;0];
//endfunction
function F = f(u,x)
    dbeta = atan(x(1)/u0)
    Cy = Cy_beta*dbeta + Cy_del_r*u(2);
    Y = Q*S*Cy;
    Cl = Cl_beta*dbeta + Cl_p*b/(2*u0)*x(2) + Cl_r*b/(2*u0)*x(3) + Cl_del_a*u(1) +  Cl_del_r*u(2);
    L = Q*S*b*Cl;
    Cn = Cn_beta*dbeta + Cn_p*b/(2*u0)*x(2) + Cn_r*b/(2*u0)*x(3) + Cn_del_a*u(1) + Cn_del_r*u(2);
    N = Q*S*b*Cn;
    F = [Y;L;N];
endfunction
function xdot = naolinear(t,x)
    u = f1(t)
    F = f(u,x); Y = F(1); L = F(2); N = F(3);
    xdot(1) = Y/m + g*sin(x(4))*cos(theta) - x(3)*u0;
    xdot(2) = (Iz*L+Ixz*N)/(Ix*Iz - Ixz^2);
    xdot(3) = (Ix*N+Ixz*L)/(Ix*Iz - Ixz^2);
    xdot(4) = x(2) + x(3)*cos(x(4))*tan(theta);
endfunction

x_nlinear = ode(x0, t(1), t, naolinear);
Titulo_grafico = ['Velocidade lateral', 'Velocidade angular na direção x', 'Velocidade angular na direção z', 'Angulo de rolagem phi']
colors = ["-k", "-k", "-k", "-k"]
for i = 1:size(x_nlinear)(1)
    scf(i)
    plot(t,x(i,:),colors(i), 'LineWidth', 2)
    title(Titulo_grafico(i),'fontsize',4)
    xgrid()
end

