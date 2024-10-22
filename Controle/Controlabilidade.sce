clc()
clear()
xdel(winsid())

// Parâmetros do avião
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
theta0 = 4.6*%pi/180;

//// Simulacao por Espaco de estados
// Matrizes do Espaco de estados
// Matriz A clássica
A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,...
-0.0380033,-1.0414,0;0,1,0,0];
//// Matriz A com adição de um elemento de correção
//A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,...
//-0.0380033,-1.0414,0;0,1,tan(theta0),0];
B = [0, 10.4276, 0.291586;0.932639, 0.745253,0.0708916;-0.0242993, ...
-3.1475, -0.0200573;0, 0, 0]
B1 = [0.291586; 0.0708916; -0.0200573; 0];
B2 = [0, 10.4276 ;0.932639, 0.745253; -0.0242993, -3.1475; 0, 0];

C = eye(A);
//C = [1 0 0 0];
D = zeros(4,1);
Matriz_Contr = cont_mat(A, B2)
n_Contr_1 = rank(Matriz_Contr)
n_Contr = contr(A, B2)
Matriz_Obser = cont_mat(A', C')
n_Obser_1 = rank(Matriz_Obser)
n_Obser = contr(A', C')

if n_Contr == size(A,1) then
    disp("Sistema é completamente controlável")
else
    disp("Sistema não é controlável")
end
if n_Obser == size(A,1) then
    disp("Sistema é completamente observável")
else
    disp("Sistema não é observável")
end
//
//p = poly(A, 's')
//coef_p = coeff(p);
//l_c = length(coef_p)
//W = zeros(l_c-1,l_c-1)
//for i = 1:(l_c-1)
//    for j = i:(l_c-1)
//        W(i,j+1-i) = coef_p(j+1)
//    end
//end
//T = Matriz_Contr*W;

omegas_n_dutch = 2.45;
zetas_dutch = [0.4 sqrt(2)/2 0.9];
s_roll = -2.79;
s_spiral = -0.0137;
//K = zeros(2,4*3);
n = 1;
t = 0:0.01:300;
x0 = [3; 0.0; 0.0; 0.0];
u = zeros(t);
//u(1:2) = 100;
colors = ['r', 'b', 'y'];
for omega_n = omegas_n_dutch
    for zeta = zetas_dutch
        s_dutch = -(omega_n*zeta + %i*omega_n*sqrt(1-zeta^2));
        s_dutch_1 = -(omega_n*zeta - %i*omega_n*sqrt(1-zeta^2));
        disp(s_dutch)
        P = [s_roll s_dutch_1 s_dutch s_spiral];
        K_1 = ppol(A,B2,P);
        K(:,(4*(n-1)+1):(4*n)) = K_1;
        aviao = syslin('c',A-B2*K_1, B1, C, D);
        y = csim( u, t, aviao, x0);
        scf(1);
        plot(t,y(1,:));
//        legend('zeta=0.4','zeta=0.707','zeta=0.9');
        xgrid();
        scf(2)
        plot(t,180/%pi*y(2,:));
//        legend('zeta=0.4','zeta=0.707','zeta=0.9');
        xgrid();
        scf(3)
        plot(t,180/%pi*y(3,:));
//        legend('zeta=0.4','zeta=0.707','zeta=0.9');
        xgrid();
        scf(4)
        plot(t,180/%pi*y(4,:));
//        legend('zeta=0.4','zeta=0.707','zeta=0.9');
        xgrid();
        n = n+1;
    end
end
