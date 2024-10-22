clc
clear all
close all


% Parâmetros do avião
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
theta0 = 4.6*pi/180;

% %  Simulacao por Espaco de estados
%  Matrizes do Espaco de estados
%  Matriz A clássica
A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,...
-0.0380033,-1.0414,0;0,1,0,0];
% %  Matriz A com adição de um elemento de correção
% A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,...
% -0.0380033,-1.0414,0;0,1,tan(theta0),0];
B = [0, 10.4276, 0.291586;0.932639, 0.745253,0.0708916;-0.0242993, ...
-3.1475, -0.0200573;0, 0, 0];
B1 = [0.291586; 0.0708916; -0.0200573; 0];
B2 = [0, 10.4276 ;0.932639, 0.745253; -0.0242993, -3.1475; 0, 0];

C = eye(size(A));
% C = [1 0 0 0];
D = zeros(4,1);
s_A = size(A);

Matriz_Contr = ctrb(A, B);
Matriz_Obser = ctrb(A', C');
n_Contr = rank(Matriz_Contr);
n_Obser = rank(Matriz_Obser);
if n_Contr == s_A(1)
    disp("Sistema é completamente controlável");
else
    disp("Sistema é não controlável");
end

if n_Obser == s_A(1)
    disp("Sistema é completamente observável");
else
    disp("Sistema é não observável");
end

% omegas_n_dutch = 2.45;
% zetas_dutch = [0.4 sqrt(2)/2 0.9];
omega_n = 2.45;
zeta = sqrt(2)/2;
s_dutch = -(omega_n*zeta - j*omega_n*sqrt(1-zeta^2));
s_dutch1 = -(omega_n*zeta + j*omega_n*sqrt(1-zeta^2));
s_roll = -2.795;
s_spiral = -0.0138;
% K = zeros(2,4*3);
p = [s_spiral, s_dutch, s_dutch1, s_roll];
K = place(A, B2, p);
n = 1;
t = 0:0.01:10;
x0 = [3; 0; 0; 0];
% x0 = [0; 0; 0; 0];
u = zeros(size(t));
% u(1:2) = 100;
colors = ['r', 'b', 'y'];
% for omega_n = omegas_n_dutch
%     for zeta = zetas_dutch
%         s_dutch = -(omega_n*zeta + j*omega_n*sqrt(1-zeta^2));
%         s_dutch_1 = -(omega_n*zeta - j*omega_n*sqrt(1-zeta^2));
%         P = [s_dutch s_dutch_1 s_roll s_spiral];
%         [K_1,prec,message] = place(A,B2,P);
%         K(:,(4*(n-1)+1):(4*n)) = K_1;
%         aviao = ss(A-B2*K_1, B1, C, D);
%         y = lsim(aviao, u, t, x0);
%         figure(1);
%         plot(t,y(:,1),colors(n));
% %         legend('zeta=0.4','zeta=0.707','zeta=0.9');
%         hold on;
%         grid on;
%         
%         figure(2)
%         plot(t,180/pi*y(:,2),colors(n));
% %         legend('zeta=0.4','zeta=0.707','zeta=0.9');
%         hold on;
%         grid on;
%         
%         figure(3)
%         plot(t,180/pi*y(:,3),colors(n));
% %         legend('zeta=0.4','zeta=0.707','zeta=0.9');
%         hold on;
%         grid on;
%         
%         figure(4)
%         plot(t,180/pi*y(:,4),colors(n));
% %         legend('zeta=0.4','zeta=0.707','zeta=0.9');
%         hold on;
%         grid on;
%         n = n+1;
%     end
% end

% aviao = ss(A, B1, C, D);
% y = lsim(aviao, u, t, x0);
% figure(1);
% plot(t,y(:,1));
% grid on;
% 
% figure(2)
% plot(t,180/pi*y(:,2));
% grid on;
% 
% figure(3)
% plot(t,180/pi*y(:,3));
% grid on;
% 
% figure(4)
% plot(t,180/pi*y(:,4));
% grid on;