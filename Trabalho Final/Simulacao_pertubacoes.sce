clc()
clear()
xdel(winsid())

A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,-0.0380033,-1.0414,0;0,1,0,0];
G = [0.291586,0.,272.;0.0708916,2.68243,-2.50156;-0.0200573,0.0380033,1.0414;0,0,0];
C = eye(A);
D = zeros(4,3);

aviao = syslin('c',A,G,C,D);
x0 = [0;0;0;0];
dt = 0.01;
ti = 0;
tf = 100;
t = ti:dt:tf;
u = zeros([t;t;t]);
//// Entrada: impulso na deflexão do aileron
//u(1,1:2) = ;
//// Entrada: impulso na deflexão do leme
//u(2,1:2) = 1/dt;
// Entrada: degrau na deflexão do aileron
u(1,:) = ones(1,length(t))
//// Entrada: degrau na deflexão do leme
//u(2,:) = -1/2*ones(1,length(t))
x = csim(u,t ,aviao,x0)
//x = %pi/180*x;
Titulo_grafico = ['Velocidade lateral', 'Velocidade angular na direção x', 'Velocidade angular na direção z', 'Angulo de rolagem phi']
for i = 1:size(x)(1)
    scf(i+1)
    plot(t,x(i,:),"-r", 'LineWidth', 2)
    title(Titulo_grafico(i),'fontsize',4)
end
//scf(1)
//plzr(aviao)
