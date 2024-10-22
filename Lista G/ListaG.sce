clear
clc
xdel(winsid())
//DEFINIÇÃO DE PARAMENTROS
M=200;
J=512;
la=0.8;
lb=0.8;
ka=10000;
kb=10000;
ba=200;
bb=200;
vh=10;
td=(la+lb)/vh;
//MATRIZES
A=[0,0,1,-la;0,0,1,lb;ka/M,-kb/M,-(ba+bb)/M,(-bb*lb+ba*la)/M;ka*la/J,-kb*lb/J,(-bb*lb+ba*la)/J,-(bb*lb^2+ba*la^2)/J];
B=[-1 0;0 -1;ba/M bb/M;-ba*la/J bb*lb/J];
C=[0,0,1,0;0,0,0,1];
D=[0,0;0,0];
x0=[0 0 0 0];
//Entrada Degrau
dt=0.1;
u=ones(2,10/dt);
for t=0:dt:td-dt,
u(2,(t/dt)+1)=0;
end
t=0:0.1:10-dt;
x0=[0 0 0 0]';
sistema=syslin('c',A, B, C, D );
y=csim(u,t,sistema,x0);
vg = y(1,:)
omega = y(2,:)
//GRAFICOS ENTRADA DEGRAU
scf(1)
xset('window',1)
xtitle('Resposta ao degrau', 'tempo', 'Vg (m/s)')
plot(t,vg, 2)
xgrid
scf(2)
xtitle('Resposta ao degrau', 'tempo', 'omega (rad/s)')
plot(t,omega, 2)
xgrid
//Entrada Senoidal

t_sen=0:0.1:25;
u_sen1 = zeros([t_sen;t_sen])
u_sen2  = zeros([t_sen;t_sen])
u_sen1 (1,:) = sin(9.8995*t_sen);
u_sen1(2,:) = sin(9.8995*t_sen);
u_sen2(1,:) = sin(4.9875*t_sen);
u_sen2(2,:) = -sin(4.9875*t_sen);
[y_sen1,x_sen1]=csim(u_sen1,t_sen,sistema,x0);
[y_sen2,x_sen2]=csim(u_sen2,t_sen,sistema,x0);
vg_sen1 = y_sen1(1,:)
omega_sen1 = y_sen1(2,:)
vg_sen2 = y_sen2(1,:)
omega_sen2 = y_sen2(2,:)
GRAFICOS ENTRADA SENOIDAL
EM FASE
scf(3)
xtitle('Resposta a entrada senoidal', 'tempo', 'Vg (m/s)')
plot(t_sen,vg_sen1, 2)
scf(4)
xtitle('Resposta a entrada senoidal', 'tempo', 'omega (rad/s)')
plot(t_sen,omega_sen1, 2)
//EM OPOSIÇÃO DE FASE
scf(5)
xtitle('Resposta a entrada senoidal', 'tempo', 'Vg (m/s)')
plot(t_sen,vg_sen2, 2)
scf(6)
xtitle('Resposta a entrada senoidal', 'tempo', 'omega (rad/s)')
plot(t_sen,omega_sen2, 2)

H = ss2tf(sistema)
//print('H',H)
//bode(H)
for i = 1:2
    for j = 1:2
        scf(10*j+i)
        bode(H(i,j),"rad")
    end
end
