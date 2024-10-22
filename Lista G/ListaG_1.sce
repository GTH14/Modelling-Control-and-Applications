// sistema linear
A=[0,0,1, -la;0,0,1,lb;-ka/M,-kb/M,-(ba + bb)/M,(ba*la + bb*lb)/M;..
ka*la/J,-kb*lb/J,(ba*la - bb*lb)/J,-(ba*(la^2) + bb*(lb^2))/J]
B=[-1,0;0,-1;-ba/M,-bb/M;(ba*la)/J,(-bb*lb)/J]
C=[0,0,1,0;0,0,0,1]
D=[0,0;0,0]
G = syslin('c',A,B,C,D)
Gs = ss2tf(G)
X0 = [0; 0; 0; 0] // condição inicial
t = 0:0.01:4
t2 = 0:0.1:10
td = 1.6/vh;
//entradas
//tipo = 'unitaria' // 'unitaria' ou 'seno'
tipo = 'seno'
u = zeros (2, length(t));
u2 = zeros (2, length(t2));
u3 = zeros (2, length(t2));
if tipo == 'unitaria' then
    u(1,:) = 1 //vc
    for i = 1:length(t)
        if t(i) >= td then // vd
            u(2:1) = 1
        end
    end
elseif tipo == 'seno' then
/// u2 = zeros (2, length(t2));
// u3 = zeros (2, length(t2));
//    i = 1
//    for t2 = 0:0.1:10
    u2(1,:) = sin(9.8995*t2);
    u2(2,:) = sin(9.8995*t2);
    u3(1,:) = sin(4.9875*t2);
    u3(2,:) = -sin(4.9875*t2);
//        i = i + 1
//    end
end
x0 = [0; 0; 0; 0]
//[y , x] = csim (u , t , G , X0)
[y2, x2] = csim (u2, t2, G, x0)
[y3, x3] =csim (u3, t2, G, x0)
if tipo == 'unitaria' then
f1 = scf(1)
plot(t,y(1,:), 'r', 'linewidth',2)
xgrid()
xtitle ('Velocidade do centro de massa', 'tempo [s]', 'deslocamento [m]')
f2 = scf(2)
plot(t,y(2,:), 'b', 'linewidth',2)
xgrid()
xtitle ('Velocidade angular', 'tempo [s]', 'deslocamento [m]')
disp(Gs)
f3 = scf(3)
bode(G(1,1))
f4 = scf(4)
bode(G(1,2))
f5 = scf(5)
bode(G(2,1))
f6 = scf(6)
bode(G(2,2))
end
if tipo == 'seno' then
    f7 = scf(7)
    plot(t,y2(1,:), 'r',t2,y2(2,:),'b', 'linewidth',2)
    xgrid()
    xtitle ('Resposta a entrada senoidal 1', 'tempo [s]', 'deslocamento [m]')
    legend('velocidade do centro de massa','velocidade angular')
    f8 = scf(8)
    plot(t,y2(1,:), 'r',t2,y3(2,:),'b', 'linewidth',2)
    xgrid()
    xtitle ('Resposta a entrada senoidal 1', 'tempo [s]', 'deslocamento [m]')
    legend('velocidade do centro de massa','velocidade angular')
end
