clc()
clear()
xdel(winsid())
// Cria o Espaco de Estados
//// Simulacao por Espaco de estados
// Matrizes do Espaco de estados
A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,-0.0380033,-1.0414,0;0,1,0,0];
B = [0, 10.4276, 0.291586;0.932639, 0.745253,0.0708916;-0.0242993, -3.1475, -0.0200573;0, 0, 0]

C = eye(A);
D = zeros(4,3);

aviao = syslin('c', A,B,C,D)
G = ss2tf(aviao)
for i = 1:4
    for j = 1:3
        scf(10*j+i)
        bode(G(i,j),1e-2,10/(2*%pi),"rad")
    end
end
//s = poly(0,'s');
//G_del_a_v = syslin('c',(8.8938 +36.4803*s+6.6094*s^2)/(0.231043 +16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_del_a_p = syslin('c', (s*(4.88501 +1.18241*s+0.932639*s^2))/(0.231043 +16.9018*s+9.42994*s^2+4.01541*s^3+s^4)) 
//G_del_a_r = syslin('c', (0.165903-0.0293406*s-0.10771*s^2-0.0242993*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_del_a_phi = syslin('c', (4.88501 +1.18241*s+0.932639*s^2)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//
//G_del_r_v = syslin('c', (-69.332+2341.59*s+894.952*s^2+10.4276*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_del_r_p = syslin('c', (s*(-58.9422-7.61948*s+0.745253*s^2))/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_del_r_r = syslin('c', (-2.03363-1.88098*s-9.1799*s^2-3.1475*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_del_r_phi = syslin('c', (-58.9422-7.61948*s+0.745253*s^2)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//
//G_vg_v = syslin('c', (0.231043+16.9018*s+6.5414*s^2+0.291586*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_vg_p = syslin('c', (0.+0.023652*s^2+0.0708916*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_vg_r = syslin('c', (0.0564964*s^2-0.0200573*s^3)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
//G_vg_phi = syslin('c', (0. +0.023652*s+0.0708916*s^2)/(0.231043+16.9018*s+9.42994*s^2+4.01541*s^3+s^4))
////aviao = syslin('c',A,B,C,D);
//scf(1)
//bode(G_del_a_v ,"rad")
//scf(2)
//bode(G_del_a_p ,"rad")
//scf(3)
//bode(G_del_a_r ,"rad")
//scf(4)
//bode(G_del_a_phi ,"rad")
//
////scf(5)
////bode(G_del_r_v ,"rad")
////scf(6)
////bode(G_del_r_p ,"rad")
////scf(7)
////bode(G_del_r_r ,"rad")
////scf(8)
////bode(G_del_r_phi ,"rad")
//
////scf(9)
////bode(G_vg_v,"rad")
////scf(10)
////bode(G_vg_p,"rad")
////scf(11)
////bode(G_vg_r,"rad")
////scf(12)
////bode(G_vg_phi,"rad")
