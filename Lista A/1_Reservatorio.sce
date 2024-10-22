// Dados do problema
S = 10; // [m^2]
R = 2e8; // [Pa/(m^3/s)^2]
rho = 1000; // [kg/m^3]
g = 10; // [m/s^2]
Qe = 0.010247; // [m^3/s]
h0 = 1; //[m] 

// Funcao que define a derivada de h
function hdot = f1(h)
    hdot = 1/S*(-sqrt(rho*g*h/R)+Qe);
endfunction

// Método_ de Euler
dt_euler = 10; // Passo de integracao
ti = 0; // Tempo inicial
tf = 30000; // Tempo final
tv_euler = ti:dt_euler:tf;
n_euler = length(tv_euler);
hv_euler = zeros(1,n_euler); // Inicializa o vetor com as alturas
hv_euler(1) = h0;
for i = 1:(n_euler-1)
    hdot = f1(hv_euler(i));
    hv_euler(i+1) = hv_euler(i) + dt_euler*hdot;
end
// Método_ de Runge-Kutta
dt_rk = 10;
tv_rk = ti:dt_rk:tf;
n_rk = length(tv_euler);
hv_rk = zeros(1,n_rk); // Inicializa o vetor com as alturas
hv_rk(1) = h0;
for i = 1:(n_rk-1)
    K1 = f1(hv_rk(i));
    K2 = f1(hv_rk(i)+dt_rk/2*K1);
    K3 = f1(hv_rk(i)+dt_rk/2*K2);
    K4 = f1(hv_rk(i)+dt_rk*K3);
    hv_rk(i+1) = hv_rk(i) + dt_rk*(K1+2*K2+2*K3+K4)/6;
end
scf(0);
plot(tv_euler,hv_euler,'b-', tv_rk, hv_rk,'r-','LineWidth',2.5);
legend([ "Método de  Euler " , "Método de Runge−Kutta " ],4);
xgrid()
