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
//A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,-0.0380033,-1.0414,0;0,1,0,0];
//// Matriz A com adição de um elemento de correção
A = [-0.291586,0.,-272.,9.7784;-0.0708916,-2.68243,2.50156,0;0.0200573,-0.0380033,-1.0414,0;0,1,tan(theta0),0];
B = [0, 10.4276, 0.291586;0.932639, 0.745253,0.0708916;-0.0242993, -3.1475, -0.0200573;0, 0, 0]

C = eye(A);
D = zeros(4,3);
// Cria o Espaco de Estados
aviao = syslin('c',A,B,C,D);
// Condicoes iniciais
x0 = [1;0;0;0]
// Cria o vetor de tempos
dt = 0.1; // Discretizacao do tempo
ti = 0; // Instante inicial
tf = 50; // Instante final
t = ti:dt:tf; // vetor de tempos

//Inicializa o vetor de entradas e pertubacoes
u = zeros([t;t;t]);
//// Entrada: impulso na deflexao do aileron [rad]
//u(1,1) = 5*%pi/180; 
//// Entrada: impulso na deflexao do leme [rad]
//u(2,1) = 5*%pi/180;

// Instante inicial e final do pulso
tpi = 1; tpf = 6;
//// Entrada: pulso na deflexao do aileron [rad]
//u(1,tpi/dt:tpf/dt) = 1*%pi/180*ones(1,length(tpi/dt:tpf/dt))
//// Entrada: pulso na deflexao do leme [rad]
//u(2,tpi/dt:tpf/dt) = 1*%pi/180*ones(1,length(tpi/dt:tpf/dt))

//// Entrada: degrau na deflexao do aileron [rad]
//u(1,:) = 1*%pi/180*ones(1,length(t))
//// Entrada: degrau na deflexao do leme [rad]
//u(2,:) = 1/10*%pi/180*ones(1,length(t))

omega = 1 // Frequencia angular da funcao senoidal
//// Entrada: degrau na deflexao do aileron [rad]
//u(1,:) = %pi/180*sin(omega*t)
//// Entrada: degrau na deflexao do leme [rad]
//u(2,:) = %pi/180*sin(omega*t)

// Rajada de vento
A_g = 10 // Metade do pico da rajada[m/s]
omega_g = 2*%pi/10 // Frequencia angular do cosseno
tgi = 0; tgf = 10; // Duração da rajada
tg = tgi:dt:tgf; // Vetor de tempos da rajada
ng = length(tg) // Tamanho do vetor de tempos
// Atribuicao para a entrada da rajada
u(3,(tgi/dt+1):(tgf/dt+1)) = A_g*(1-cos(omega_g*tg))

// Simulacao em EE do sistema
x = csim(u,t ,aviao,x0)

// Geracao dos graficos das respostas do sistema
Titulo_grafico = ['Resposta da velocidade lateral', 'Reposta da velocidade angular na direção x',...
'Resposta da velocidade angular na direção z', 'Resposta do ângulo de rolagem']
colors = ["-r", "-b", "-g", "-m"]
y_label = ["Velocidade lateral (m/s)","Velocidade angular p (rad/s)","Velocidade angular r (rad/s)",...
 "Ângulo de rolagem (rad)"]
name = ["v", "p", "r", "phi"]
for i = 1:size(x)(1)
    scf(i)
    plot(t,x(i,:),colors(i), 'LineWidth', 2)
    title(Titulo_grafico(i),'fontsize',4)
    ylabel(y_label(i))
    xlabel("Tempo t(s)")
    xgrid()
//    xs2png(i, "Resp_rajada_"+name(i))
end

//// Simulacao por Matriz de Transicao
I = eye(A) // Matriz identidade
delta_t = 0.1; // Intervalo de tempo (s)
nt = length(t)-1; // Numero de intervalos de tempo
n = 4  // Numero de termos da expansao de Taylor
Phi = zeros(A) // Inicializa a Matriz de Transicao
[T,lambda] = spec(A) // Calcula
// Cálculo das matrizes pela Expansao de Taylor
for i = 0:(n-1)
    Phi = Phi+A^i*delta_t^i/factorial(i)
end
Gamma = inv(A)*(Phi-I)

// Simulacao com a Matriz de Transição
t_transicao = 0:dt:nt*delta_t// Vetor de tempos para a simulacao com a Matriz de Transicao
x_transicao = zeros(4,length(t_transicao)) // Inicializa o vetor de estados em funcao do tempo
x_transicao(:,1) = x0 // Condicao incial
for i = 2:nt+1 //Simulacao
   Conv = Gamma*B*u(:,i) //Termo de convolucao
   x_transicao(:,i) = Phi*x(:,i-1) + Conv 
end

//// Geracao dos gráficos das respostas obtidas pela Matriz de Transicao
//for i = 1:size(x_transicao)(1)
//    scf(i)
//    plot(t_transicao,x_transicao(i,:),"k", 'LineWidth', 1)
//    legend(["Simulação em EE", "Simulação pela matriz de transição"],4)
////    xs2png(i, "Comparacao_matriz_transicao_"+name(i))
//end

//// Geracao dos Diagramas de Bode do Sistema
G = ss2tf(aviao)
for i = 1:4
    for j = 1:3
        scf(10*j+i)
        bode(G(i,j),1e-2,10/(2*%pi),"rad")
    end
end


//// Simulacao com equacoes não lineares
// Entrada: nula
function [u] = f1(t)
    u = [0;0];
endfunction
// Funcao que calcula a forca e momentos do movimento latero-direcional
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
// Funcao com as equacoes nao lineares
function xdot = naolinear(t,x)
    u = f1(t)
    F = f(u,x); Y = F(1); L = F(2); N = F(3);
    xdot(1) = Y/m + g*sin(x(4))*cos(theta0) - x(3)*u0;
    xdot(2) = (Iz*L+Ixz*N)/(Ix*Iz - Ixz^2);
    xdot(3) = (Ix*N+Ixz*L)/(Ix*Iz - Ixz^2);
    xdot(4) = x(2) + x(3)*cos(x(4))*tan(theta0);
endfunction 
// Simulacao do modelo nao linear
x_nlinear = ode(x0, t(1), t, naolinear);

///Geracao dos gráficos da resposta do modelo nao linear
//colors_n = ["-k", "-k", "-k", "-k"]
//for i = 1:size(x_nlinear)(1)
//    scf(i)
//    plot(t,x_nlinear(i,:),colors_n(i), 'LineWidth', 1.)
//    title(Titulo_grafico(i),'fontsize',4)
//    xgrid()
//    legend(["Linear", "Não Linear"])
////    xs2png(i, "Comparação_não_linear_corrigida_"+name(i))
//end

//// Geracao do Diagrama de Polos
scf(0)
plzr(aviao)

//// Geracao do gráfico da rajada de vento
//scf(5)
//plot(t,u(3,:), 'LineWidth', 2)
//title("Rajada de vento em função do tempo",'fontsize',4)
//xlabel("Tempo t(s)")
//ylabel("Velocidade lateral da rajada (m/s)")
//xgrid()
////xs2png(5, "Modelo rajada de vento")
