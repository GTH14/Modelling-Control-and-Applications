clc()
clear()
xdel(winsid())

// Matrizes jacobianas do EE
A = [0, 1;-100, 0];
B = [0;10];
C = [1,0];
D = [0];

I = eye(A) // Matriz identidade
dt = 0.2; // Intervalo de tempo (s)
nt = 20; // Número de intervalos de tempo
n = 4  // Número de termos da expansão de Taylor
Phi = zeros(A) // Inicializa a Matriz de Transição
Gamma = zeros(A) // Inicializa a Matriz Gamma
// Cálculo das matrizes pela Expansão de Taylor
for i = 0:(n-1)
    Phi = Phi+A^i*dt^i/factorial(i)
    Gamma = Gamma+A^i*dt^(i+1)/factorial(i+1)
end
//Gamma = dt*Gamma
// Cálculo do termo de convolução para u = 1(t)
Conv = Gamma*B*1

// Simulação com a Matriz de Transição
t = 0:dt:nt*dt// Vetor de tempos para a simulação com a Matriz de Transição
x0 = [0;0] // Vetor de estados inicial
x = zeros(2,length(t)) // Inicializa o vetor de estados em função do tempo
x(:,1) = x0
for i = 2:nt+1 //Simulação
   x(:,i) = Phi*x(:,i-1)+Conv
end

// Simulação utilizando o csim com maior discretização
// praticamente igual a resposta exata do sistema
t1 = 0:0.01:nt*dt // Vetor de tempos para a simulação "Exata"
sistema = syslin('c',A,B,C,D); // Cria o sistema linear do EE
x1 = csim("step",t1,sistema,x0) // Simulação

// Gera o gráfico de x1
scf(1)
plot(t, x(1,:), "-b", t1, x1, "-r", 'LineWidth', 2.5)
xlabel('t [s]','fontsize',4)
ylabel('x1','fontsize',4)
title('Resposta do sistema para a variável x1','fontsize',4)
legend (["Aproximação pela Matriz de Transição","Solução Exata"],3)
xgrid();
