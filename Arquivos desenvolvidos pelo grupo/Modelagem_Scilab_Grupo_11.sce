clear
xdel( winsid() ) // Fecha janelas gráficas

xe = -19 // distancia do CG ate a empenagem em metros (negativo porque é para trás), tirado do desenho em planta
xa = -0.5 //distancia do CG ate a asa em metros (negativo porque é para trás), tirado do desenho em planta
m = 50000 // mmasa da aeronave em kg, tirado do paper
Aa = 92.5 // area de referencia da asa em m², tirado do paper
ARa = 8 // Aspect Ratio da asa, tirado do paper
Ae = 26 // area de referencia da empenagem em m², tirado do paper
ARe = 5.6 // Aspect Ratio da empenagem, calculado no XFLR5
Vcz = 830/3.6 // Velocidade de cruzeiro tirada da wikipedia, convertida de km/h para m/s
h = 10000 // Altura de cruzeiro em metros, estimada com base na altura máxima tirada de diversas fontes
ea = 1 // Fator de correção de Oswald da asa
ee = 1 // Fator de correção de Oswald da empenagem

I = 2157603 // Inércia do avião em kg*m² estimada no CAD
g = 9.8 // aceleração da gravidade ao nível do mar em m/s²
ca = 3.57 // corda média da asa em metros, calculada no XFLR5
ce = 2.24 // corda média da empenagem em metros, calculada no XFLR5
W = m*g // Peso da aeronave em Newtons
//Vmax = 890/3.6 // km/h -> m/s

//Calculando rho do ar, com dados e fórmula tirada do PPT

rho0 = 1.225 // Densidade do ar no nível do mar em kg/m³
lambda = -6.5/1000 // Constante que entra na equação em K/m
R = 287 // Constante universal dos gases em J/kgK
T0 = 288.15 // Temperatura de referência em K
rho = rho0*(1+lambda*h/T0)^(-(g/(lambda*R))-1) // Densidade do ar na altura em cruzeiro em kg/m³

//Thrust maximo (82.3kN tirado da General Electric)
Tmax = 2*82.3*10^3*rho/rho0 // N corrigindo pela densidade na altura de cruzeiro

//Coeficientes aerodinamicos com base nas contas preliminares, foi selecionado um perfil NACA, de onde foram tirados os coeficientes do XFLR5 e do database AirfoilTools.com

//Cl
cldeltae = -0.876
cl0a  = 0.358
cl0e =  -cl0a
clalphaa = 0.088*180/%pi // para ficar em radianos~]
clalphae = clalphaa

//Cd
cd0 = 0.008 // cd0 do perfil de asa
//cd0fus = (Tmax-0.5*rho*cd0*Vmax^2*(Aa+Ae))/(0.5*rho*Vmax*Vmax*Af)
cd0fus = 0.1*%pi*3.35*3.35*0.25/Aa //0.1 foi tirado do forum de aviacao para area frontal
cd0e = cd0 // cd0 empenagem
cd0a = cd0+cd0fus // cd0 da asa mais a fuselagem

//Cm para o centro aerodinamico (1/4 de corda, conforme calculado no XFLR5)
cmaca = -0.05
cmace = 0.05

//Valores tirados do EES
Tcz = 30415.4924 // Thrust de equilíbrio em Newtons
alphacz = 0.0288763294 // angulo de ataque de equilibrio em radianos
deltaequilibrio = -0.157216933 // angulo do elevador na condição de cruzeiro em radianos

//Encontrando as matrizes do sistema linear


//Definindo variaves para economizar escrita

sacz = sin(alphacz) // seno do alpha de cruzeiro
cacz = cos(alphacz) // cosseno do alpha de cruzeiro
PdAa = 0.5*rho*Aa*Vcz^2 // Pressao dinamica * Areaa
DPdAa = rho*Aa*Vcz // Derivada da pressão dinâmica * Aa na velocidade
PdAe = 0.5*rho*Ae*Vcz^2 // Pressao dinamica * Areae
DPdAe = rho*Ae*Vcz // Derivada da pressão dinâmica * Ae na velocidade
DenCda = %pi*ARa*ea // Denominador do Cda
DenCde = %pi*ARe*ee // Denominador do Cde

//Fazendo derivadas parciais de cada termo, bem como forças e momentos, tudo na condição de equilibrio calculada no EES

//Sustentação Asa
cla = clalphaa*alphacz+cl0a
La = PdAa*cla // Sustentação gerada pela asa em newtons
DLaDv = DPdAa*cla // Derivada de La na velocidade
DLaDgama = -PdAa*clalphaa // Derivada de La em gama
DLaDtheta = PdAa*clalphaa // Derivada de La em theta
DLaDdelta = 0

//Sustentação empenagem
cle = +clalphae*alphacz+cl0e+cldeltae*deltaequilibrio
Le = PdAe*cle // Sustentação gerada pela empenagem em newtons
DLeDv = DPdAe*cle // Derivada de Le na velocidade
DLeDgama = -PdAe*clalphae//talvez sinal errado -> Acho que tava certo // Derivada de Le em gama
DLeDtheta = PdAe*clalphae//talvez sinal errado -> Acho que tava certo // Derivada de Le em theta
DLeDdelta = PdAe*cldeltae

//Arrasto Asa
cda = cd0a+(cla^2)/DenCda
Da = PdAa*cda // Arrasto gerado pela asa em newtons
DDaDv = DPdAa*cda // Derivada de Da na velocidade
DDaDgama = PdAa*(clalphaa^2*(-2*alphacz)-2*clalphaa*cl0a)/DenCde // Derivada de Da em gama
DDaDtheta = PdAa*(clalphaa^2*(2*alphacz)+2*clalphaa*cl0a)/DenCda // Derivada de Da em theta
DDaDdelta = 0

//Arrasto empenagem
cde = cd0e+(cle^2)/DenCde
De = PdAe*cde // Arrasto gerado pela empenagem em newtons
DDeDv = DPdAe*cde // Derivada de De na velocidade
DDeDgama = PdAe*(clalphae^2*(-2*alphacz)-2*clalphae*cl0e-2*clalphae*cldeltae*deltaequilibrio)/DenCde//talvez sinal do termo que multilpica dois cls tenha que inverter -> Acho que arrumei// Derivada de Dem em gama
DDeDtheta = PdAe*(clalphae^2*(2*alphacz)+2*clalphae*cl0e+2*clalphae*cldeltae*deltaequilibrio)/DenCde//talvez sinal do termo que multilpica dois cls tenha que inverter -> Acho que arrumei // Derivada de De em theta
DDeDdelta = PdAe*(2*cl0e*cldeltae+2*clalphae*cldeltae*alphacz + 2*(cldeltae^2)*deltaequilibrio)/DenCde

//Momentos
Ma = PdAa*ca*cmaca // Momento gerado pela distribuicao de pressoes na asa em Nm
Me = PdAe*ce*cmace // Momento gerado pela distribuicao de pressoes na empenagem em Nm

//Conferindo se derivadas do equilibrio sao nulas

Vdoteq = (Tcz*cacz-Da-De)/m
Gamadoteq = (Tcz*sacz+La+Le-W)/(m*Vcz)
ThetaDDeq = (La*xa*cacz+Le*xe*cacz+Da*xa*sacz+De*xe*sacz+Ma+Me)/I

//Derivadas parciais das equações

//Velocidade dot
DvdotDv = 1*(-DDaDv-DDeDv)/m
DvdotDgama = 1*(Tcz*sacz-DDaDgama-DDeDgama-W)/m
DvdotDtheta = 1*(-Tcz*sacz-DDaDtheta-DDeDtheta)/m
DvdotDT = cacz/m
DvdotDdelta = 1*(-DDeDdelta)/m

//Gama dot
DgamadotDv = ((1*(DLaDv+DLeDv)*(m*Vcz))-m*(Tcz*sacz+La+Le-W))/((m*Vcz)^2)
DgamadotDgama = 1*(-Tcz*cacz+DLaDgama+DLeDgama)/(m*Vcz)
DgamadotDtheta = 1*(Tcz*cacz+DLaDtheta+DLeDtheta)/(m*Vcz)
DgamadotDT = sacz/(m*Vcz)
DgamadotDdelta = 1*(DLeDdelta)/(m*Vcz)

//Theta dot dot
DthetaddDv = 1*(DLaDv*xa*cacz+DLeDv*xe*cacz+DDaDv*xa*sacz+DDeDv*xe*sacz+DPdAa*ca*cmaca+DPdAe*ce*cmace)/I
DthetaddDgama = 1*(DLaDgama*xa*cacz+La*xa*sacz+DLeDgama*xe*cacz+Le*xe*sacz+DDaDgama*xa*sacz-Da*xa*cacz+DDeDgama*xe*sacz-De*xe*cacz)/I
DthetaddDtheta = 1*(DLaDtheta*xa*cacz-La*xa*sacz+DLeDtheta*xe*cacz-Le*xe*sacz+DDaDtheta*xa*sacz+Da*xa*cacz+DDeDtheta*xe*sacz+De*xe*cacz)/I
DthetaddDT = 0
DthetaddDdelta = 1*(DLeDdelta*xe*cacz+DDeDdelta*xe*sacz)/I

//Montando o espaco de estados, fazendo x1 = V-Vcz, x2 = Gama, x3 = thetadot, x4 = theta-thetaeq, u1 = T - Tcz, u2 = delta - deltaequilibrio

//Modulos das perturbacoes
PertV = Tcz*0.1 // Forca em Newtons
PertGama = Tcz*0.1 // Forca em Newtons
PertThetad = 60000 // Momento em Nm

A = [DvdotDv,DvdotDgama,0,DvdotDtheta;DgamadotDv,DgamadotDgama,0,DgamadotDtheta;DthetaddDv,DthetaddDgama,0,DthetaddDtheta;0,0,1,0]
B = [DvdotDT, DvdotDdelta, PertV/m;DgamadotDT,DgamadotDdelta,PertGama/(m*Vcz);DthetaddDT,DthetaddDdelta,PertThetad/(I);0,0,0]
C = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]
D = [0,0,0;0,0,0;0,0,0;0,0,0]

//Montando EE
aviao=syslin('c',A,B,C,D)

//Diagrama de polos e zeros
//scf()
//plzr(aviao)

tf = 2000
//Simulando
p=1 // passo de integracao
t=0:p:tf // discretizando o vetor de tempo

//Entradas
EntradaT = ones(1,length(t))*0
//EntradaT(1,100:500) = Tmax-Tcz //Pulso
EntradaDelta = ones(1,length(t))*0

//seno
/*for i = 1:length(t)
    EntradaDelta(1,i) = sin(2*%pi*9/1000*t(i))*1*%pi/180
end
*/

EntradaPert = zeros(1,length(t)) // inicializando

//EntradaPert em pulsos
InicioEntradaP = [100,1100] // em segundos
FimEntradaP = [400,1400] // em segundos
SinalEntradaP = [0,0]
for i = 1:length(InicioEntradaP)
    EntradaPert(1,InicioEntradaP(i)/p:FimEntradaP(i)/p) = SinalEntradaP(i)
end

//EntradaPert em seno
//for i = 1:length(t)
//   EntradaPert(1,i) = sin(2*%pi*1/10*t(i))
//end

udegrau=[EntradaT;EntradaDelta;EntradaPert] // Vetor de entradas
V0 = 0 // V-Vcz
Gamma0 = 0*%pi/180 // Gamma
ThetaDot0 = 0 // ThetaDot
Theta0 = 0*%pi/180 // Theta-ThetaEquil (theta-alphacz)

x0 = [V0;Gamma0;ThetaDot0;Theta0] //CIs

// Simulando
[yc,xc]=csim(udegrau,t,aviao,x0)
    Velocidade = xc(1,:) + Vcz
    Gama = xc(2,:)
    Thetadot = xc(3,:)
    Theta = xc(4,:) + alphacz
    Alpha = Theta - Gama

//Calculando posicoes    
Z(1,1)=h
X(1,1)=0
for i = 1:length(t)-1
    Z(:,i+1) = Z(:,i) + p*Velocidade(i)*sin(Gama(i))
    X(:,i+1) = X(:,i) + p*Velocidade(i)*cos(Gama(i))
end
Vnorm = (Velocidade-Vcz)/(Velocidade(length(t))-Vcz)
//Plotando graficos
/*scf(99) // Grafico de velocidade normalizada
plot2d(t,(Velocidade-Vcz)/(Velocidade(length(t))-Vcz),1,rect=[0,0,2000,2])
T=list("Resposta do sistema no tempo","Tempo (s)","Variação da velocidade normalizada em relação à de cruzeiro","V^*")
xtitle(T(1),T(2),T(3))
xgrid()*/

scf(0)// Variaveis do sistema e trajetoria
subplot(221)   
plot2d(t,(Velocidade),1)//,rect=[0,230,200,240])
T=list("Resposta do sistema no tempo","Tempo (s)","Velocidade (m/s)","Velocidade")
xtitle(T(1),T(2),T(3))
xgrid()

subplot(222)
plot2d(t,Gama*180/%pi,2)//,rect=[0,-5,1000,7])
plot2d(t,Theta*180/%pi,3)
plot2d(t,Alpha*180/%pi,4)
T=list("Resposta do sistema no tempo","Tempo (s)","Resposta", "Gama (graus)", "Theta (graus)","Alpha (graus)");
xtitle(T(1),T(2),T(3))
legends([T(4),T(5),T(6)],[2,3,4],1);
xgrid()

subplot(223)
plot2d(t,Thetadot*180/%pi,5,rect=[0,-10,1000,10])
T=list("ThetaDot no Tempo","Tempo (s)","ThetaDot (graus/s)");
xtitle(T(1),T(2),T(3))
xgrid()

subplot(224)
plot2d(X,Z,6)//,rect=[0,9900,300000,12400])
xtitle("Posicao do avião no espaço","Deslocamento (m)","Altura(m)")
xgrid()



scf(1) // Variaveis do sistema e entradas
subplot(231)   
plot2d(t,Velocidade,1)//,rect=[0,214,2000,247])
T=list("Resposta do sistema no tempo","Tempo (s)","Velocidade (m/s)","Velocidade")
xtitle(T(1),T(2),T(3))
xgrid()

subplot(232)
plot2d(t,Gama*180/%pi,1,rect=[0,-2,2000,7])
plot2d(t,Theta*180/%pi,2)
plot2d(t,Alpha*180/%pi,3)
T=list("Resposta do sistema no tempo","Tempo (s)","Resposta", "Gama (graus)", "Theta (graus)","Alpha (graus)");
xtitle(T(1),T(2),T(3))
legends([T(4),T(5),T(6)],[1,2,3],1);
xgrid()

subplot(233)
plot2d(t,Thetadot*180/%pi,5)//,rect=[0,-10,2000,10])
T=list("ThetaDot no Tempo","Tempo (s)","ThetaDot (graus/s)");
xtitle(T(1),T(2),T(3))
xgrid()

subplot(234)
plot2d(t,EntradaT,6)
T=list("Entradas Thrust no tempo","Tempo(s)","Entrada Thrust(Newtons)");
xtitle(T(1),T(2),T(3))
xgrid()

subplot(235)
plot2d(t,EntradaDelta*180/%pi,7)
T=list("Entradas Delta no tempo","Tempo(s)","Entrada Delta (rad)");
xtitle(T(1),T(2),T(3))
xgrid()

subplot(236)
plot2d(t,EntradaPert*PertV,1,rect=[0,-4000,2000,9000])
plot2d(t,EntradaPert*PertGama,2)
plot2d(t,EntradaPert*PertThetad,3)
T=list("Entradas perturbacao","Tempo(s)","Entradas", "PertV(N)", "PertGama(N)","PertThtad(Nm)");
xtitle(T(1),T(2),T(3))
legends([T(4),T(5),T(6)],[1,2,3],1);
xgrid()

/*scf(2) // Posicoes do sistema linearizado no tempo
subplot(221)
plot2d(t,Z,1,rect=[0,9900,2000,10400])
xtitle("Altura do avião ao longo do tempo","Tempo (s)","Altura (m)")
xgrid()

subplot(222)
plot2d(t,X,1,rect=[0,0,2000,300000])
xtitle("Posicao longitudinal do avião ao longo do tempo","Tempo (s)","Deslocamento (m)")
xgrid()

subplot(223)
plot2d(X,Z,1,rect=[0,9900,300000,10400])
xtitle("Posicao do avião no espaço","Deslocamento (m)","Altura(m)")
xgrid()
*/

//Solução do sistema não linear
//x1 = V, x2 = Gamma, x3 = thetadot, x4 = theta, u1 = T, u2 = delta
function xdot = Sistema(t,x)
    //Vetor com entradas
    u = [EntradaT(1)+Tcz;EntradaDelta(1)+deltaequilibrio]
    //Parâmetros que serão usados nos cálculos a seguir 
    alpha = x(4)-x(2) //ângulo de ataque
    cla = cl0a + clalphaa*alpha //Coeficiente de sustentação da asa
    cle = cl0e + clalphae*alpha + (cldeltae)*u(2) //Coeficiente de sustentação da empenagem + elevador
    cda = cd0a + (cla^2)/(%pi*ARa*ea) // coeficiente de arrasto da asa
    cde = cd0e + (cle^2)/(%pi*ARe*ee) // coeficiente de arrasto da empenagem
    
    La = rho*Aa*x(1)^2*cla/2 // Força de sustentação da asa
    Le = rho*Ae*x(1)^2*cle/2 // Força de sustentação da empenagem + elevador
    Da = rho*Aa*x(1)^2*cda/2 // Força de arrasto da asa
    De = rho*Ae*x(1)^2*cde/2 // Força de arrasto da empenagem
    Ma = rho*Aa*x(1)^2*cmaca*ca/2 // Momento gerado pela asa em seu centro 
    Me = rho*Ae*x(1)^2*cmace*ce/2 // Momento gerado pela empenagem em seu centro
    
    //Equações diferenciais do sistema
    xdot(1) = (u(1)*cos(alpha)-Da-De-W*sin(x(2)))/m //Vdot
    xdot(2) = (u(1)*sin(alpha)+La+Le-W*cos(x(2)))/(m*x(1)) //GammaDot
    xdot(3) = (La*xa*cos(alpha)+Le*xe*cos(alpha)+Da*xa*sin(alpha)+De*xe*sin(alpha)+Ma+Me)/I //ThetaDotDot
    xdot(4) = x(3) //ThetaDot
endfunction

//Condições iniciais
y0 = [Vcz+V0;Gamma0;ThetaDot0;alphacz+Theta0]
//Simulação das equações diferenciais
NaoLinear = ode(y0,0,t,Sistema)

//Cálculo dos deslcamentos
Z_NaoLinear(1,1)=h //Deslocamento vertical
X_NaoLinear(1,1)=0 //Deslocamento horizontal
for i = 1:length(t)-1
    Z_NaoLinear(:,i+1) = Z_NaoLinear(:,i) + p*NaoLinear(1,i)*sin(NaoLinear(2,i))
    X_NaoLinear(:,i+1) = X_NaoLinear(:,i) + p*NaoLinear(1,i)*cos(NaoLinear(2,i))
end

/*
//Gráficos comparando com o sistema linearizado
scf(3) // linear x nao linear 
subplot(2,2,1)
plot2d(t,NaoLinear(1,:),1,rect = [0,200,500,400])
plot2d(t,Velocidade,2)
xtitle("Comparação do modelo não linear com o linearizado","Tempo (s)","Velocidade (m/s)")
legends(["Não Linear","Linearizado"],[1,2],1)
xgrid()
subplot(2,2,2)
plot2d(t,NaoLinear(2,:)*180/%pi,1,rect = [0,-2,500,2])
plot2d(t,Gama*180/%pi,2)
xtitle("Comparação do modelo não linear com o linearizado","Tempo (s)","Ângulo de Trajetória (Graus)")
legends(["Não Linear","Linearizado"],[1,2],1)
xgrid()
subplot(2,2,3)
plot2d(t,NaoLinear(3,:)*180/%pi,1,rect = [0,-1,500,1])
plot2d(t,Thetadot*180/%pi,2)
xtitle("Comparação do modelo não linear com o linearizado","Tempo (s)","Taxa de variação do ângulo de arfagem (graus/2")
legends(["Não Linear","Linearizado"],[1,2],1)
xgrid()
subplot(2,2,4)
plot2d(t,NaoLinear(4,:)*180/%pi,rect = [0,0,500,4])
plot2d(t,Theta*180/%pi,2)
xtitle("Comparação do modelo não linear com o linearizado","Tempo (s)","ângulo de arfagem (graus)")
legends(["Não Linear","Linearizado"],[1,2],1)
xgrid()
*/

//Obtenção das funções de transferência e diagramas de Bode
aviao_TF = ss2tf(aviao)
//Nomenclatura da FT: i_j, onde i é a variável do espaço de estados e j é a entrada
V_T = aviao_TF(1,1)
V_Delta = aviao_TF(1,2)
V_Pert = aviao_TF(1,3)
Gamma_T = aviao_TF(2,1)
Gamma_Delta = aviao_TF(2,2)
Gamma_Pert = aviao_TF(2,3)
ThetaD_T = aviao_TF(3,1)
ThetaD_Delta = aviao_TF(3,2)
ThetaD_Pert = aviao_TF(3,3)
Tehta_T = aviao_TF(4,1)
Tehta_Delta = aviao_TF(4,2)
Tehta_Pert = aviao_TF(4,3)
//Listas para automatizar os plots
TFs = [V_T, V_Delta, V_Pert, Gamma_T, Gamma_Delta, Gamma_Pert, tf(1,1), ThetaD_Delta, ThetaD_Pert, Tehta_T, Tehta_Delta, Tehta_Pert]
Names = ['V_T', 'V_Delta', 'V_Pert', 'Gamma_T', 'Gamma_Delta', 'Gamma_Pert', 'ThetaD_T', 'ThetaD_Delta', 'ThetaD_Pert', 'Tehta_T', 'Tehta_Delta', 'Tehta_Pert']
//Plotando os diagramas de Bode
scf(5)
for j = 1:6
    subplot(2,3,j)
    bode(TFs(1,j),10^-3,10)
    title(Names(1,j))
end
scf(6)
for j = 7:12
    subplot(2,3,j-6)
    bode(TFs(1,j),10^-3,10)
    title(Names(1,j))
end

//Outra forma de obter FTs
[Ds,NUM,chi] = ss2tf(aviao)
FT = ss2tf(aviao)
Den = chi(1)

polos = roots(chi(1))

//Criterio de Routh
CoeficientesDaEqCarac = coeff(chi(1))
R(1,1) = CoeficientesDaEqCarac(5)
R(2,1) = CoeficientesDaEqCarac(4)
R(1,2) = CoeficientesDaEqCarac(3)
R(2,2) = CoeficientesDaEqCarac(2)
R(1,3) = CoeficientesDaEqCarac(1) 
R(2,3) = 0
R(3,1) = (R(2,1)*R(1,2)-R(1,1)*R(2,2))/R(2,1)
R(3,2) = (R(2,1)*R(1,3)-R(1,1)*R(2,3))/R(2,1)
R(3,3) = 0
R(4,1) = (R(3,1)*R(2,2)-R(3,2)*R(2,1))/R(3,1)
R(4,2) = (R(3,1)*R(2,3)-R(2,1)*R(3,3))/R(3,1)
R(5,1) = (R(4,1)*R(3,2)-R(4,2)*R(3,1))/R(4,1)

//Cálculo da Matriz de Transição (phi) para Δt = passo, até n = 4
//GAMMA é a matriz Γ que compõe o termo forçante

phi = eye(A) + A*p + (p^2*A^2)/2 + (p^3*A^3)/6 + (p^4*A^4)/24
GAMMA = A^-1*(phi-eye(A))

//Solução numérica pela matriz de transição
Transicao(:,1) = x0
u = [EntradaT(1);EntradaDelta(1); 0]
for i = 2:length(t)
    Transicao(:,i) = phi*Transicao(:,i-1) + GAMMA*B*u
end

scf()//Comparacao do linear com o metodo usando matriz de transicao
subplot(221)
plot2d(t,Transicao(1,:)+Vcz,rect=[0,220,100,250])
plot2d(t,Velocidade,2)
title("Comparação entre as respostas de Velocidade")
ylabel("Velocidade (m/s)")
xlabel("Tempo (s)")
legends(["Sistema obtido da matriz de transição","Sistema linearizado simulado"],[1,2],3)
xgrid()
subplot(222)
plot2d(t,Transicao(2,:)*180/%pi,rect=[0,-2,100,6])
plot2d(t,Gama*180/%pi,2)
title("Comparação entre as respostas de Gamma")
ylabel("Ângulo de trajetória (graus)")
xlabel("Tempo (s)")
legends(["Sistema obtido da matriz de transição","Sistema linearizado simulado"],[1,2],1)
xgrid()
subplot(223)
plot2d(t,Transicao(3,:)*180/%pi,rect=[0,-7,100,7])
plot2d(t,Thetadot*180/%pi,2)
title("Comparação entre as respostas de Theta Dot")
ylabel("Taxa de variação do ângulo de arfagem (graus/s)")
xlabel("Tempo (s)")
legends(["Sistema obtido da matriz de transição","Sistema linearizado simulado"],[1,2],3)
xgrid()
subplot(224)
plot2d(t,(Transicao(4,:)+alphacz)*180/%pi,rect=[0,-2,100,7])
plot2d(t,Theta*180/%pi,2)
title("Comparação entre as respostas de Theta")
ylabel("Ângulo de arfagem (graus)")
xlabel("Tempo (s)")
legends(["Sistema obtido da matriz de transição","Sistema linearizado simulado"],[1,2],3)
xgrid()
//
//Conferindo estabilidade com Transformação linear
[T,diagevals]=spec(A)
MatrizTransformada = T^-1*A*T
[RMT,diagevalsMT]=spec(MatrizTransformada)
