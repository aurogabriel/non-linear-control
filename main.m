%% Variables
clc
clear
close all

%Valores Tabelados:
Rs=4.3;         %Resistência do stator
Ld=27e-3;       %Indutância no eixo d
Lq=67e-3;       %Indutância no eixo q
np=4;           %Numero de Polos
Lambda_m=0.279; %PM Flux Linkage 
J=1.79e-4;      %Momento de Inercia


%Valores não tabelados:
Kf=0.001;        %Coeficiente de atrito 
T_load=1;       %Torque carga

%Fatores de controle
q=100;
p=(Ld^2/Lq^2)*q;
r=0;

%% Diagrama de phases Malha Aberta

% Malha aberta, u=0
syms x1 x2 x3
a = 0.2;
[x1, x2, x3] =meshgrid(-2e3:4e2:2e3, -2e3:4e2:2e3, -1700:340:1700);
% x3=1700;
% Calcula as derivadas
dx1 = (-Rs/Ld)*x1 + np*(Lq/Ld)*x2.*x3;
dx2 = (-Rs/Lq)*x2 - np*(Ld/Lq)*x1.*x3 - np*(Lambda_m/Ld)*x3;
dx3 = (-Kf/J)*x3-(T_load/J);

% % Plota o campo vetorial
figure
quiver3(x1, x2, x3, dx1, dx2, dx3);
% streamslice(x1,x2,dx1,dx2)
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Diagrama de Phases em Malha Aberta')
hold on
x0 = [2000 2000 300];
[t, x] = ode45(@MalhaAberta,[0 5e-2], x0);
plot3(x(:,1),x(:,2),x(:,3));
legend('Diagrama de Fases','Plot 3d do Caminho')

figure
plot(t,x(:,1),t,x(:,2),t,x(:,3))
legend('x1','x2','x3')
title('Variaveis de Estado no Tempo')
grid on

%% Diagrama de phases Malha Fechada

syms x1 x2 x3
% Parâmetros do grid
a = 0.8;
[x1, x2] = meshgrid(-2e3:4e2:2e3, -2e3:4e2:2e3);

% Cálculo das variáveis de controle u(x)
u1 = (-p/Ld)*x1;
u2 = (-q/Lq)*x2;
u3 = (r/J)*x3; 
x3=1700;
% Equações diferenciais
dx1 = (-Rs/Ld)*x1 + np*(Lq/Ld)*x2.*x3 + (u1/Ld);
dx2 = (-Rs/Lq)*x2 - np*(Ld/Lq)*x1.*x3 - np*(Lambda_m/Ld)*x3 + (u2/Lq);
% dx3 = (-Kf/J)*x3 - (T_load/J) + (u3/J);
dx3=0;
% % Plota o campo vetorial usando quiver3

figure
streamslice(x1,x2,dx1,dx2)
xlabel('x1');
ylabel('x2');
title('Diagrama de Fases em Malha Fechada de X1, X2', 'para velocidade Nominal')


[x1, x2, x3] = meshgrid(-2e3:4e2:2e3, -2e3:4e2:2e3, -1700:340:1700);
% Cálculo das variáveis de controle u(x)
u1 = (-p/Ld)*x1;
u2 = (-q/Lq)*x2;
u3 = (r/J)*x3; 

% Equações diferenciais
dx1 = (-Rs/Ld)*x1 + np*(Lq/Ld)*x2.*x3 + (u1/Ld);
dx2 = (-Rs/Lq)*x2 - np*(Ld/Lq)*x1.*x3 - np*(Lambda_m/Ld)*x3 + (u2/Lq);
dx3 = (-Kf/J)*x3 - (T_load/J) + (u3/J);

figure
quiver3(x1, x2, x3, dx1, dx2, dx3); 
title('Diagrama de Fases em Malha Fechada de X1, X2, X3')
hold on
x0 = [200 200 200];
[t, x] = ode45(@MalhaFechada,[0 5e-9], x0);
plot3(x(:,1),x(:,2),x(:,3));
legend('Diagrama de Fases','Plot 3d do Caminho')

%% Verificar dissipatividade estrita global
pvar x1 x2 x3 u1 u2 u3;
dpvar Q R S rho;

vars = [x1;x2;x3;u1;u2;u3];
prog = sosprogram(vars);

x = [x1;x2;x3];
U = [u1;u2;u3];

% Definição de g como um vetor de polinômios
g1 = -3.0074e10;
g2 = -7.4627e10;
g3 = 0;
g = [g1; g2; g3];


f = [(-Rs/Ld)*x1+np*(Lq/Ld)*x2*x3;(-Rs/Lq)*x2-np*(Ld/Lq)*x1*x3-np*(Lambda_m/Ld)*x3;(-Kf/J)*x3-(T_load/J)];

% Monomios   
h = x;
prog = sosprogram(vars); %inicializa o programa sos

[prog,V] = sospolyvar(prog,monomials(x1,2:4),'wscoeff');
[prog,T] = sospolyvar(prog,monomials(x,1:2),'wscoeff');

Vx1 = diff(V,x1);
Vx2 = diff(V,x2);
Vx3 = diff(V,x3);

GradV = [Vx1 Vx2 Vx3];
prog = sosdecvar(prog,(rho));
[prog, Q] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), 
length(h)], 'symmetric');
[prog, S] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), 
length(U)]);
[prog, R] = sospolymatrixvar(prog, monomials(vars, 0), [length(U), 
length(U)], 'symmetric');

expr = -((GradV*(f+g.*U)+ T) -h'*Q*h - 2*h'*S*U - U'*R*U );

prog = sosineq(prog,expr);
solver_opt.solver = 'sedumi';

prog = sossolve(prog,solver_opt);
q = sosgetsol(prog,Q);
r = double(sosgetsol(prog,R));
s = sosgetsol(prog,S);
delta = s*inv(r)*s'-q;

disp('Delta=')
eig(double(delta))
disp('V=')
sosgetsol(prog,V)
disp('T=')
sosgetsol(prog,T)

%% Function defines

function dx= MalhaFechada(t,x)

%Valores Tabelados:
Rs=4.3;         %Resistência do stator
Ld=27e-3;       %Indutância no eixo d
Lq=67e-3;       %Indutância no eixo q
np=4;           %Numero de Polos
Lambda_m=0.279; %PM Flux Linkage 
J=1.79e-4;      %Momento de Inercia


%Valores não tabelados:
Kf=0.02;        %Coeficiente de atrito 
T_load=1;       %Torque carga

%Fatores de controle
q=100;
p=Ld^2/Lq^2;
r=1e1;

u1 = (-p/Ld)*x(1);
u2 = (-q/Lq)*x(2);
u3 = (r/J)*x(3);
dx=[(-Rs/Ld)*x(1) + np*(Lq/Ld)*x(2)*x(3) + (u1/Ld); 
    (-Rs/Lq)*x(2) - np*(Ld/Lq)*x(1)*x(3) - np*(Lambda_m/Ld)*x(3) + (u2/Lq);
    (-Kf/J)*x(3) - (T_load/J) + (u3/J)];
end 


function dx= MalhaAberta(t,x)
%Valores Tabelados:
Rs=4.3;         %Resistência do stator
Ld=27e-3;       %Indutância no eixo d
Lq=67e-3;       %Indutância no eixo q
np=4;           %Numero de Polos
Lambda_m=0.279; %PM Flux Linkage 
J=1.79e-4;      %Momento de Inercia


%Valores não tabelados:
Kf=0.02;        %Coeficiente de atrito 
T_load=1;       %Torque carga

%Fatores de controle
q=100;
p=(Ld^2/Lq^2)*q;
r=0;

dx=[(-Rs/Ld)*x(1) + np*(Lq/Ld)*x(2)*x(3) ; 
    (-Rs/Lq)*x(2) - np*(Ld/Lq)*x(1)*x(3) - np*(Lambda_m/Ld)*x(3);
    (-Kf/J)*x(3) - (T_load/J)  ];
end 
