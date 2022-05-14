function [Sol] = dkinetics(theta)
%clear all; 

t1=theta(1);
t2=theta(2);
t3=theta(3);
t4=theta(4);
t5=theta(5);
t6=theta(6);

d1 = 47;
a1 = 30;
d2 = 45;
a2 = 120; 
a4 = 20;
d5 = 130;
d6 = 29;

%Determinar a matriz de parametros DH, tendo em conta as frames definidas
DH_matrix = zeros(8,4); 
DH_matrix(1,:) = [0,0,d1,t1];
DH_matrix(2,:) = [0,-a1,d2,0];
DH_matrix(3,:) = [sym(pi)/2,0,0,t2];
DH_matrix(4,:) = [0,a2,0,(sym(pi)/2)+t3];
DH_matrix(5,:) = [-(sym(pi)/2),a4,0,t4];
DH_matrix(6,:) = [0,0,d5,0];
DH_matrix(7,:) = [sym(pi)/2,0,0,t5];
DH_matrix(8,:) = [-sym(pi)/2,0,d6,t6];

%Elaboracao da matriz de transformacao 0->6
T_matrix=zeros(4,4,8);

for i = 1:8 
    m=zeros(4,4);
    m(1,:) = [cos(DH_matrix(i,4)) -sin(DH_matrix(i,4)) 0 DH_matrix(i,2)];
    m(2,:) = [sin(DH_matrix(i,4))*cos(DH_matrix(i,1)) cos(DH_matrix(i,4))*cos(DH_matrix(i,1)) -sin(DH_matrix(i,1)) -sin(DH_matrix(i,1))*DH_matrix(i,3)];
    m(3,:) = [sin(DH_matrix(i,4))*sin(DH_matrix(i,1)) cos(DH_matrix(i,4))*sin(DH_matrix(i,1)) cos(DH_matrix(i,1)) cos(DH_matrix(i,1))*DH_matrix(i,3)];
    m(4,:) = [0 0 0 1];
    T_matrix(:,:,i)=m;
   
   
end

T_final = T_matrix(:,:,1)*T_matrix(:,:,2)*T_matrix(:,:,3)*T_matrix(:,:,4)*T_matrix(:,:,5)*T_matrix(:,:,6)*T_matrix(:,:,7)*T_matrix(:,:,8);



% arrendondar p/ zero entradas mto pequenas?
%Matriz de rotacao, obtida atraves da matriz de transformacao
Rota = zeros(3,3);
Rota(1,:) = T_final(1,1:3);
Rota(2,:) = T_final(2,1:3);
Rota(3,:) = T_final(3,1:3);

%Definicao da orientacao, seguindo o sentido ZYX
c_beta = sqrt(Rota(1,1)^2 + Rota(2,1)^2); 
if ~(c_beta == 0)
    beta = atan2(-Rota(3,1), c_beta);
    alpha = atan2(Rota(2,1), Rota(1,1));
    gamma = atan2(Rota(3,2), Rota(3,3));
else
    alpha = 0;
    if -Rota(3,1) == 1
         beta = sym(pi)/2;
         gamma = atan2(Rota(1,2), Rota(2,2));
    else
        beta = -sym(pi)/2;
        gamma = -atan2(Rota(1,2), Rota(2,2));
    end
end

%Posicao
x = T_final(1,4);
y = T_final(2,4);
z = T_final(3,4);

Sol = [x, y, z, alpha, beta, gamma];

end




