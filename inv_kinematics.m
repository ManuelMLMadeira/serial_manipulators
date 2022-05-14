function [theta_f] = invkinetics(V)
%clear all; 

if ~(length(V)==6)
    disp ('6 arguments are required.')

else
    x = V(1);
    y = V(2);
    z = V(3);
    a = V(4);
    b = V(5);
    g = V(6);


    d1 = 47;
    a1 = 30;
    d2 = 45;
    a2 = 120; 
    a4 = 20;
    d5 = 130;
    d6 = 29;



    %matriz de rotacao obtida atraves de alpha, beta, gamma, seguindo a
    %convencao usada na cinetica direta - ZYX
    Rota = [cos(a)*cos(b), cos(a)*sin(b)*sin(g)-sin(a)*cos(g), cos(a)*sin(b)*cos(g)+sin(a)*sin(g); 
        sin(a)*cos(b), sin(a)*sin(b)*sin(g)+cos(a)*cos(g), sin(a)*sin(b)*cos(g)-cos(a)*sin(g);
        -sin(b), cos(b)*sin(g), cos(b)*cos(g)];

    %Definicao das coordenadas da origem da frame 5, em relacao ao referencial
    %da base
    p0_6 = [x, y, z];
    p5_6 = [0, 0, d6];
    p0_5 = (-Rota * p5_6' + p0_6')';
    
    [x5, y5, z5] = deal(p0_5(1), p0_5(2), p0_5(3));


    %Obtencao de theta1, theta2 e theta3 pelo metodo geometrico



    %Existem 4 solucoes para os valores de theta2 e theta3

    t123_m = zeros(4,3);
    
    if x5 == 0 && y5==0 && x == 0 && y == 0
        disp('SINGULARITY: There are infinite solutions for the value of theta1 and theta6. By default, theta1 will assume the values 0 and pi. The value of theta6 must take into account other angles, and will be defined accordingly.')
    end
    
    % Right-armed Configuration
    % (definicao de Theta 1 - primeiras 2 solucoes)
    t1 = atan2(y5 ,x5);
    
    valid = true;
    while valid

        %definicao de direcoes e angulos comuns
        xy2_5 = sqrt((x5^2)+(y5^2)) + a1; %atentar no sinal de a1
        z2_5 = z5-(d1+d2);
        d2_5 = sqrt(xy2_5^2+z2_5^2);
        t_aux1 = atan2(z2_5,xy2_5);

        d3_5 = sqrt(d5^2 + a4^2);
        phi = atan2 (d5,a4);
        
        %seguindo a regra de carnot
        c1 = ((d3_5^2)+(a2^2)-(d2_5^2))/(2*d3_5*a2);
        if abs (c1)>1
            string1 = ['for theta1 =', num2str(t1),', the coordinates inserted are out of reach.'];
            disp(string1)
            t123_m(1,:) = [NaN, NaN, NaN];
            t123_m(2,:) = [NaN, NaN, NaN];
            break
        end
        t_aux2 = acos (c1);
        
        c2 =(-(d3_5^2)+(a2^2)+(d2_5^2))/(2*a2*d2_5);
        if abs (c2)>1
            string2 = ['for theta1 =', num2str(t1),', the coordinates inserted are out of reach.'];
            disp(string2)
            t123_m(1,:) = [NaN, NaN, NaN];
            t123_m(2,:) = [NaN, NaN, NaN];
            break
        end
        t_aux3 = acos(c2);
        
        

        %Primeira solucao
        t2 = t_aux1 - t_aux3; 
        t3 = sym(pi)/2 - phi - t_aux2;

        t123_m(1,:) = [t1, t2, t3];

        %Segunda solucao
        t2 = t_aux1 + t_aux3;
        t3 = sym(pi)/2 - phi + t_aux2;

        t123_m(2,:) = [t1, t2, t3];
        
        valid = false;
        
    end 
    
    % Left-armed Configuration
    % (definicao de Theta 1 - ultimas 2 solucoes)
    t1 = double(sym(pi)) + t1;
    
    valid = true;
    while valid
        
        %definicao de direcoes e angulos comuns
        xy2_5 = -(sqrt(x5^2+y5^2) - a1); %atentar no sinal de a1
        z2_5 = z5-(d1+d2);
        d2_5 = sqrt((xy2_5)^2+(z2_5)^2);
        t_aux1 = atan2(z2_5,xy2_5);

        d3_5 = sqrt((d5^2) + (a4^2));
        phi = atan2 (d5,a4);
        
        %seguindo a regra de carnot
        c1 = (d3_5^2+a2^2-d2_5^2)/(2*d3_5*a2);
        if abs (c1)>1
            string1 = ['for theta1 =', num2str(t1),', the coordinates inserted are out of reach.'];
            disp(string1)
            t123_m(3,:) = [NaN, NaN, NaN];
            t123_m(4,:) = [NaN, NaN, NaN];
            break
        end
        t_aux2 = acos (c1);
        
        c2 =(-d3_5^2+a2^2+d2_5^2)/(2*a2*d2_5);
        if abs (c2)>1
            string2 = ['for theta1 =', num2str(t1),', the coordinates inserted are out of reach.'];
            disp(string2)
            t123_m(3,:) = [NaN, NaN, NaN];
            t123_m(4,:) = [NaN, NaN, NaN];
            break
        end
        t_aux3 = acos(c2);


        %Terceira solucao
        t2 = t_aux1 - t_aux3;
        t3 = sym(pi)/2 - t_aux2 - phi;

        t123_m(3,:) = [t1, t2, t3];

        %Quarta solucao
        t2 = t_aux1 + t_aux3;
        t3 = sym(pi)/2 + t_aux2 - phi;

        t123_m(4,:) = [t1, t2, t3];

        valid = false;
    end


    % Metodo Algebrico para obtencao de theta4, theta5, theta6

    theta_f = zeros(8,6);

    for i=1:4
        if ~isnan(t123_m(i,1))

            %Definir t2, t2, t3
            t1 = t123_m(i,1);
            t2 = t123_m(i,2);
            t3 = t123_m(i,3);

            % Matriz de transformacao homogenea 0->3
            DH_matrix = zeros(4,4); 
            DH_matrix(1,:) = [0,0,d1,t1];
            DH_matrix(2,:) = [0,-a1,d2,0];
            DH_matrix(3,:) = [sym(pi)/2,0,0,t2];
            DH_matrix(4,:) = [0,a2,0,sym(pi)/2+t3]; 

            T_matrix=zeros(4,4,4);
            for j = 1:4 
                m=zeros(4,4);
                m(1,:) = [cos(DH_matrix(j,4)) -sin(DH_matrix(j,4)) 0 DH_matrix(j,2)];
                m(2,:) = [sin(DH_matrix(j,4))*cos(DH_matrix(j,1)) cos(DH_matrix(j,4))*cos(DH_matrix(j,1)) -sin(DH_matrix(j,1)) -sin(DH_matrix(j,1))*DH_matrix(j,3)];
                m(3,:) = [sin(DH_matrix(j,4))*sin(DH_matrix(j,1)) cos(DH_matrix(j,4))*sin(DH_matrix(j,1)) cos(DH_matrix(j,1)) cos(DH_matrix(j,1))*DH_matrix(j,3)];
                m(4,:) = [0 0 0 1];
                T_matrix(:,:,j)=m;
            end

            T0_3 = T_matrix(:,:,1)*T_matrix(:,:,2)*T_matrix(:,:,3)*T_matrix(:,:,4);

            % Inversa de T0_3 pelas suas propriedades de mat de rotacao
            R0_3 = T0_3(1:3,1:3);
            P0_3 = T0_3(1:3,4);
            inv_T0_3 = zeros(4,4);
            inv_T0_3(1:3,1:3) = R0_3';
            inv_T0_3(1:3,4) = -R0_3'*P0_3;
            inv_T0_3(4,4) = 1;

            % Matriz de transformacao homogenea 0->6
            T0_6 = zeros(4,4);
            T0_6(1:3,1:3) = Rota;
            T0_6(1:3,4) = p0_6;
            T0_6(4,4) = 1;

            % Produto da inversa de transf. hom. 0->3 com transf. hom. 0->6
            T_prod = inv_T0_3*T0_6;

            % Determinar t4, t5 e t6
            
            zeros_bool = false;
            t5_1 = acos(T_prod(2,3));
            t5_2 = -t5_1;
            if T_prod(2,3) == 1 || T_prod(2,3) == -1
                zeros_bool = true;
                disp('SINGULARITY: There are infinite solutions for the value of theta4 and theta6. These angles will then assume symmetrical values.')
            end
            
            % Se  T_prod(2,2)==0 && T_prod(2,1)==0 && ~zeros_bool =>
            % => sin(t6)=0 e cos(t6)=0.
            % Trata-se de um erro na computacao. O valor para t6 vai ser,
            % por default, 0.
            t6_1 = atan2(-T_prod(2,2), T_prod(2,1));
            t6_2 = atan2(T_prod(2,2),-T_prod(2,1));
            
            % Se  T_prod(3,3)==0 && T_prod(1,3)==0 && ~zeros_bool =>
            % => sin(t4)=0 e cos(t4)=0.
            % Trata-se de um erro na computacao. O valor para t4 vai ser,
            % por default, 0.
            t4_1 = atan2(T_prod(3,3), -T_prod(1,3));
            t4_2 = atan2(-T_prod(3,3), T_prod(1,3));
            
            
            % Formar solucoes
            theta_f(2*i-1,:) = [t1, t2, t3, t4_1, t5_1, t6_1];
            theta_f(2*i,:) = [t1, t2, t3, t4_2, t5_2, t6_2];



        else
            theta_f(2*i-1,:)=[NaN,NaN,NaN,NaN,NaN,NaN];
            theta_f(2*i,:)=[NaN,NaN,NaN,NaN,NaN,NaN];

        end
    end
end

end






