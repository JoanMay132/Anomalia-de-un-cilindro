%*******************************************************************%
%    DESCRIPCIÓN: INVERSIÓN DE LA ANOMALÍA DE UN CILINDRO           %
%                   ENTERRADO MEDIANTE MINIMOS CUADRADOS            %
%                       MÉTODO DE INVERSIÓN LOCAL                   %
%                                                                   %
%               ELABORADO POR: MAY DÍAZ JOAN CARLOS                 % 
%                FECHA DE MODIFICACIÓN: DICIEMBRE 2022              %
%*******************************************************************%

%% LIMPIEZA DE VARIABLES, GRAFICOS Y CONSOLA DE TRABAJO
tic
close all;
clear all;
clc


%% ANOMALIA DE GRAVEDAD PARA UN  CILINDRO (RESPUESTA TEORICA)
x= [-50:1:50]';     % Distancia-Posicion
z= 8; r= 5; drho= 0.9; % Parametros a buscar
k= 1e-20;      % Cte. de amortiguamiento

cte= 41.9e-9; % Constantes auxiliar
dg= (cte*drho*r^2)./(z.*(1+(x./z).^2)); % Anomalia teorica de un cilindro



%% DATOS SINTETICOS DE ANOMALIA DE GRAVEDAD
% Gráfica de anomalias

subplot(2,1,1) 
plot(x,dg,'b-'); title('Anomalia gravimétrica de un cilindro');
xlabel('Distancia-posición [m]'); ylabel('\Delta g [mGal]'); %Etiquetas
hold on;
grid on;

%% Construcción de la matriz jacobiana-kernel

z0=4; r0=3; drho0=0.5; % Modelo inicial

for i=1:20
    
dg1=(cte*drho0*r0^2)./(z0.*(1+(x./z0).^2));  % Anomalía para nuevo modelo (modelo de mejor ajuste)
    %%Jacobiano
dgdrho=(cte*(r0^2))./(z0*(1+(x./z0).^2)); %Derivada parcial dg/drho
dgdr=(2*cte*drho0*r0)./(z0*(1+(x./z0).^2)); %Derivada parcial dg/dr
dgdz=((-cte*drho0*(r0.^2)).*((z0.^2)-(x.^2)).*(((z0.^2)+(x.^2)).^-2)); %Derivada parcial dg/dz
       
          
G=[dgdrho dgdr dgdz]; %Construccion del kernell


dif=(dg-dg1); %Diferencia entre la respuesta teórica y el modelo
l2(i)=dif'*dif; %Cálculo del error con la norma L2
   
    
delta=((G'*G)+k*eye(length(G'*G)))^-1*G'*dif; %Vector de perturbaciones
    

drho0=drho0+delta(1); %Ajuste del nuevo parametro drho
r0=r0+delta(2);    %Ajuste del nuevo parametro r0
z0=z0+delta(3);    %Ajuste del nuevo parametro z0
 
plot(x,dg1,'c-'); pause(0.1); %Gráfica de la anomalía ajustada

legend('Respuesta teórica','Datos ajustados');
    
    if (l2(i)<=1e-25)   %tolerancia para error
        disp('Termino');
        break
    end
 
end
parametros = [z,z0; r,r0; drho,drho0]

%Gráfica de convergencia
subplot(2,1,2)
plot([1:i],l2,'*g-'); title('Gráfica del error con la norma L2'); xlabel(...
    'N° iteración'); ylabel('Norma L2'); legend('Convergencia'); grid on
toc
