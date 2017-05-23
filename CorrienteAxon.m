function CorrienteAxon(X, T, h, k, lambda, tao, I, s)
%CorrienteAxon Algoritmo para resolver la EDP de un axón en reposo al que se
% le aplica una corriente durante un lapso de tiempo siguiendo la
%Teoria del cable
%   EL PVI a resolver es
%   [  lambda^2(V_xx)=tao(V_t)+V         x in [0,B], t in T,
%   [  V(x, 0)  = V0 (-60mv),              x in [0,B]
%   [  V(0,t)= I si t<seg (corriente aplicada en mV durante una cantidad de segundos)
% 
%
%   Con condiciones de frontera Neuman homogeneas
%   Ejemplo: CorrienteAxon([0,1],[0,1],.1,.0001,2, 1, 10, .05)
%       Axón de largo 1, al que se le aplica una corriente de 10mv 5
%       milesimas de segundo

%Primero lo primero
close
tpausa = (k);

xm = X(1); xM = X(2);
x = xm:h:xM;
M =length(x);

if length(T) == 2
    t0 = T(1); tF = T(2);
else
    t0 = 0; tF = T;
end
pasos = ceil((tF - t0)/k);
Vviejo=zeros(1,M);
Vnuevo=zeros(1,M);

%Implementamos el método, tomando en cuenta que en x=0 estamos aplicando
%una corriente y en xM hay condiciones de Neuman nulas
for i=1:pasos
    for j=1:M
        % En x=0 vamos a aplicar una corriente I una cantidad seg de
        % segundos
        if(j==1)
            segundo=i*k
            if(segundo<s)
            Vnuevo(1)=I; %se aplica la corriente
            else Vnuevo(1)=Vviejo(2); %Si no, la frontera es neuman
            end
        else if(j==M)
            Vnuevo(M)=Vnuevo(M-1); %La frontera es neuman
            else
              %El esquema
            Vnuevo(j)=(lambda.^2/tao)*(k/h.^2)...
            *(Vviejo(j+1)-2*Vviejo(j)+Vviejo(j-1))... 
                +Vviejo(j)*(1-k/tao);
            end
        end
    end %Acaba el cálculo de Vnuevo

    
plot(x,Vnuevo)
axis([xm,xM,0,I])
pause(tpausa)
for j=1:M
Vviejo(j)=Vnuevo(j);
end
end
end