clear all
close all 
clc

%%  definizione input  
a_1 = 1.8; % autovalore 1 (DEFINIRA' L ANDAMENTO DELLA COMPONENTE X1 )
a_2 = 2.2; % autovalore 2 (DEFINIRA' L ANDAMENTO DELLA COMPONENTE X2 )
A = [a_1 0 ; 0 a_2] ; % matrice di stato 
B = [1;1]; % matrice d'ingresso
v_max =15; % passi ni 
col = 3; % inserire qui l'indice 
% indice   angolo considerato
%   1  ->  0
%   2  ->  pi/6
%   3  ->  pi/4
%   4  ->  pi/3
%   5  ->  pi/2
%   6  ->  3*pi/4
%   7  ->  -pi

[row_A,cols_A] = size(A); % definisco numero di righe e colonne di A
angoli = ([0, pi/6, pi/4, pi/3, pi/2, 3*pi/4, -pi]); % angoli su la circonferenza unitaria che voglio considerare
[row_angoli , cols_angoli] = size(angoli);
x_des = zeros(row_A, cols_angoli);

% definisco un vettore di ingressi la cui 
% PRIMA COMPONENTE e il coseno dell angolo 
% SECONDA COMPONENTE il seno dell angolo 

for r = 1:cols_angoli
    x_des(:,r) = ([cos(angoli(r)); sin(angoli(r))]);
end
fprintf("vettore di stato che il sistema deve raggiungere è  %f ,%f\n",x_des(1,col),x_des(2,col))
x = [x_des(1,col) ; x_des(2,col)];
rit = Plot_Andamento_U_X(A,B,v_max,x);
if(rit==-1)
    return;
end

Plot_Andamento_J(A,B,v_max,x);
% porzione di codice che stampa la circonferenza e gli angoli 
figure(4);
title("Stati considerati");
hold on;
grid on;
for r = 1:cols_angoli
    plot(x_des(1,r), x_des(2,r), 'o'); % disegna i i punti su la circonferenza presa
end
rad = 0:pi/64:pi;
plot(cos(rad), sin(rad), 'k');
xlabel("Coordinata X");
ylabel("Coordinata Y");




