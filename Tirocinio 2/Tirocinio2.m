clear all 
close all 
clc
%% RIchieste 

% scegliere i seguenti autovalori 
% 1) due autovalori reali e distinti, entrambi negativi, 
%                gamma1 = -2;  
%                gamma2 = -3;
% autovalori (-1 , -2)  
% gramiana non complessa OK
% TEMPI DI SIMULAZIONE ACCETTABILI

% 2) due autovalori reali e distinti, uno negativo e uno non negativo
%       gamma1 = 3;  
%       gamma2 = -1; 
% autovalori (1,30 , -2.3O)   
% gramiana non complessa OK


% 3) due autovalori complessi coniugati 
%    gamma1 = -1; % 
%    gamma2 = -1; %
% autovalori (-0.5+0.8660i , -0.5-0.8660i )  
% gramiana non complessa OK
% TEMPI DI SIMULAZIONE POCO PIU DI 1 MINUTO

%% Dichiarazione delle matrici 
gamma1 = -2; % 
gamma2 = -3; %
syms t tau

A = [0 1 ; gamma1 gamma2] % matrice di stato 
B = [0;1] % matrice di ingressi 

opt_u = 0 ; % se settato a 1 abilita il calcolo per tutto l insieme dei tempi 
opt_x = 0 ;  % se settato a 1 abilita il calcolo per tutto l insieme dei tempi 
Maxima_COMPARE = 1; % Abilita il confronto con la matrice gramiana calcolata in MAXIMA
% ingresso definito dall utente  
u0 = sin(tau);
%u0 = 10*tau+5

T = 0 ; % se viene specificato un tempo T diverso da 0 verrà valutata la matrice gramiana per quel tempo specifico 
T_2 = 7; % tempo scelto dall utente per calcolare il confronto con t_bar
t_bar =7; % tempo in cui devo trovare il vettore di stato x_bar (stato finale )
max_time =14; % valore massimo dell'insieme dei tempi 
ts = 2; % distanza di campionamento fra un campione ed un altro (distanza per ogni Ti) 
ts_1 = 0.1; % tempo di campionamento per la simulazione richiesta la punto 4


%% Calcolo Autovalori dell amatrice A
[row_A, cols_A] = size(A);
if (row_A ~= cols_A)
        fprintf("la matrice non è quadrata\n");
        return   
end
disp("la Matrice A ha i seguenti autovalori: ");
aut_val = eig(A)

%% calcolo (PUNTO 1 )
G_t = CalcoloGramiana(A,B,T) % calcolo della gramiana 

%% MATRICI GRAMIANE CALCOLATE IN MAXIMA 
if(Maxima_COMPARE == 1) % gramiana prese da maxima Decommenteare il caso specifico 
    % caso con autovalori negativi (gamma : -2 , -3)
%G_t_1= [1/12-(exp(-4*t)*(6*exp(2*t)-8*exp(t)+3))/12, (exp(-4*t)*(exp(2*t)-2*exp(t)+1))/2;
 %       (exp(-4*t)*(exp(2*t)-2*exp(t)+1))/2, 1/6-(exp(-4*t)*(3*exp(2*t)-8*exp(t)+6))/6];
 % Nuova gramiana 
    G_t_1 = [1/12-(exp(-4*t)*(6*exp(2*t)-8*exp(t)+3))/12, (exp(-4*t)*(exp(2*t)-2*exp(t)+1))/2;
            (exp(-4*t)*(exp(2*t)-2*exp(t)+1))/2, 1/6-(exp(-4*t)*(3*exp(2*t)-8*exp(t)+6))/6];

 % caso con autovalori negativo e positivo (gamma: 3 , -1)
%   G_t_1 = [((((13^(3/2)+13)*exp(2*sqrt(13)*t)+312*exp(sqrt(13)*t)-13^(3/2)+13)*exp(-sqrt(13)*t-t))/12-169/6)/169, ((13*exp(2*sqrt(13)*t)-26*exp(sqrt(13)*t)+13)*exp(-sqrt(13)*t-t))/338;
%           ((13*exp(2*sqrt(13)*t)-26*exp(sqrt(13)*t)+13)*exp(-sqrt(13)*t-t))/338, (((13^(3/2)-13)*exp(2*sqrt(13)*t)-312*exp(sqrt(13)*t)-13^(3/2)-13)*exp(-sqrt(13)*t-t)+338)/676];
 % Nuova gramiana  
%G_t_1 = [((((13^(3/2)+13)*exp(2*sqrt(13)*t)+312*exp(sqrt(13)*t)-13^(3/2)+13)*exp(-sqrt(13)*t-t))/12-169/6)/169, ((13*exp(2*sqrt(13)*t)-26*exp(sqrt(13)*t)+13)*exp(-sqrt(13)*t-t))/338;
 %       ((13*exp(2*sqrt(13)*t)-26*exp(sqrt(13)*t)+13)*exp(-sqrt(13)*t-t))/338, (((13^(3/2)-13)*exp(2*sqrt(13)*t)-312*exp(sqrt(13)*t)-13^(3/2)-13)*exp(-sqrt(13)*t-t)+338)/676];


% caso con autovalori negativo e positivo (gamma: -1, -1)
% Nuova gramiana 
% G_t_1 = [(4*(3/8-(exp(-t)*(sqrt(3)*sin(sqrt(3)*t)-cos(sqrt(3)*t)+4))/8))/3, -(exp(-t)*(sqrt(3)*cos(sqrt(3)*t)-sqrt(3)))/3^(3/2);
 %       -(exp(-t)*(sqrt(3)*cos(sqrt(3)*t)-sqrt(3)))/3^(3/2), ((exp(-t)*(3^(3/2)*sin(sqrt(3)*t)+3*cos(sqrt(3)*t)-12))/2+9/2)/9];

 % vecchia gramaiana con unità immaginarie 
 %G_t_1 = [((((3^(3/2)*i+3)*exp(2*sqrt(3)*i*t)-24*exp(sqrt(3)*i*t)-3^(3/2)*i+3)*exp(-sqrt(3)*i*t-t))/4+9/2)/9, -((3*exp(2*sqrt(3)*i*t)-6*exp(sqrt(3)*i*t)+3)*exp(-sqrt(3)*i*t-t))/18;
    % -((3*exp(2*sqrt(3)*i*t)-6*exp(sqrt(3)*i*t)+3)*exp(-sqrt(3)*i*t-t))/18, (18-((3^(3/2)*i-3)*exp(2*sqrt(3)*i*t)+24*exp(sqrt(3)*i*t)-3^(3/2)*i-3)*exp(-sqrt(3)*i*t-t))/36];

end

%% Calcolo (PUNTO 2)
% il nostro sistema parte da 0 quindi non ha un evoluzione libera 
disp("lo stato finale del sistema è: ")
x_bar = rispostaForzata(A,B,u0,t_bar) % trovare il vettore di stato dato un ingresso e un tempo 

%% Calcolo (PUNTO 3)
%NOTA : il calcolo del indice di costo verrà eseguito per ultimo in quanto 
% richiede un tempo maggiore 

time= timevector(ts,max_time,t_bar); % Calcolo insieme dei tempi 

% Calcolo dei vettori di controllo ottimo passando la matrice Gramiana 
% di seguito calcoliamo  [u(tau)= B'*e^(t-tau)*beta] con [beta = G_t*x_bar]

u_star_T=TrovaControlloTC(A,B,G_t,time,x_bar); % CALCOLO vettori di controllo(nella variabile tau)
u_T_2=TrovaControlloTC(A,B,G_t,T_2,x_bar);

if(Maxima_COMPARE == 1)
     u_T_3=TrovaControlloTC(A,B,G_t_1,T_2,x_bar);
else 
    u_T_3 = 0; 
end



%% Graficazione dei vettori di controllo 
AndamentoU(u_star_T,time,t_bar,u0,ts_1 , opt_u);

%% Simulazione dell dell andamento del vettore di stato

AndamentoX(A,B,u_star_T,time,t_bar,u0,ts_1,opt_x,T_2,u_T_2,u_T_3);

%% calcolo dell indice di costo e graficazione 
disp("calcolo indice di costo") % utile per capire dove il codice e arrivato in run 
J_val = Calcolo_e_AndamentoJ(u_star_T,time,t_bar,u0,u_T_2,T_2)


