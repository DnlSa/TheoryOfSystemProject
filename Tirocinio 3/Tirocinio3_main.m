clc
clear all 
close all 

%% NOTE 
% il sistema con due carrelli e vincolato allo scorrimento orizzontale 
% quindi cio che varia l unico vettore di stato x


%% AUTOVALORI COMPLESSI CONIUGATI  (CASO 1)
%   m1 = 2;
%   m2 = 2;
%   k_el = 4;
%   c = 2;

%% AUTOVALORI IMMAGINARI PURI (+- 2i)
   m1 = 1;
   m2 = 1;
   k_el = 2;
   c = 0;


%% AUTOVALORI REALI 
% m1 = 1;
% m2 = 1;
% k_el =2;
% c = 3;
%% Paramentri di simulaizone del sistema 

ts = 0.1 ; % sampling time 
max_time_simulation = 10 ; % tempo massimo di simulaizone 


%% VETTORE DI STATO INIZIALE DEL SISTEMA 
L = 5; % lunghezza della molla
q1 = 0;
q2 = 15;

x1 =  q1; % posizione iniziale carrello 1 (SAREBBE LA q1)
x2 = q2 - L;  % posizione iniziale carrello 2 (SAREBBE LA q2)
v1 = 0 ; % velocità iniziale carrello 1
v2  = 0 ; % velocità iniziale carrello 2
xE = [0 L]'; % posizione di equilibrio 
x0 = [x1; x2; v1 ; v2 ]; % vettore di stato inizile
e0 = x0 ; % stima dell'errore iniziale 


%% MATRICI DEL SISTEMA DINAMICO
A = [ 0    0     1     0; 
      0    0     0     1;
    -k_el/m1 k_el/m1 -c/m1  c/m1;
    k_el/m2 -k_el/m2  c/m2 -c/m2];

autovalori_A = eig(A)

B = [0 0 1/m1 0 ]'; 
M = [0 0 0 1/m2]';
C = [1 0 0 0]; 
D = 0; 

%% graficare andamento del sistema 

max_simulation = 20; % tempo massimo di simulazione 
ts=0.1 ; % tempo di campione 
t = 0:ts:max_simulation; % tempo di simulazione


% ingresso sinusoidale al sistema  
Amplitude= 1;
freq = 1;
u = Amplitude*sin(freq*t);

% ingresso costante al sistema 
%count = ones(size(t));
%u = 2*count; % ingresso costante

StatiUscita(A,B,C,D,t,u,x0)


%% Creazione delle matrici incerte 


d_1_k = 0.1 ; % incertezza su costante elastica 
d_1_c = 0;  % incertezza su attrito viscoso dello smorzatore
d_1_massa_1 = -0.01; % supponiamo che la massa del carrello 1 venga sottostimate di 0.5gr
d_1_massa_2  = 0.008;% suppongo che la massa del carrello 2 venga sovrastimata di 0.8 kg


m1_d = m1+d_1_massa_1;
m2_d = m2+d_1_massa_2;
k_el_d = k_el+d_1_k;
c_d = c+d_1_c;

A_d = [ 0    0     1     0; 
      0    0     0     1;
    -k_el_d/m1_d k_el_d/m1_d -c_d/m1_d  c_d/m1_d;
    k_el_d/m2_d -k_el_d/m2_d  c_d/m2_d -c_d/m2_d];

autovalori_A = eig(A)

B_d = [0 0 1/m1_d 0 ]'; 
M_d = [0 0 0 1/m2_d]';
C_d = [1 0 0 0]; 
D_d = 0; 


%% definizione di vettore di autovalori che dovrà avere (A-VC) e (A-BF)
 %Aut_val_distance = -0.5 ; % osservatore LENTO 
 %Aut_val_distance = -10 ; % osservatore VELOCE
 %Aut_val_distance = +1 ; % osservatore INSTABILE
 %aut_val_des = [-2,-3,-4,-1]; % se si desidera si possono impostare autovalori
 
 %all_real = 1 ; % 0 se vogliamo anche autovalori complessi 
                % 1 se autovalori devono essere tutti reali 

 %offset_autval = 1 ; % 0 se si desiderano autovalori impostati dall utente 
                     % 1 se si desiderano autovalori di offset impostati da Aut_val_distance 
 %aut_val_des = CalcoloAutovaloriDesiderati(A,Aut_val_distance,aut_val_des,offset_autval,all_real);
 
 %direttamente a mano 

 % impongo una distanza tra autovalori della matrice di stato A e quelli 
 % che deidero questi determineranno la dinamica del sistema
 % SI TENGA CONTO CHE : 
 % 1) se Re(lambda)<0 -> Sistema Asinototicamente Stabile
 % 2) se Re(lambda)<= 0 -> Sistema Stabile semplicemente 
 % 3) se Re(lambda)> 0 -> Sistema Instabile 
 

%% Progetto OSSERVATORE LENTO 
fprintf("\n\nProgetto osservatore LENTO\n\n");
Aut_val_distance = -0.5 ; % osservatore LENTO 
aut_val_des = [-1,-1.7,-2,-3]; % se si desidera si possono impostare autovalori 
all_real = 0 ;  % 0 se vogliamo anche autovalori complessi 
                % 1 se autovalori devono essere tutti reali 

offset_autval = 1 ; % 0 se si desiderano autovalori impostati dall utente 
                     % 1 se si desiderano autovalori di offset impostati da Aut_val_distance 

aut_val_des_slow = CalcoloAutovaloriDesiderati(A,Aut_val_distance,aut_val_des,offset_autval,all_real);

[F_slow,V_slow] = ProgettoOsservatore(A,B,C,aut_val_des_slow);

AVC_slow = (A-V_slow*C);
figure(1)
title("Dinamica d'errore con OSSERVATORE LENTO")
AndamentoErroreDiStima(AVC_slow,e0,ts,max_time_simulation)

%%  Progetto OSSERVATORE VELOCE
fprintf("\n\nProgetto osservatore VELOCE\n\n");

Aut_val_distance = -5 ; % osservatore VELOCE
aut_val_des = [-2,-3,-6,-3]; % se si desidera si possono impostare autovalori 

all_real = 0 ;  % 0 se vogliamo anche autovalori complessi 
                % 1 se autovalori devono essere tutti reali 

offset_autval = 1 ; % 0 se si desiderano autovalori impostati dall utente 
                     % 1 se si desiderano autovalori di offset impostati da Aut_val_distance 

aut_val_des_fast = CalcoloAutovaloriDesiderati(A,Aut_val_distance,aut_val_des,offset_autval,all_real);

[F_fast,V_fast] = ProgettoOsservatore(A,B,C,aut_val_des_fast);
AVC_fast = (A-V_fast*C);
figure(2)
title("Dinamica d'errore con OSSERVATORE VELOCE")
AndamentoErroreDiStima(AVC_fast,e0,ts,max_time_simulation)




%% PROGETTIAMO ADESSO UN OSSERVATORE OTTIMO

fprintf("\n\nProgetto osservatore OTTIMO\n\n");
% parametri di progettazione
rho = 0.1; % indice di incertezza piu e piccolo e piu mi fido 
gamma = 0.1 ; % incertezza
[K,S,CLP] = ProgettoOsservatoreOTTIMO(A,B,C,rho,gamma);  % Calcolata con i duali 

%[F_ott,V_ott] = ProgettoOsservatore(A,B,C,CLP'); 

V_ott = K';
F_ott = K;
AVC_ott= (A-V_ott*C);
aut_val_ottimi = eig(AVC_ott)
CLP
figure(3)
title("Dinamica d'errore con OSSERVATORE OTTIMO")
AndamentoErroreDiStima(AVC_ott,e0,ts,max_time_simulation)

%K =  è la matrice dei guadagni del controller LQR.
%S =  è la soluzione dell'equazione di Riccati associate al controller LQR.
%CLP =  contiene gli autovalori del sistema controllato.

%% grafici autovalori 
figure(4)
title("Grafico Autovalori su piano Complesso")
grafici(A,aut_val_des_slow,aut_val_ottimi,aut_val_des_fast)







