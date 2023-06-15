 clc
 clear all 
 close all 

%% NOTE 
% il sistema con due carrelli e vincolato allo scorrimento orizzontale 
% quindi cio che varia l unico vettore di stato x

    
%1)  se i coefficenti di smorzamento sono maggiori delle masse questi avranno 
% dei coefficenti REALI  indice che i 2 carrelli avranno un moto relativo
% NON oscillante dovuto alla molla e allo smorzatore che li collega 
% (CASO 3)

%2) se le masse sono minori  della costante elastica e presentano un
% coefficente di smorzamento c diverso 0 si genereranno dei moti
% oscillatori smorzati . Tale smorzamento e prettamente dovuto allo
% smorzatore che introduce dei poli che sono discostanti dall asse
% immaginario

% 3) se invece ho una C = 0 e moto relativo fra i 2 carrelli sarà
% prettamente oscillatorio presentando dei autovalori con sola parte
% immaginaria 


%% Definizioni delle variabili simboliche 

syms t s a 
fig=0 ; % tiene traccia delle figure
%% AUTOVALORI COMPLESSI CONIUGATI  (CASO 1)
%   m1 = 2;
%   m2 = 2;
%   k_el = 4;
%   c = 2;

%% AUTOVALORI IMMAGINARI PURI (+- 2i)
%   m1 = 2;
%   m2 = 2;
%   k_el = 1;
%   c = 0;

%%

m1 = 1;
m2 = 1;
k_el =2;
c = 3;


%% Paramentri di simulaizone del sistema 

ts = 0.1 ; % sampling time 
max_time_simulation = 10 ; % tempo massimo di simulaizone 

%% segnale in ingresso CASO 1  

M1  = 2 ; 
M0  = 3 ; 
r_t_1 = M1*t^2+M0; 
r_s_1 =laplace(r_t_1);

[pole_r1,N_1] = poles(r_s_1);

dim_pole_r1 = size(pole_r1); 
for k = 1 : dim_pole_r1
    if(pole_r1(k,:)==0)
        types_r1 =  cast(N_1(k,:),'uint8');
    end
end

%% segnale in ingresso CASO 2

OM = 1.5 ;  % omega di pulsazione (conosco)
M  = 1 ; % ampiezza sinusoide 
f  = 1 ; % sfasato di 45°  (scelta arbitraria )

r_t_2 = M*sin(OM*t+f);
r_s_2 = expand(laplace(r_t_2)); % trasformata di laplace del segnale 

[pole_r2,N_2] = poles(r_s_2);
dim_pole_r2 = size(pole_r2);
for k  = 1 : dim_pole_r2
    if(pole_r2(k,:)==0)
        types_r2 =  cast(N_2(k,:),'uint8');
    end
end



%% VETTORE DI STATO INIZIALE DEL SISTEMA 
L = 5; % lunghezza della molla
q1 = 0;
q2 = 10;

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

disp("Autovalori della matrice A: ")
disp(eig(A))
B = [0 0 1/m1 0 ]'; 
C = [1 0 0 0]; 
D = 0; 


%% matrici con incertezza 
d_1_k = 0.8 ; % incertezza su costante elastica 
d_1_c = 0.2;  % incertezza su attrito viscoso dello smorzatore
d_1_massa_1 =  0.5; % supponiamo che la massa del carrello 1 venga sottostimate di 0.5gr
d_1_massa_2  = 0.2;% suppongo che la massa del carrello 2 venga sovrastimata di 0.8 kg


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

[Nsys_unc,Dsys_unc] = ss2tf(A_d,B_d,C_d,D_d) ;

%% definizione di vettore di autovalori che dovrà avere (A-VC) e (A-BF)
 Aut_val_distance = -0.5 ; % osservatore LENTO 
 %Aut_val_distance = -10 ; % osservatore VELOCE
 %Aut_val_distance = +1 ; % osservatore INSTABILE 
 
 aut_val_des = [-2,-2,-1,-1]; % se si desidera si possono impostare autovalori
 %direttamente a mano 
 all_real = 1 ; % 0 se vogliamo anche autovalori complessi 
                % 1 se autovalori devono essere tutti reali 

 offset_autval = 0 ; % 0 se si desiderano autovalori impostati dall utente 
                     % 1 se si desiderano autovalori di offset impostati da Aut_val_distance 

 % impongo una distanza tra autovalori della matrice di stato A e quelli 
 % che deidero questi determineranno la dinamica del sistema
 % SI TENGA CONTO CHE : 
 % 1) se Re(lambda)<0 -> Sistema Asinototicamente Stabile
 % 2) se Re(lambda)<= 0 -> Sistema Stabile semplicemente 
 % 3) se Re(lambda)> 0 -> Sistema Instabile 
 
%% Calcolo del Del numero di poli in 0 presenti nela funzione d'impianto

[Nsys,Dsys] = ss2tf(A,B,C,D) ; % numeratore e denominatore della fdt del sistema 
types = 0 ; 
pole_sys = roots(Dsys); % devo determinare il tipo del processo 
num_pole = size(pole_sys);
for k = 1:num_pole(1)
    ret_im = round(imag(pole_sys(k,1)));
    ret_re = round(real(pole_sys(k,1)));
    if(ret_im == 0 && ret_re ==0)
        types= types+1;
    end
end
fprintf("la funzione di trasferimento del sistema e di tipo %d\n",types);

%% PROGETTAZIONE DEL COMPENSATORE NEL CASO TRAIETTORIA r(t) = M1*t+M0   

%PROGETTAZIONE DI Cm(s)
if(types <= types_r1)
    disp("non occorre progettare una Cm(s) stabilizzante ")
else
    disp("occorre implementare una Cm(s) o cambiare il segnale di riferimento")
    return;
end
%In questo caso l'inseguimento di traiettoria è del tipo r(t) = M1*t+M0
%quindi r(s) ha polo 0 con molteplicità 2. 
%La funzione di trasferimento P(s) del sistema ha già quel polo con la
%stessa molteplicità, quindi Cm(s) = 1.
    

% utilizzando il vecchio script del tirocinio 3 e immediato
% calcolare le due matrici in quanto precedentemente ho implementato
% tutti i possibili casi 

[F,V] = ProgettoOsservatore(A,B,C,aut_val_des);
Ac_1 = A-(V*C)+(B*F)-(V*D*F);
Bc_1 = V;
Cc_1 = -F;
Dc_1 = 0;
[Ncom_1,Dcom_1] = ss2tf(Ac_1,Bc_1,Cc_1,Dc_1);
%fig=StabilityTest (Nsys,Dsys,Ncom_1,Dcom_1,r_s_1,fig);

%% PROGETTAZIONE DEL COMPENSATORE CASO riferimento r(t) = M*sin(OM*t+f) 
% (PUNTO 4)
% mi tocca progettare una Cm(s) = 1/fi(s) 
% trasformata del sengnale di riferimento è OM^2/(s^2+OM^2)
% da questo segnale so che per la progettazione del blocco stabilizzante 
% mi occorre un controllore lamento di tipo 2 se ne e sprovvista la
% funzione di impianto (G(s)). la forma di fi(s) = (s^2+OM^2);
% inoltre ci occorre evitare che avvengano cancellazioni illecite

fprintf("\n\nProgetto compensatore con segnale(t) = M*sin(OM*t+f) \n\n");


% Progetto della Cm(s)
% definisco Cm(s) = 1/fi(s) con fi(s) = (s^2 + OM^2).

Aut_val=eig(A);
Aut_val_dim = size(Aut_val) ; % numero di autovalori di A

for k = 1: k
    Aut_val_imag = imag(Aut_val(k,1)); % prendo parte immaginaria
    if(OM ==Aut_val_imag || -OM ==Aut_val_imag  ) % discrimino omega e il suo complesso coniugato 
        disp("cancellazione illecita considerare una OM diversa")
        return 
    end
end
% scrivo le matrici definite al PASSO 2
% fi(s) = (s^2+OM^2) -> s^2 +  0    + OM^2*s^0
% fi(s) = (s^2+OM^2)->  b^m  + b^(1)  + b^0





% m_seg = 2 % polinomio di secondo grado
% n = m_seg * p = 2 * 1 = 2 
% h=2
% REALIZZAZIONE IN FORMA CANONICA DI RAGGIUNGIBILITA'

dim_B = size(B); % calcolo dimensione degli ingressi 
p = dim_B(2);
Ip = eye(p); % matrice identità di dimensione p

% dichiaro coefficenti matrice Am
beta_1 = 0  ;
beta_0 = OM^2;

% calcolo coefficenti di MARKOV
W_0 = D;

% W_h = C*A^(h-1)*B  con h= 1,2,3,.....
W_1 = C*(A^0)*B; 
W_2 = C*(A^1)*B;

% matrici del blocco stabilizzante 
Am = [  0      -beta_0*Ip    ; 
        Ip     -beta_1*Ip   ];


 Bm = [Ip ;
       0 ];

 Cm = [ W_1 W_2];
 
 Dm = W_0;



 %ora si trova una rappresentazione dello spazio di stato per la serie
 %di Cm(s) e P(s)
    dim_Am = size(Am);
    dim_A = size(A);
    mat = zeros(dim_Am(1),dim_A(2));  % completamento 

    % forma A_bar = [Am , 0 ; B*Cm , A]
    % B ha dimensione 4                                                                                         
    A_bar = [Am , mat; B*Cm , A];
    B_bar = [Bm;
             B*Dm];
    C_bar = [D*Cm,C];
    D_bar = D*Dm;

 % PROGETTAZIONE DI Cs(s)

 aut_val_des_A_bar = CalcoloAutovaloriDesiderati(A_bar,Aut_val_distance,aut_val_des,offset_autval,all_real);
    
 % devo rimodificare sia l osservatore che la matrice di retro 
 [F1,V1] = ProgettoOsservatore(A_bar,B_bar,C_bar,aut_val_des_A_bar); 
 As = A_bar - (V1*C_bar) + (B_bar*F1) - (V1*D_bar*F1);
 Bs = V1;
 Cs = -F1;
 Ds = 0;

   %ora si scrive la rappresentazione nello spazio di stato per C(s) che è
   %la serie di Cs(s) e Cm(s)
   dim_As = size(As);
   mat = zeros(dim_As(1),dim_Am(2)); 
   Ac_2 = [As mat;Bm*Cs Am];
   Bc_2 = [Bs;Bm*Ds];
   Cc_2 = [Dm*Cs Cm];
   Dc_2 = Dm*Ds;
   [Ncom_2,Dcom_2] = ss2tf(Ac_2,Bc_2,Cc_2,Dc_2); % Da matrici dello spazio di stato a funzione di trasferimento 
%    fig = StabilityTest (Nsys,Dsys,Ncom_2,Dcom_2,r_s_2,fig);
% 
%   

