function x = rispostaForzata(A,B,u,t_1)

if(t_1<0)
    fprintf("hai inserito un tempo negativo per favore inserisci un tempo>0")
    return 
end

% calcolo della risposta (APPLICO QUANTO SPIEGAGO NEL PROBLEMA PIU GENERALE DLE CONTROLLO )
% il calcolo avverrà completamente nel dominio di LAPLACE così da rendere
% tutto piu semplice 
% PASSO 1 -> calcolamo la matrice esponenziale e^At
syms s t tau
[row_A, ~]  = size(A);  % Calcolo dimensione di A 
Imat= eye(row_A); % creo una matrice identità della stessa dimensione 

% porzione di codice identica a quanto visto nel calcolo della gramiana 
rit = det((s*Imat-A)); 
if(rit == 0)% verifica dell effettiva invertibilità
    fprintf("la matrice non è invertibile quindi non possiamo andare avanti");
    return ; 
end

subs(u,tau,t); 
A_inv = inv((s*Imat-A)); % tocca anti-trasformarla con laplace (sI-A)^-1
u_s = laplace(u); % trasformo con laplace il vettore di controllo 

% calcoliamo la risposta forzata adesso in quanto la libera sarà nulla in quanto 
% il nostro sistema parte da una condizione inziale nulla come la traccia
% specifica 

% sappiamo gia da fodamenti di controlli che il calcolo della risposta 
% di un sistema e data da X(s)= ((sI-A)^-1)*x0 + ((sI-A)^-1)*B*u(s)

% ((sI-A)^-1)*x0 -> risposta libera del sistema nulla in quanto il nostro
% sistema parte nel suo punto di equilibrio x0=0

% ((sI-A)^-1)*B*u(s) -> risposta forzata che definisce il vettore di stato
% del nostro sistema 
x_s_l = 0 ; % risposta libera 0 per completezza didattica e preferibile includerla 
x_s_f= (A_inv*B)*u_s; % risposta forzata del nostro sistema 
X_s = x_s_l+x_s_f;
X_t = ilaplace(X_s) ; % antitrasformo per trovare la risposta nel dominio del tempo 
X = subs(X_t ,t,t_1); % valuto la risposta al tempo passatomi in input
%fprintf("la risposta forzata vale : ");
x = round(X,6);
end
