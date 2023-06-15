%% calcolo dell amatrice gramiana 
% Script di calcolo della matrice GRAMIANA DI RAGGIUNGIBILITA' 
% il calcolo di una matrice esponenziale dobbiamo avvalerci della
% trasformata di laplace quindi dobbiamo calcolare
% (laplace^-1){(s*I-A))}

 
%% CALCOLIAMO LA FUNZIONE INTEGRANDA PER IL CALCOLO DELL AMATRICE GRAMIANA 

function G_t = CalcoloGramiana(A,B,T)

syms s t tau  % variabili simboliche per il calcolo della matrice esponenziale 
[n,~] =size(A); 
Imat = eye(n);

rit = det((s*Imat-A)); 
if(rit == 0)% verifica dell effettiva invertibilità
    fprintf("la matrice non è invertibile quindi non possiamo andare avanti");
    G = -1; 
    return ; 
end

A_inv = inv((s*Imat-A)); % tocca anti-trasformarla con laplace (sI-A)^-1 
eAt= ilaplace(A_inv); % calcolo della matrice esponenziale 
eAtB = (eAt*B);  % primo pezzo della funzione integranda per il calcolo della gramiana 

A_inv_trasp = inv((s*Imat-A')); % tocca anti-trasformarla con laplace (sI-A')^-1 

eAt_trasp= ilaplace(A_inv_trasp); % calcolo della matrice esponenziale 
eAtB_trasp = ((B')*eAt_trasp);  % parte trasposta della  funzione integranda per il calcolo della gramiana 



%NOTA : ho invertito t con tau in quanto t mi viene restituita
%automaticamente dall antitrasformata di laplace . inoltre tau viene spesso
%adottato come notazione alternativa a t per definire dei tempi allora si e
%deciso di  invertire queste 2 notazioni .
% quindi adesso integreremo la funzione in t i cui estremi vanno da 0 a tau
fun =  (eAtB*eAtB_trasp); % funzione integrande in 

if(T==0)
    % calcolo diretto dell integrale 
    G_t = int(fun,t , 0 , tau );
    G_t = subs(G_t,tau,t); % sostituisco a t con tau cosi mi torna una matrice gramiana espressa nella variabile t 
else
    G_t = int(fun,t , 0 ,T);
end
% calcolo con trasformata di laplace 
%G_s= laplace(int(fun,t , 0 , tau ));
%G_t = ilaplace(G_s); % calcolo integrale in t e con estremo di integrazione tau

end

