 
function u = TrovaControllo(A,B,v_max,x)
% la funzione size mi restiuisce 2 valori 
% 1 -> il numero di righe 
% 2 -> il numero di Colonne
[~ , cols_A]  = size(A); % dimensione della matrice A
[row_B , cols_B] = size(B); % dimensione della matrice B
v = v_max; % inutile ma ho uniformato la firma degli argomenti  

if(cols_B>2) % controllo su la matrice di ingresso
    disp("la matrice di ingresso presenta PIU DI 2 INGRESSI ")
    u =-1;
    return
end

% contorllo su l'esatta dimensione delle matrici 
if(row_B  ~= cols_A )
  disp("Passa una matrice ..." + ...
        "A con dimensione n*n..." + ...
        "B con dimensiuone p*n di cui p = numero di ingressi ")
   u =-1;
  return 
end

% Creazione della matrice di raggiungibilità
% P= [B|AB|A^2B|....|A^k-1B];
P=[]; % diachiaro un matrice vuota altrimenti mi da errore
for k=0:1:cols_A-1 % ciclo for per calcolare le matrici 
    mat_partial=(A^k)*B; % calcolo delle matrici parziali 
    P=horzcat(P,mat_partial); % concatenazione delle matrici 
end
q = rank(P); % calcolo il rango della matrice P
% ne verifico la raggiungibilità
% Caso 1 -> il sistema non e completamente raggiungibile quindi procederemo
% a verificare singolarmente se ogni vettore di stato e raggiungibile 
% Caso 2 -> il sistema è completamente raggiungibile sapendo che rank(P)=n 
if( q ~= cols_A) 
    disp("il sistema non e completamente raggiungibile " + ...
        "verificare la raggiungibilià dei vettori di input");
end


% DA QUI IN POI ADOTTEREMO L'ALGORITMO SPIEGATO PER TROVARE IL VETTORE DI
% CONTROLLO CON IL COSTO MINORE (spiegato nella lezione del 20/03/2023)


% PUNTO 1 calcolo di una matrice P_bar
P_bar = [] ; 
for i = 0:v-1 % da 0 a ni-1
   P_bar = horzcat(P_bar,(A^i)*B);
end

% verifica dell esistenza della soluzione per il controllo a minima energia
s=rank(P_bar);
P_bar_x = [P_bar,x];  %concateno vettore di stato alla matrice P(creo la matrice completa  )
r = rank(P_bar_x) ;

% il rango della matrice deve rimanere invariato
if(r == s ) % rank(P_bar)==rank(P_bar_x) 
   % fprintf("il sistema ammette %d soluzioni quindi lo stato finale " + ...
    %    "è raggiungibile  \n",a);
else 
    disp("il sistema non ammette soluzioni ");
    u =-1;
    return ;
end

% PUNTO 2-BIS calcolo una matrice G quadrata 
% ricordiamo che Im(P)= Im(P*P')

% punto A) 
G = P_bar*(P_bar'); 

% punto B) calcolare : [ beta = G_inv*x ] 
% beta = vettore di controllo 
% x = vettore di stato 
% mi occorre calcolare la matrice inversa a tal proposito 
% adotto una matrice pseudoinversa in quanto non è detto che la matrice
% G possa essere invertibile 

G_inv  = inv(G);  % calcolo pseudoinversa
beta = G_inv*x; % trovo il vettore beta

% punto C) calcoliamoci [ w = P_bar'*beta] 
w = (P_bar')*beta;

% PUNTO 3 ->   partizioniamo w in sottovettori di controllo 
% u(.) con p componenti 

% calcolo dimensione di w , mi serviranno solamente le righe di w
[row_w , ~] = size(w); 
j = row_w;
u = zeros(cols_B,1); % definisco un prevettore di controllo con p componenti 

% ricordiamo che i sotto vettori verranno scritti al contrario 
% partendo da u(v-1) fino all ultimo che sarà u(0)
% cio va fatto perche inizialmente ho posto P = [B|AB|....|A^(n-1)*B]
if(cols_B ==2) % ogni sotto vettore avra p colonne e una sola riga
    for i = 0:1:v-1
            u = horzcat(u,w(j-1:j));
            j = j-2;
    end
else  % altrimenti cols_B ha 1 sola colonna i sottovettori u avranno p=1 componenti
    for i = 0:v-1
        u=horzcat(u,w(j));
        j = j-1;
    end
end
%eliminazione della prima colonna cosi si ottiene una matrice che ha
%per colonne i segnali di controllo u da applicare in istanti diversi
u(:,1) = [];

% adesso verificheremo l'esattezza del risultato 
% si applica la solita verifica usata nell esercizio 5.2.2
% x(k+1)= A*x(k)+B*u(k); 

x_k = zeros(cols_A,1); % x e preso in input 
for  b= 1:1:v % calcolo vettori di stato per ogni passo 
     x_k = (A*x_k)+(B*u(:,b)); 
end

c1 = (x_k(1) - x(1));
c2 = (x_k(2) - x(2));
c1 = round(c1,4); % arrotonda fino alla 4 cifra decimale
c2 = round(c2,4);

% prendo solamente i valori interi 
a = zeros(cols_A, 1); 
if(c1 ~=0 || c2 ~= 0)
      fprintf("i vettori di stato calcolati sono\n%d - %d\n%d - %d\n" + ...
         "non è stato possibile raggiungere\nlo stato desiderato in %d passi\n",x_k(1),x_k(2),x(1),x(2),v_max);
      u=-1;
     return
 end
  %fine verifica vettori
  %fprintf("il controllo u(.) calcolato permette di raggiungere lo stato x\n");

end