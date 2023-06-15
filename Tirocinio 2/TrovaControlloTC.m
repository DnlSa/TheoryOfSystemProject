function u_star_T = TrovaControlloTC(A,B,G_t,time,x_bar)
   % la funzione restiuira un vettore di controlli 
   if(time(1,1) == 0)
       disp("Inserire un tempo diverso da 0 ")
       u_star_T = 0;
       return;
   end

   len_time = size(time);
   [row_A, cols_A]  = size(A); % deve avere dimensione nxn 
   I = eye(row_A); % matrice identità

   %sfruttiamo il teorema della raggiungibilità e sappiamo che anche a
   %tempo continuo Xr = Im(G(t)) = Im(P) e che il sistema e raggiungibile 
   % se e solo se rank(P) = n la cui conseguenza è che
   % rank(G(t)) = n per ogni t != 0 per semplicità computazionale
   %possimo calcolare P. 

   % posso adoperare 
   P=[]; 
   for k=0:1:cols_A-1 % ciclo for per calcolare le matrici 
     mat_partial=(A^k)*B; % calcolo delle matrici parziali 
     P=horzcat(P,mat_partial); % concatenazione delle matrici 
   end
   P_xbar = horzcat(P,x_bar); % affianco P|x_bar
   rank_P_xbar = rank(P_xbar); % calcolo il rango della matrice P
   rank_P= rank(P);
   rank_A = rank(A); 
   if(rank_A~=rank_P) % ragiungibile generale 
       fprintf("il sistema non è completamente raggiungibile\n")
   end

    % Adesso ci calcoleremo tutti i stati [beta= G_t*x_bar]
    G_inv = inv(G_t); % non dovrebbe darmi problemi in quanto la matrice Gt e sempre invertibile 
    beta = G_inv*x_bar; % x_bar e il nostro stato finale desiderato 
    syms t tau s 
    A_inv_trasp = inv((s*I-A')); % (sI-A')^-1
    eAt_trasp= ilaplace(A_inv_trasp); % e^(A't)
    eAtau_trasp = subs(eAt_trasp,t,-tau);
    eA_t_tau = eAtau_trasp*eAt_trasp; % e^(A'(t-tau))
    u_star = (B')*(eA_t_tau)*beta; % calcolo controllo
    for k = 1:len_time(2) % calcolo dei controlli u_star_T dicendenza di tau
        u_star_T(k,1) = subs(u_star,t,time(1,k)); 
    end 
end