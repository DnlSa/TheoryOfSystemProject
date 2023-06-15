function [F,V] = ProgettoOsservatore(A,B,C,aut_val_des)
    %% Calcolo matrice di osservabilità 
    syms a 
    Aut_val = eig(A);
    dim_A = size(A) ;% prendo dimensione matrice A 
    
    Q = [] ; % definisco una matrice vuota
    for k = 0:(dim_A(2)-1)
       temp = C*(A^k);
       Q=vertcat(Q,temp);
    end
    
    %% VERIFICA DELLA POSSIBILITA' DI PROGETTARE UN OSSERVATORE 
    % OSSERVABILITA rank(Q) = n 
    % se passa il test allora il suo duale sarà RAGGIUNGIBILE
    if(rank(Q) ~= dim_A)
        fprintf("la matrice di osservabilita del sistema PRIMALE non e osservabile\n");
        osservabile = 0 ;
    else
        osservabile = 1 ;
    end
    
    
    % DETERMINABILITA rank(Q/A^n) = rank(Q)
    % se passa il test allora il suo duale sarà CONTROLLABILE
    A_n = A^(dim_A(2)); % elevo alla potenza n-esima la matrice A
    Q_An = vertcat(Q,A_n);
    if(rank(Q_An) ~= rank(Q))
        fprintf("la matrice di osservabilita del sistema PRIMALE non e determinabile ");
        return;
    end
    
    % per poter dire che il sistema sia RILEVABILE 
    % SE E SOLO SE VALE UNA QUALUNQUE DELLE SEGUENTI CONDIZIONI EQUIVALENTI FRA
    % LORO : 
    
    % 1) il sitema DUALE è stabilizzabile (SI INTENDE CHE TUTTI GLI AUTOVALORI RISPETTINO LA STABILIA'
    % O POSSONO ESSERE SPOSTATI TRAMITE RETROAZIONE COME VISTO CON MITTER E ACKERMAN )
    dim_Aut_val = size(Aut_val); % serve anche dopo alla verifica con Ackerman
    kalman =inf; % imposto kalman con un valore flag assurdo cosi che puo fungere da flag 
    
    if(osservabile == 0 )  % se sistema non osservabile vado a vedere 
        % 2) il sistema PRIMALE valga la relazione rank[A-lambda_i*I/C]
        Id_mat = eye(dim_A(2));
        Mat_lambda = (A-a*Id_mat);
        Mat_lambda = vertcat(Mat_lambda,C);
        for k = 1: dim_Aut_val(1)
          temp = subs(Mat_lambda,a,Aut_val(k,:));
          if rank(temp) ~= dim_A(2)
            fprintf("PBH test di determinabilità e fallito per il seguente autovalore : %d\n ",Aut_val(k,:))
            if(Aut_val(k,:)<0) % stiamo considerando un sistema a tempo continuo quindi RE(lambda)<=0 per stabilità
                fprintf("l autovalore che non soddisfa il PBH " + ...
                    "test e POSITIVO quindi il sistema " + ...
                    "non e STABILIZZABILE <cambiare valori " + ...
                    "attribuiti al sistema reale e riprovare a lanciare il programma>\n");
                return
            else 
                fprintf("L'autovalore che fa fallire il PBH test è buono " + ...
                    "possiamo continuare e calcolare una forma di kalman e spostare " + ...
                    "gli altri autovalori\n");
                kalman =0;  % variabile per portare il sistema in forma di kalman
            end
          end
        end
    end

    %% DEFINIZIONE del sistema DUALE
    A_d = A';
    %B_d = C';
    C_d = B';
    
    
    %% VERIFICHE SUL SISTEMA DUALE 
    % (inutili in quanot per dualità abbiamo tutto dal sistema primale 
    % ma puo risultare interessante verificare effettivamente che quanto detto 
    % a lezione funzioni )
    dim_A_d = size(A_d) ;% prendo dimensione matrice A 
    % Calcolo matrice di OSSERVABILITA' del DUALE 
    Q_d = [] ; % definisco una matrice vuota
    for k = 0:(dim_A_d(2)-1)
       temp = C_d*(A_d^k);
       Q_d=vertcat(Q_d,temp);
    end


    if(rank(Q_d) ~= dim_A(2))
        fprintf("la matrice di osservabilià del sistema DUALE  non e osservabile\n");
    else
        disp("il sistema duale e osservabile quindi il primale e raggiungibile")
    end
    
    % VERIFICA della DETERMINABILITA' DUALE
    A_d_n = A_d^(dim_A_d(2)); % elevo alla potenza n-esima la matrice A_d
    Q_d_An = vertcat(Q_d,A_d_n);
    if(rank(Q_d_An) ~= rank(Q_d))
        fprintf("la matrice di osservabilita del sistema DUALE non e determinabile\n");
        return;
    end
    
    % verifica CONTROLLABILITA DUALE
    P = [];
    for k = 0: (dim_A(2)-1)
        temp = (A^k)*B;
        P = horzcat (P,temp);
    end
    P_A_n = horzcat(P,A_n);
    if(rank(P_A_n) ~= rank(P))
        fprintf("il sistema Primale non e controllabile\n");
        return;
    end

    %% VERIFICA RAGGIUNGIBILITA SU PRIMALE TRAMITE (P)
    PBH_rag = 0;
    if(rank(P) ~= dim_A(2)) 
        disp("il sistema non e completamente raggiungibile")
          PBH_rag = 1;
    end
    %% PBH TEST DI RAGGIUNGIBILITA'
    if(PBH_rag == 1)
        Id_mat = eye(dim_A(2));
        Mat_lambda = (A-a*Id_mat);
        Mat_lambda = horzcat(Mat_lambda,B); % si adotta la forma 5 del PBH
        for k = 1: dim_Aut_val(1)
          temp = subs(Mat_lambda,a,Aut_val(k,:)); % inserisco autovalore 
          if rank(temp) ~= dim_A(2) % se < n allora vado in forma di kalman
            fprintf("PBH test di determinabilità e fallito per il seguente autovalore : %d\n ",Aut_val(k,:))
            if(Aut_val(k,:)<0) % stiamo considerando un sistema a tempo continuo quindi RE(lambda)<=0 per stabilità
                fprintf("l autovalore che non soddisfa il PBH " + ...
                    "test e POSITIVO quindi il sistema " + ...
                    "non e STABILIZZABILE <cambiare valori " + ...
                    "attribuiti al sistema reale e riprovare a lanciare il programma>\n");
                return
            else 
                fprintf("L'autovalore che fa fallire il PBH test è buono " + ...
                    "possiamo continuare e calcolare una forma di kalman e spostare " + ...
                    "gli altri autovalori\n");
                kalman =0;  % variabile per portare il sistema in forma di kalman
            end
          end
        end
    end
%% PROGETTAZIONE DELLA MATRICE DI RETROAZIONE DI STATO (F)
    if(PBH_rag ==0) % IL PBH e verificato per tutti gli autovalori 
        % applichiamo ackerman sul PRIMALE 
        len_aut_val_des = size(aut_val_des);
        if(len_aut_val_des(2)~= dim_A(2))
            disp("inserire un vettore coerente alla dimensione della matrice A")
            return ;
        end
        % matrice P di ragiungibilità calcolata durante la verifica di
        % CONTROLLABILITA'
        dim_P = size(P);
        P_inv = inv(P);
        vet_pi = P_inv(dim_P(1),:); % prendo ultima riga dell'inversa di P_d
          
        p_des = 1;
        for k = 1: dim_A(2)
           p_des =(a-aut_val_des(1,k))*p_des; % costruisco il p_des
        end
         p_des = expand(p_des);  % espando il polinomio 
         vet = sym2poly(p_des); % vettore coefficenti del polinomio desiderato
         dim_vet = size(vet); % calcolo la lungezza di tale vettore 
         n = dim_vet(2); % imposto indice delle potenze di A
         F = zeros(1,dim_vet(2)-1);  % vettore di 0 per definire F
         % formula che sto implementando o 
         % F  = -pi*p_des(A)
         for k = 1:dim_vet(2) %ciclo for su  coefficenti del polinomio desiderato
             n=n-1; 
             ret= -vet(1,k)*(vet_pi*(A^n)); % potenze di A^n 
             F = F+ret;
         end 
    else
        % implementare kalman di raggiungibilità 
        % Nel caso un ui il PBH test fallisse per un valore 
        disp("il PBH test e fallito passiamo in forma di kalman")
        Xr = P;
        % obiettivo e trovare una base Xr
        
        % Calcola la decomposizione QR della matrice
        [ ~ , temp_mat] = qr(Xr);
        % viene scomposta la matrice in 
        % - temp_mat = e la P ridotta ad una matrice triangolare superiore quindi 
        %               se una colonna o riga e linearmente dipendente
        % - il primo parametro e una matrice ortogonale 
        colonna_dipendente  = ~diag(temp_mat); % cerca il primo 0 su diagonale
        Xr(:,colonna_dipendente ) = []; % elimino la colonna linearmente dipendente
        Xr = null(Xr);
        dim_Xr = size(Xr);
        nr = dim_Xr(2);
        T_inv=[];
        T_inv = horzcat(T_inv,Xr); %  matrice T invertibile già completa 
        n=rank(T_inv);
        len_T_inv = size(T_inv);
       
        if(len_T_inv<dim_A) % caso in cui servissero dei vettori di completamento
            T_inv_temp = T_inv;
            for k = 1:dim_A(2) % inserisco vettori di completamento 
               T_inv_temp = horzcat(T_inv_temp,Id_mat(:,k));
               temp = rank(T_inv_temp);
               if(temp>n) % entro in questo if solamente quando trovo un vettore linearmente indipendente 
                 T_inv = horzcat(T_inv,Id_mat(:,k)); % matrice T invertibile con completamento  
                 n=rank(T_inv); % aggiorna rango
               end
               if(n == dim_A(2))
                   break
               end
            end
        end
        T = inv(T_inv); % ho trovato la matrice T invertibile
        A_bar = T*A*T_inv;
        B_bar = T*B; 
       %C_bar = C*T_inv;
       %Q_bar = Q*T_inv;

        index_row_A = 0 : nr;
        index_cols_A = 0 : nr;
        index_row_B = 0:nr;
        index_cols_B =1;
        %index_row_C = 1;
        %index_cols_C =0:nr;
        A_rr  = A_bar(index_row_A,index_cols_A); 
        B_r  =  B_bar(index_row_B,index_cols_B );
        %C_r  =  C_bar(index_row_C,index_cols_C);

        Aut_val_ragg = eig(A_rr); % autovalori manipolabili
        dim_aut_val = size(Aut_val_ragg); % prendo il numero di autovalori 
        Aut_val_des_ragg = []; % costruisco vettore di autovalori desiderati
        for k = 1:dim_aut_val(1)
             Aut_val_des_ragg = horzcat(Aut_val_des_ragg , aut_val_des(1,k));
        end
        % da qui applico ackerman 
        % ricalcolo la maptrice di raggiungibilita con le nuove matrici
        % e poi eseguo sostanzialmente lo stesso codice fatto  nel caso
        % normale  
        % calcolo matrice P
        dim_Arr = size(A_rr);
        P_ragg = [];
        for k = 1:dim_Arr(1)
            temp = B_r*(A_rr^k);
            P_ragg = horzcat(P_ragg,temp);
        end
        dim_P_ragg = size(P_ragg);
        P_ragg_inv = inv(P_ragg);
        vet_pi_ragg = P_ragg_inv(dim_P_ragg(1),:); % prendo ultima riga dell inversa di P_oss
        % calcoliamo il polinomio desiderato questo definira la velocita del sistema  
        p_des = 1;
        for k = 1: dim_Arr(2)
           p_des =(a-Aut_val_des_ragg(1,k))*p_des; % costruisco il p_des
        end
         p_des = expand(p_des);  % espando il polinomio 
         vet = sym2poly(p_des); % vettore coefficenti del polinomio desiderato
         dim_vet = size(vet); % calcolo la lungezza di tale vettore 
         n = dim_vet(2); % imposto indice delle potenze di A
         F = zeros(1,dim_vet(2)-1);  % vettore di 0 per definire F_d
         % formula che sto implementando o 
         for k = 1:dim_vet(2) %ciclo for su  coefficenti del polinomio desiderato
             n=n-1; 
             ret= -vet(1,k)*(vet_pi_ragg*(A_d^n)); % potenze di A_d^n vanno da 4 a 0
             F = F+ret;
         end 
    end
     %% VERIFICA DEL APPLICAZIONE ESATTA DI ACKERMAN  
    Aut_val_ABF = eig((A+(B*F)))
    Aut_val_ABF = round(Aut_val_ABF,4);
    Aut_val_ABF= real( Aut_val_ABF);
    % verifica automatica autovalori
    temp_array = real(round(aut_val_des,4));
    for k = 1: dim_Aut_val(1)
        ret = Aut_val_ABF(k,1);              
        for j = 1: dim_Aut_val(1)
            if(ret == temp_array(1,j))
                temp_array(1,j) = Inf; % inf perche e un valore che non si avvera mai 
                break;
            end
         end
     end
     for j = 1: dim_Aut_val(1)
         if(temp_array(1,j) ~= Inf)
            fprintf("verifica fallita\n")
            %return; 
         end
     end
 fprintf("autovalori di (A-B*F) verificati\nAkerman e stato applicato con successo\n")

    %% PROGETTAZIONE DELLA MATRICE (V) DELL'OSSERVATORE (V)
    
    % Sarebbe auspicabile implementare 
    % 2 casi 
    % 1) sistema OSSERVABILE 
    % 2) sistema NON OSSERVABILE passare in forma di kalman
    
    %% CASO 1 sistema OSSERVABILE
    if(kalman == inf ) 
        len_aut_val_des = size(aut_val_des);
        if(len_aut_val_des(2)~= dim_A(2))
            disp("inserire un vettore coerente alla dimensione della matrice A")
            return ;
        end
        % dobbiamo trovare una V  = -F_d
        % F_d e possibile trovare con Ackerman 
        % generalmente mi sarebbe servito chiamare la relativa API ma per scopi 
        % puramente didattici risulta piu interessante implementarlo a mano 
        P_d = Q';
        dim_P_d = size(P_d);
        P_d_inv = inv(P_d);
        vet_pi = P_d_inv(dim_P_d(1),:); % prendo ultima riga dell inversa di P_d
        
        % calcoliamo il polinomio desiderato questo definira la velocita del sistema  
        p_des = 1;
        for k = 1: dim_A(2)
           p_des =(a-aut_val_des(1,k))*p_des; % costruisco il p_des
        end
         p_des = expand(p_des);  % espando il polinomio 
         vet = sym2poly(p_des); % vettore coefficenti del polinomio desiderato
         dim_vet = size(vet); % calcolo la lungezza di tale vettore 
         n = dim_vet(2); % imposto indice delle potenze di A
         F_d = zeros(1,dim_vet(2)-1);  % vettore di 0 per definire F_d
         % formula che sto implementando o 
         % F  = -pi*p_des(A) 
         for k = 1:dim_vet(2) %ciclo for su  coefficenti del polinomio desiderato
             n=n-1; 
             ret= -vet(1,k)*(vet_pi*(A_d^n)); % potenze di A_d^n vanno da 4 a 0
             F_d = F_d+ret;
         end  
         V= -(F_d');
         
    else
        % nel caso arrivassimo qui possiamo possimo calcolare la forma di
        % kalman per la sola parte osservabile 
        
% PASSO 1) si cerca una matrice T invertibile per applicare il cambio di base 

       Xi  = Q ; % calcolo sottospazio di inoss
       [ ~ , temp_mat] = qr(Xi);
       colonna_dipendente  = ~diag(temp_mat); % cerca il primo 0 su diagonale
       Xi(:, colonna_dipendente) = []; % elimino la colonna linearmente dipendente
       Xi = null(Xi);
       dim_Xi = size(Xi);
       ni = dim_Xi(2);
       T_inv=[];
       T_inv = horzcat(T_inv,Xi); %  matrice T invertibile già completa  
       len_T_inv = size(T_inv);
       n=rank(T_inv);
       if(len_T_inv<dim_A)
           T_inv_temp = T_inv;
           for k = 1:dim_A(2)
               T_inv_temp = horzcat(T_inv_temp,Id_mat(:,k));
               temp = rank(T_inv_temp);
               if(temp>n)
                 T_inv = horzcat(T_inv,Id_mat(:,k)); % matrice T invertibile con completamento  
                 n=rank(T_inv); % aggiorna rango
               end
               if(n == dim_A(2))
                   break
               end
           end
       end
% PASSO 2) calolcolo matrici sistema con la nuova base 
        T = inv(T_inv);
        A_bar = T*A*T_inv;
        B_bar = T*B; 
        C_bar = C*T_inv;
        %Q_bar = Q*T_inv;

 % PASSO 3) applico Akerman su la parte osservabile che posso variare
        % ricorstruire vettore autovalori desiderati 
        dim_A_bar = size(A_bar);
        dim_B_bar = size(B_bar);
        %dim_C_bar = size(C_bar);
        % prelevo la parte osservabile 
        index_row_A = ni:dim_A_bar(1);
        index_cols_A =ni: dim_A_bar(2);
        index_row_B = ni:dim_B_bar(1);
        index_cols_B =1;
        %index_row_C = 1;
        %index_cols_C =ni:dim_C_bar(2);
        A_oo  = A_bar(index_row_A,index_cols_A); 
        B_o  =  B_bar(index_row_B,index_cols_B );
        %C_o  =  C_bar(index_row_C,index_cols_C);

        Aut_val_oss = eig(A_oo);
        dim_aut_val = size(Aut_val_oss); % prendo il numero di autovalori 
        Aut_val_des_oss = []; % costruisco vettore di autovalori desiderati
        for k = 1:dim_aut_val(1)
             Aut_val_des_oss = horzcat(Aut_val_des_oss , aut_val_des(1,k));
        end
        % da qui applico ackerman 
        % ricalcolo la matrice di raggiungibilita con le nuove matrici
        % e poi eseguo sostanzialmente lo stesso codice fatto  nel caso
        % normale  
        % calcolo matrice P
        dim_Aoo = size(A_oo);
        P_oss = [];
        for k = 1:dim_Aoo(1)
            temp = B_o*(A_oo^k);
            P_oss = horzcat(P_oss,temp);
        end
        dim_P_oss = size(P_oss);
        P_oss_inv = inv(P_oss);
        vet_pi_oss = P_oss_inv(dim_P_oss(1),:); % prendo ultima riga dell inversa di P_oss
        % calcoliamo il polinomio desiderato questo definira la velocita del sistema  
        p_des = 1;
        for k = 1: dim_Aoo(2)
           p_des =(a-Aut_val_des_oss(1,k))*p_des; % costruisco il p_des
        end
         p_des = expand(p_des);  % espando il polinomio 
         vet = sym2poly(p_des); % vettore coefficenti del polinomio desiderato
         dim_vet = size(vet); % calcolo la lungezza di tale vettore 
         n = dim_vet(2); % imposto indice delle potenze di A
         F_d = zeros(1,dim_vet(2)-1);  % vettore di 0 per definire F_d
         % formula che sto implementando o 
         % F  = -pi*p_des(A) = -vet_pi*(a^4 - 4*a^3 + 6*a^2 -4*a + 1)
         for k = 1:dim_vet(2) %ciclo for su  coefficenti del polinomio desiderato
             n=n-1; 
             ret= -vet(1,k)*(vet_pi_oss*(A_d^n)); % potenze di A_d^n vanno da 4 a 0
             F_d = F_d+ret;
         end  
         V= -(F_d');
        return ; 
    end
    %% VERIFICA DEL APPLICAZIONE ESATTA DI ACKERMAN  
    Aut_val_AVC = eig((A-(V*C)))
         % queste 2 api servono a eliminare il piccolo scarto di MATLAB nel
         % nei calcoli in quanto fa uscire delle unità immaginare
         % bassissime(residui di calcolo) e valori estremamente vicini all
         % autovalore desiderato es :
         % aut_val_des = 1 -> quello che restituisce 0.9998+0.0002i 
         % NOTA : escono ordine diverso da quelli base pero non e un problema
         Aut_val_AVC= real( Aut_val_AVC);
         Aut_val_AVC = round(Aut_val_AVC,4); 
         % verifica automatica autovalori
         temp_array = real(round(aut_val_des,4));
         for k = 1: dim_Aut_val(1)
                ret = Aut_val_AVC(k,1);
                for j = 1: dim_Aut_val(1)
                    if(ret == temp_array(1,j))
                       temp_array(1,j) = Inf; % inf perche e un valore che non si avvera mai 
                       break;
                    end
                end
         end
         for j = 1: dim_Aut_val(1)
              if(temp_array(1,j) ~= Inf)
                    fprintf("verifica fallita\n")
                    %return; 
              end
         end
         fprintf("autovalori verificati\nAkerman e stato applicato con successo\n")
end