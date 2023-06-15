%%
% funzionei indice di costo 
% di seguito riportate le formule per il calcolo 
% scritte nei 3 modi equivalenti fra loro 
% sum(from h=0 to v-1) {(||u(h)||2)^2}
% sum(from h=0 to v-1) {u(h)'*u(h)} -> useremo questa forma
% w'*w dove ( w = P_bar'*beta)

function Jv = IndiceCosto(A,B,v_max,x)
    % dobbiamo calcolare il controllo ottimale u(.) per ogni indice v fino
    % a v_max .
    if (v_max<2) % conrollo sul parametro v_max passato 
        disp("l indice v deve avere almeno un valore >=2")
        Jv = -1;
        return
    end

    Jv=[]; % dichiaro una matrice vuota
    for  v = 2:1:v_max %calcolo i vettori di controllo per ogni passo 
        u = TrovaControllo(A,B,v,x);
        if(u == -1) % controllo sul ritorno di TrovaControllo 
            Jv = -1;
            return
        end
      
        val_J = 0; % indice di costo 
        for h = 1:1:v    % sum(from h=0 to v-1)
            ret = (u(:,h)')*(u(:,h)); %{u(h)'*u(h)}
            val_J = val_J+ret ;     
        end
        Jv = horzcat(Jv , val_J); % Jv e un array il cui ritorno verra adottato per fare il plot 
        
    end
    fprintf("Jv per v = 2...%d : %f\n",v_max,val_J);
end

