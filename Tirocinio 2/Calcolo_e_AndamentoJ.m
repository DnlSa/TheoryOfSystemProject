function J_val = Calcolo_e_AndamentoJ(u_star_T,time,t_bar,u0,u_T_2,T_2)
                                    
syms tau
m =size(time); 
m= m(2); % aggiorno il numero di elementi componente il vettore di tempi
% Calcolo indice di costo tramite int(u(tau)*u(tau)',tau,0,Ti))
J_val = zeros(1,(m));
for k= 1:(m)  
    u= u_star_T(k,1)*(u_star_T(k,1))';
    J_val(:,k) = int(u,tau,0,time(1,k)); % il calcolo di questo integrale impiega tanto tempo
    if(time(1,k)==t_bar)
        u_star_compare = u_star_T(k,1);
    end
end

% calcolo indice di costo del controllo passato in input dall utente 
u_0=  u0*u0'; % norma 2 al quadrato del vettore di controllo noto
J_u0 = int(u_0,tau,0,t_bar) % calcolo indice di costo di controllo e tempo noto 

%%
if(T_2 ~= 0)
    u_1=  u_T_2*u_T_2'; % norma 2 al quadrato del vettore di controllo noto
    J_u1 = round(int(u_1,tau,0,T_2),4) % calcolo indice di costo di controllo e tempo noto 
    figure(5)
    grid on 
    hold on 
    plot(time,J_val,"r-*");
    plot(t_bar,J_u0,"b-*");
    plot(T_2,J_u1,"m-*");
    title("andamento dell'indice di costo J_Ti");
    xlabel("tempo t");
    ylabel("indice di costo J");
    legend("JTi(u*)","J_tbar(u0)");
    hold off
else
   % grafico dell andamento J
    figure(5)
    grid on 
    hold on 
    plot(time,J_val,"r-*");
    plot(t_bar,J_u0,"b-*");
    title("andamento dell'indice di costo J_Ti");
    xlabel("tempo t");
    ylabel("indice di costo J");
    legend("JTi(u*)","J_tbar(u0)");
    hold off
 end
end
