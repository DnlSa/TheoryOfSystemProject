function AndamentoU(u_star_T,time,t_bar,u0,ts_1,opt)

syms tau
m =size(time); 
m= m(2); % aggiorno il numero di elementi componente il vettore di tempi

% Grafico dei controlli ottimi trovati 
% per graficare il controllo abbiamo bisogno di valutare u*Ti(.) che e
% espressa nella variabile tau
% verra simulato un grafico per veder
 u_valued = [];
 u_valued_compare = [];
if(opt == 1)    
    for k= 1:m % scorre il vettore dei tempi  
         legendInfo{k} = ['time(' num2str(time(1,k)) '[sec])'];
        for j = 0:ts_1:time(1,k) % vettore di tempi i-esimo simulato a passi di ts_1 
            u_val = subs(u_star_T(k,1),tau,j);
            if(time(1,k) == t_bar) % prelevo il controllo ottimo di t_bar
                u_valued_compare = horzcat(u_valued,u_val);
            end
            u_valued = horzcat(u_valued,u_val);
        
        end
        figure(1)
        grid on 
        hold on 
        time_1 =0:ts_1:time(1,k);
        plot(time_1,u_valued);
        xline(time(1,k));
        title("andamento di u*(.)");
        xlabel("tempo t");
        legend(legendInfo);
        u_valued=[];
    end
else
     for k= 1:m % scorre il vettore dei tempi   
        if(time(1,k) == t_bar) % prelevo il controllo ottimo di t_bar
            for j = 0:ts_1:time(1,k) % vettore di tempi i-esimo simulato a passi di ts 
                u_val = subs(u_star_T(k,1),tau,j);
                u_valued_compare = horzcat(u_valued_compare,u_val);
            end
            break;
        end
     end
end
time_1 = 0;

%graficazione del vettore di controllo passato dall utente 
for j = 0:ts_1:t_bar % vettore di tempi i-esimo simulato a passi di ts 
   u_val = subs(u0,tau,j);
   u_valued = horzcat(u_valued,u_val);
end
 figure(2)
 grid on 
 hold on
 time_1 =0:ts_1:t_bar;
 plot(time_1,u_valued);
 plot(time_1,u_valued_compare);
 title("confronto andamento di u0(.) con u*(.)");
 xlabel("tempo t");
 legend("controllo NON ottimo u0(.)","Controllo ottimo u*(.)");
 hold off
end