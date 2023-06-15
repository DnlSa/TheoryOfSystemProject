function AndamentoX(A,B,u_star_T,time,t_bar,u0,ts_1,opt,T_2,u_T_2,u_T_3)

syms tau
m =size(time); 
m= m(2); % aggiorno il numero di elementi componente il vettore di tempi
x_star_T=[];
x_star_T_compare=[];
x_star_T_compare_2=[];

if(opt == 1) % calcola tutti i vettori 
    for k= 1:m % scorre il vettore dei tempi   
        legendInfo{k} = ['time(' num2str(time(1,k)) '[sec])']; % Crea un array di stringhe per la legenda
        for j = 0:ts_1:time(1,k) % vettore di tempi i-esimo simulato a passi di ts 
            temp=  rispostaForzata(A,B,u_star_T(k,1),j);% clcolo risposta forzata 
            temp_1 = double(temp);
            if(time(1,k)== t_bar)
                x_star_T_compare =horzcat(x_star_T,temp_1);
            end
            x_star_T =horzcat(x_star_T,temp_1); % concateno ad un vettore bidimensionale
        end
        % stampa dopo aver fatto la simulazione con il controllo k-esimo per il
        % tempo di simulazione assegnato 
         time_1 =0:ts_1:time(1,k);
         x1 = x_star_T(1,:); % componene 1 del vettore di stato 
         x2 = x_star_T(2,:); % componente 2  del vettore di stato 
     
        figure(3)
        hold on
        grid on 
        xline(time(1,k));
        plot(time_1,x1);
        plot(time_1,x2);
        title("andamento di x*(.)");
        xlabel("tempo t");
        legend("1° componente x*(.)","2° componente x*(.)");
        

        figure(7)
        hold on
        grid on 
        plot(x1,x2,"LineWidth",1);
        title("traiettorie del vettore di stato OTTIMO");
        legend(legendInfo);
        hold off

        % sanifico variabili  per allocare la prossima simulazion
        x1=0 ;
        x2=0 ;
        time_1=0;
        x_star_T = [];
    end
else
    if(T_2>ts_1)
        for j = 0:ts_1:T_2 % vettore di tempi i-esimo simulato a passi di ts
                temp=  rispostaForzata(A,B,u_T_2,j);% clcolo risposta forzata 
                temp_1 = double(temp);
                x_star_T_compare =horzcat(x_star_T_compare,temp_1);
                if(u_T_3 ~= 0 )
                    temp2=  rispostaForzata(A,B,u_T_3,j);% clcolo risposta forzata 
                    temp_2 = double(temp);
                    x_star_T_compare_2 =horzcat(x_star_T_compare_2,temp_2);
                end
        end
    else
     for k= 1:m % scorre il vettore dei tempi
        if(time(1,k)== t_bar)
            for j = 0:ts_1:time(1,k) % vettore di tempi i-esimo simulato a passi di ts 
                temp=  rispostaForzata(A,B,u_star_T(k,1),j);% clcolo risposta forzata 
                temp_1 = double(temp);
                x_star_T_compare =horzcat(x_star_T_compare,temp_1);
            end
        end
     end
    end
end

x_star_T = [];
for j = 0:ts_1:t_bar % vettore di tempi i-esimo simulato a passi di ts 
      temp=  rispostaForzata(A,B,u0,j);% clcolo risposta forzata 
      temp_1 = double(temp);
      x_star_T =horzcat(x_star_T,temp_1); % concateno ad un vettore bidimensionale
end
    time_1 =0:ts_1:t_bar;
    time_2 =0:ts_1:T_2;
    x1 = x_star_T(1,:); % componene 1 del vettore di stato 
    x2 = x_star_T(2,:); % componente 2  del vettore di stato 
    x1_c = x_star_T_compare(1,:); % componene 1 del vettore di stato 
    x2_c = x_star_T_compare(2,:); % componente 2  del vettore di stato
    
    if(u_T_3~=0)
        x1_r = x_star_T_compare_2(1,:); % componene 1 del vettore di stato MAXIMA 
        x2_r = x_star_T_compare_2(2,:); % componente 2  del vettore di stato MAXIMA
    end
    figure(4) % STAMPA LA TRAIETTORIA DDELLO STATO IN BASE AL TEMPO E DIVIDENDO LE COMPONENTI X Y   
    grid on 
    hold on; 
    plot(time_1,x1,"LineWidth",1);
    plot(time_1,x2,"LineWidth",1);
    plot(time_2,x1_c,"LineWidth",1);
    plot(time_2,x2_c,"LineWidth",1);
    legend("1° componente x0(.)","2° componente x0(.)","1° componente x*(.)","2° componente x*(.)");
    title("andamento di confronto con x0(.) e x*(.)");
    xlabel("tempo t");
    hold off;
    
    figure(6) % STAMPA LA TRAIETTORIA DI STATO DEL SISTEMA SU CCORDINATE X Y 
    hold on 
    grid on 
    plot(x1,x2,"LineWidth",1);
    plot(x1_c,x2_c,"LineWidth",1);
    if(u_T_3 ~= 0)
        plot(x1_r,x2_r,"LineWidth",1);
    end
    xlabel("coordinata - X");
    ylabel("coordinata - Y");
    legend("vettore di stato NON ottimo ","Vettore di stato OTTIMO matlab","Vettore di stato OTTIMO maxima");
    title("Traiettoria dello stato x");
    hold off
end