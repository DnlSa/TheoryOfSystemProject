 function Plot_Andamento_J(A,B,v_max,x)
    v = 2:v_max; % vettore contenete i valori di v (dato dal testo dell esercizio il range)
    J =  IndiceCosto(A,B,v_max,x);% inserisco tutti i costi in un vettore J
    if(J == -1)
        return; 
    end
    figure(3)
    plot(v,J,"b-o"); % ne faccio il plot vedere le stampe
    title("Andamento J"); % titolo grafico 
    xlabel("tempo discreto v[s]"); % etichetta su asse x
    ylabel("indice di costo");% etichetta su asse y
      grid on
end