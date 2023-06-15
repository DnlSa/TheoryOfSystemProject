 
% una funzione che grafica l andamento dei vettori di controllo 

function rit = Plot_Andamento_U_X(A,B,v_max,x)
    
  [~ , cols_A] = size(A); % dimensione della matrice A

  if(cols_A ~= 2)
     fprintf("le colonne di A come le righe devono avere dimensioene 2\n");
     return
  end

  u = TrovaControllo(A,B,v_max,x);
  if(u==-1 )
     rit = -1;
     return
  else 
      rit = 1;
  end
  u = horzcat(u,0); % devo impostare questo vettore di controllo cosi altrimenti mi da errore 
  v_1 = 0:1:v_max; 
  % adesso vediamo l'andamento di x(.) tra [0 , v]
  X = zeros(cols_A,1);
  x_f = X ; % conterrà tutti i stati x raggiunti e verra graficato 
  for(i=1:1:v_max)
     X = (A*X)+(B*u(:,i)); % calcolo vettori di stato 
     x_f = horzcat(x_f,X); % li concateno dentro un vettore x_f
  end
  x_f
  size(x_f)
  x1 = x_f(1,:);
  x2 = x_f(2,:);
  % adesso verrà costruito il grafico dell andamento di stato x
  figure(1);

  plot(v_1,u,"b:o");
  legend("u(.) -> vettore di controllo ");
  grid on
  
  figure(5);
    hold on
  plot(v_1,x1,"r-x");
   plot(v_1,x2,"g--*");
  title ("Andamento di u(.) e di x(.)");
  xlabel("tempo discreto v");
  grid on
  legend("x1 componente_1 " ,"x2 omponente_2 ");
  hold off

  figure(2)
  hold on
  grid on
  rad = 0:pi/64:pi;
  plot(cos(rad), sin(rad), 'k');
  plot(x1,x2,"b--*");
  xlabel("Coordinata X");
  ylabel("Coordinata Y");
  title ("andamento del vettore di stato del sistema ");
  legend("c.unitaria","coordiante dello stato x finale");
  hold off
end