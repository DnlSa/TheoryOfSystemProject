%% implementazione andamento errore di stima
% si assume che la stima dell errore per lo stato iniziale sia 0 
% QUINDI e(0) = x(0)
% e(t) = exp^((A-V*C)*t)*e(0)
% se A-VC = H presenta tutti poli minori di 0 questi genereranno dei modi 
% convergenti a 0 per lim(t->inf)e(t)=0
function AndamentoErroreDiStima(AVC,e0,ts,max_time_simulation)
    syms s t
    time = 0:ts:max_time_simulation; %tempo di simulazione  
    dim_AVC = size(AVC); %calcolo dimensione AVC
    Id_mat  = eye(dim_AVC(1)); 
    
    temp_mat = inv((s*Id_mat-AVC)); % (sI-(A-V*C))^(-1)
    e=ilaplace(temp_mat)*e0; %antitrasformata di laplace per dar vita a  exp^((A-V*C)*t)*e(0)
    
    e1 = []; % 4 stime in quanto la matrice (A-V*C) e una matrice 4*4 
    e2 = [];
    e3 = [];
    e4 = [];
    for k = 0:ts:max_time_simulation
        temp = subs(e,t,k);
         e1 = horzcat(e1,temp(1));
         e2 = horzcat(e2,temp(2));
         e3 = horzcat(e3,temp(3));
         e4 = horzcat(e4,temp(4));
    end
    hold on 
    grid on 
    plot(time,e1,"b-");
    plot(time,e2,"g-");
    plot(time,e3,"y-");
    plot(time,e4,"r-");
    legend("e1","e2","e3","e4"); 
    hold off
end
