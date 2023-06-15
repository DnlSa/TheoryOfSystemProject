function fig = StabilityTest (Nsys,Dsys,Ncom_1,Dcom_1,r_s_1,fig)
 Plant = tf(Nsys,Dsys); % funzione di trasferimento del sistema fisico  
   r=r_s_1;
   Controller = tf(Ncom_1,Dcom_1); % controllore sistema
   L = Plant*Controller; % funzione d anello 
   Wyr_0 = L/(1+L);
   fig= fig+1;
   figure(fig)
    rlocus(L)
    grid on 
    fig= fig+1;
   figure(fig)
    margin(L)
    grid on 
    legend
    fig= fig+1;
   figure(fig)
    nyquist(L)
    grid on 
    fig= fig+1;
   figure(fig) % verifica prestazioni tramite la funzione step
    hold on 
    step(Wyr_0)
    grid on 
    legend
end