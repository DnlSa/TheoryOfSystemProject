function time = timevector(ts,max_time,t_bar)

% calcola vettore dei tempi 

% calcola tempi consecutivi 
time_t = ts:ts:max_time; % creo vettore dei tempi
m =size(time_t); 
m= m(2); % estrapolo l'indice che mi serve 
time=[];
j=1;
for k= 1:1:m
    if(t_bar==time_t(1,k)&& j==1)
        time = horzcat(time, time_t(1,k));
        j=0;
    elseif(t_bar<time_t(1,k)&& j==1)
        time=horzcat(time,t_bar);
        time=horzcat(time,time_t(1,k));
        j=0;
    else
        time = horzcat(time, time_t(1,k));
    end
end
end
