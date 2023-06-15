function StatiUscita(A,B,C,D,t,u,x0)

sys = ss(A,B,C,D); 
[y,T,x] = lsim(sys,u,t,x0);

figure(5)
title("Andamento ingresso-uscita")
plot(t,x);
grid on 
legend("posCarr1","posCarr2","velCarr1","velCarr2");

figure(6)
title("Andamento stati del sistema ")
hold on ;
plot(t,y);
plot(t,u);
hold off;
grid on 
legend("uscita y(t)","ingresso u(t)");