function [K,S,CLP] = ProgettoOsservatoreOTTIMO(A,B,C, rho,gamma)
    %% DEFINIZIONE del sistema DUALE
    A_d = A';
    B_d = C';
    %C_d = B';
    dim_A = size(A);
    dim_B = size(B);
    R = rho*eye(dim_B(2)) ; % matrice del peso dell ingresso 
    Q_lqr = eye(dim_A(2)); % matrice del peso dello stato 
    A_d_mod = A_d + gamma*eye(dim_A(2));
    [K,S,CLP] = lqr(A_d_mod,B_d,Q_lqr,R);
    
end
  