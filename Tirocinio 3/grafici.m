function grafici(A,aut_val_1,aut_val_2,aut_val_3)
    hold on 
    xline(0);
    yline(0);
    plot(real(aut_val_1), imag(aut_val_1), 'r*','MarkerSize',8);
    aut_val_A = eig(A);
    plot(real(aut_val_A), imag(aut_val_A), 'b*','MarkerSize',8);
    plot(real(aut_val_3), imag(aut_val_3), 'g*','MarkerSize',8);
    plot(real(aut_val_2), imag(aut_val_2), 'c*','MarkerSize',8);
    legend("asse REALE","asse IMMAGINARIO","autovalori OSS. LENTO" , "autovalori ORIGINALI di A", ...
        "autovalori OSS. VELOCE","autovalori OSS. OTTIMO ")
    title('Autovalori su piano Complesso');
    xlabel('Asse reale');
    ylabel('Asse immaginario');
    grid on;
    hold off;
end