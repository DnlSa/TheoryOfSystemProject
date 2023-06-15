function aut_val = CalcoloAutovaloriDesiderati(A,Aut_val_distance,aut_val_des,offset_autval,all_real)
    Aut_val = eig(A);
    dim_Aut_val = size(Aut_val);
    dim_aut_val_des = size(aut_val_des);
    if( offset_autval == 1 && dim_Aut_val(1) >= dim_aut_val_des(2)) % calcola autovalori considerando l'offset
        aut_val  = [];
        for k = 1:dim_Aut_val(1)
            temp = round(Aut_val(k,1)+Aut_val_distance,4);
            if(all_real == 1)
                temp = real(temp);
            end
            aut_val = horzcat(aut_val,temp);  % li metto reali perche almeno ho meno fenomeni oscillatori 
            if(temp >= 0)
                if (temp > 0 )
                    fprintf("L'autovalore %d darà luogo ad un sistema INSTABILE\n",temp)
                else
                    fprintf("L'autovalore %d darà luogo ad un sistema STABILE SEMPLICEMENTE\n",temp)
                end
            end
        end
    else    
       if(dim_Aut_val(1) < dim_aut_val_des(2) ) % nel tirocinio 4 non si realizza ma nel 3 e probabile nelle forme di kalman
            for k= 2  : dim_aut_val_des(2)
                if(aut_val_des(:,k-1) <= aut_val_des(:,k))
                    temp  = k; % salvo la posizione del vettore da eliminare   
                end
            end
           aut_val_des(temp) = []; % elimina elemento k esimo piu sensibile
           aut_val = aut_val_des;
       else % inserisco un completamento del vettore di autovalori desiderati
           aut_val = aut_val_des;
            for k = 1:(dim_Aut_val(1)- dim_aut_val_des(2))
                temp = round(Aut_val(k,1)+Aut_val_distance,2);
                if(all_real == 1)
                    temp = real(temp);
                end
                aut_val = horzcat(aut_val,temp); 
                if(temp >= 0)
                    if (temp > 0 )
                        fprintf("L'autovalore %d darà luogo ad un sistema INSTABILE\n",temp)
                    else
                        fprintf("L'autovalore %d darà luogo ad un sistema STABILE SEMPLICEMENTE\n",temp)
                    end
                end  
            end  
       end
    end
end