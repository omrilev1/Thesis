function [Estimated_U,Estimated_X,Estimated_L] = gArikan_SC_Decoder(Lambda,U,N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)

    if(N==2)
           
        Estimated_U = zeros(1,2);
        
        Estimated_L = zeros(1,2);
                
        Estimated_L(1) = (max(0,Lambda(1)+Lambda(2)) + log(1+exp(-abs(Lambda(1)+Lambda(2))))) - (max(Lambda(1),Lambda(2)) + log(1+exp(-abs(Lambda(1)-Lambda(2)))));
        
        %if(abs(Estimated_L(1))<10^-10)
            
            %Estimated_L(1) = 2*atanh(tanh(Lambda(1)/2)*tanh(Lambda(2)/2));

            %if(Estimated_L(1)==0)

                %Estimated_L(1) = min(Lambda(1),Lambda(2));

            %end

        %end
        
        if(Is_Frozen_Bit_Index_Vec(1))
            Estimated_U(1) = U(1);
        else
            Estimated_U(1) = floor((1-sign(Estimated_L(1)))/2);
        end
        
        if(Is_Genie_Aided)
            Estimated_L(2) = ((-1)^U(1))*Lambda(1)+Lambda(2);
        else
            Estimated_L(2) = ((-1)^Estimated_U(1))*Lambda(1)+Lambda(2);
        end
        
        if(Is_Frozen_Bit_Index_Vec(2))
            Estimated_U(2) = U(2);
        else
            Estimated_U(2) = floor((1-sign(Estimated_L(2)))/2);
        end
        
        if(Is_Genie_Aided)
            Estimated_X = [mod(U(1) + U(2),2),U(2)];
        else
            Estimated_X = [mod(Estimated_U(1) + Estimated_U(2),2),Estimated_U(2)];
        end
        
    else
        
        next_N = N/2;
                           
        Estimated_U = zeros(1,N);
                
        Estimated_X = zeros(1,N);
        
        Estimated_L = zeros(1,N);

        temp_odd = Lambda(1,1:2:end);
        temp_even = Lambda(1,2:2:end);
        temp_sum = temp_odd + temp_even;
        temp_diff = temp_odd - temp_even;
        L = (max(0,temp_sum) + log(1+exp(-abs(temp_sum)))) - (max(temp_odd,temp_even) + log(1+exp(-abs(temp_diff))));
        
        %msk = abs(L)<10^-10;
        %L(msk) = 2.*atanh(tanh(temp_odd(msk)./2).*tanh(temp_even(msk)./2));
        %msk = L == 0;
        %L(msk) = min(temp_odd(msk),temp_even(msk));
        
        [Estimated_U(1,1:next_N),Estimated_X(1:2:end),Estimated_L(1:next_N)] = gArikan_SC_Decoder(L,U(1:next_N),next_N,Is_Frozen_Bit_Index_Vec(1:next_N),Is_Genie_Aided);       

        L = ((-1).^Estimated_X(1:2:end)).*temp_odd+temp_even;

        [Estimated_U(1,next_N+1:N),Estimated_X(2:2:end),Estimated_L(next_N+1:N)] = gArikan_SC_Decoder(L,U(next_N+1:N),next_N,Is_Frozen_Bit_Index_Vec(next_N+1:N),Is_Genie_Aided);
                    
        Estimated_X(1:2:end) = mod(Estimated_X(1:2:end) + Estimated_X(2:2:end),2);
         
    end 
    
end