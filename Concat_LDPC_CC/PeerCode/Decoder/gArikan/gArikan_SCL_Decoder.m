function [Estimated_U,Estimated_X,Estimated_L0,Estimated_L1,Indexes] = gArikan_SCL_Decoder(LL0,LL1,L,U,N,Is_Frozen_Bit_Index_Vec)

%     LL_Matrix = [LL0;LL1];
%     max_LL = max(LL_Matrix(:));
%     
%     LL0 = LL0 - max_LL;
%     LL1 = LL1 - max_LL;

    if(isempty(find(~Is_Frozen_Bit_Index_Vec,1)))
        
        [Estimated_U,Estimated_X,Estimated_L0,Estimated_L1,Indexes] = gArikan_SCL_Decoder_Frozen(LL0,LL1,U,N);
        
    elseif(N==2)
        
        [Estimated_U,Estimated_X,Estimated_L0,Estimated_L1,Indexes] = gArikan_SCL_Decoder_2(LL0,LL1,L,U,Is_Frozen_Bit_Index_Vec);
                
    else
        
        next_N = N/2;
        
        next_L = min(L,(2^sum(~Is_Frozen_Bit_Index_Vec))*size(LL0,1));
                                           
        Estimated_X = zeros(next_L,N);

        temp_odd_0 = LL0(:,1:2:end);
        temp_even_0 = LL0(:,2:2:end);
        temp_odd_1 = LL1(:,1:2:end);
        temp_even_1 = LL1(:,2:2:end);
        
        temp_sum_00 = temp_odd_0 + temp_even_0;
        temp_sum_11 = temp_odd_1 + temp_even_1;

        temp_sum_01 = temp_odd_0 + temp_even_1;
        temp_sum_10 = temp_odd_1 + temp_even_0;

        L0 = Add_LL(temp_sum_00,temp_sum_11);
        L1 = Add_LL(temp_sum_01,temp_sum_10);
        
        [Estimated_U_1,Estimated_X_1,Estimated_L0_1,Estimated_L1_1,Indexes_1] = gArikan_SCL_Decoder(L0,L1,L,U(1:next_N),next_N,Is_Frozen_Bit_Index_Vec(1:next_N));       
        
        temp_odd_0 = temp_odd_0(Indexes_1,:);
        temp_even_0 = temp_even_0(Indexes_1,:);
        temp_odd_1 = temp_odd_1(Indexes_1,:);
        temp_even_1 = temp_even_1(Indexes_1,:);
        
        [L0,L1] = Calc_LL(Estimated_X_1==1,temp_odd_0,temp_odd_1,temp_even_0,temp_even_1);
        
        [Estimated_U_2,Estimated_X_2,Estimated_L0_2,Estimated_L1_2,Indexes_2] = gArikan_SCL_Decoder(L0,L1,L,U(next_N+1:N),next_N,Is_Frozen_Bit_Index_Vec(next_N+1:N));
        
        Indexes = Indexes_1(Indexes_2);
        
        Estimated_U = [Estimated_U_1(Indexes_2,:),Estimated_U_2];
        Estimated_L0 = [Estimated_L0_1(Indexes_2,:),Estimated_L0_2];
        Estimated_L1 = [Estimated_L1_1(Indexes_2,:),Estimated_L1_2];
        
        Estimated_X(:,1:2:N) = mod(Estimated_X_1(Indexes_2,:)+Estimated_X_2,2);
        Estimated_X(:,2:2:N) = Estimated_X_2;
         
    end 
    
end

function [Estimated_U,Estimated_X,Estimated_L0,Estimated_L1,Indexes] = gArikan_SCL_Decoder_Frozen(LL0,LL1,U,N)

    g_Arikan = [1,0;1,1];
    g = g_Arikan;
    for i=1:1:log2(N)-1
        g = kron(g,g_Arikan);
    end
    g = bitrevorder(g);

    temp_LL_length = size(LL0,1);
    
    Estimated_U = ones(temp_LL_length,1) * U;
    Estimated_X = ones(temp_LL_length,1) * mod(U*g,2);
    Estimated_L0 = LL0;
    Estimated_L1 = LL1;
    Indexes = (1:temp_LL_length)';
    
end

function [Estimated_U,Estimated_X,Estimated_L0,Estimated_L1,Indexes] = gArikan_SCL_Decoder_2(LL0,LL1,L,U,Is_Frozen_Bit_Index_Vec)

    temp_LL_length = size(LL0,1);
        
    temp_odd_0 = LL0(:,1);
    temp_even_0 = LL0(:,2);
    temp_odd_1 = LL1(:,1);
    temp_even_1 = LL1(:,2);

    temp_sum_00 = temp_odd_0 + temp_even_0;
    temp_sum_11 = temp_odd_1 + temp_even_1;

    temp_sum_01 = temp_odd_0 + temp_even_1;
    temp_sum_10 = temp_odd_1 + temp_even_0;

    L0 = Add_LL(temp_sum_00,temp_sum_11);
    L1 = Add_LL(temp_sum_01,temp_sum_10);

    [Estimated_U_1,Indexes_1,Estimated_L0_1,Estimated_L1_1] = Find_Best_Paths(U(1),L,L0,L1,Is_Frozen_Bit_Index_Vec(1),temp_LL_length);
    
    temp_LL_length = size(Estimated_U_1,1);

    temp_odd_0 = temp_odd_0(Indexes_1,:);
    temp_even_0 = temp_even_0(Indexes_1,:);
    temp_odd_1 = temp_odd_1(Indexes_1,:);
    temp_even_1 = temp_even_1(Indexes_1,:);

    [L0,L1] = Calc_LL(Estimated_U_1==1,temp_odd_0,temp_odd_1,temp_even_0,temp_even_1);

    [Estimated_U_2,Indexes_2,Estimated_L0_2,Estimated_L1_2] = Find_Best_Paths(U(2),L,L0,L1,Is_Frozen_Bit_Index_Vec(2),temp_LL_length);

    Indexes = Indexes_1(Indexes_2);

    Estimated_U = [Estimated_U_1(Indexes_2),Estimated_U_2];

    Estimated_X = [mod(Estimated_U_1(Indexes_2)+Estimated_U_2,2),Estimated_U_2];
    
    Estimated_L0 = [Estimated_L0_1(Indexes_2),Estimated_L0_2];
    
    Estimated_L1 = [Estimated_L1_1(Indexes_2),Estimated_L1_2];

end

function [Estimated_U,Indexes,Estimated_L0,Estimated_L1] = Find_Best_Paths(U,L,L0,L1,Is_Frozen_Bit,temp_LL_length)

    if(Is_Frozen_Bit)
        
        Estimated_L0 = L0;
        Estimated_L1 = L1;

        Indexes = (1:temp_LL_length)';

        Estimated_U = ones(temp_LL_length,1) * U;

    else

        temp_LL = [L0;L1];
        temp_next_L = length(temp_LL);

        if(temp_next_L>L)
            best_Indexes = find(temp_LL>=median(temp_LL));
            best_Indexes = best_Indexes(1:L);
        else
            best_Indexes = (1:temp_next_L)';
        end

        Indexes = mod((best_Indexes-1),temp_LL_length)+1;

        Estimated_L0 = L0(Indexes);
        Estimated_L1 = L1(Indexes);

        Estimated_U = (best_Indexes > temp_LL_length);

    end

end

function [L] = Add_LL(temp_sum_0,temp_sum_1)

    L = (max(temp_sum_0,temp_sum_1) + log(1+exp(-abs(temp_sum_0 - temp_sum_1))));
    
%     temp = log(1+exp(-abs(temp_sum_0 - temp_sum_1)));
%     temp(isnan(temp)) = log(2);
%     L = max(temp_sum_0,temp_sum_1) + temp;

end

function [L0,L1] = Calc_LL(Estimated,odd_0,odd_1,even_0,even_1)

    L0 = bsxfun(@times,Estimated,odd_1)+bsxfun(@times,~Estimated,odd_0)+even_0;
    L1 = bsxfun(@times,Estimated,odd_0)+bsxfun(@times,~Estimated,odd_1)+even_1;

%     L0 = zeros(size(Estimated));
%     L1 = zeros(size(Estimated));
% 
%     L0(Estimated) = odd_1(Estimated);
%     L0(~Estimated) = odd_0(~Estimated);
%     L0 = L0 + even_0;
%     
%     L1(Estimated) = odd_0(Estimated);
%     L1(~Estimated) = odd_1(~Estimated);
%     L1 = L1 + even_1;

end