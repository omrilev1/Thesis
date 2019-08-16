function [Estimated_Bhattacharyya,SER_Vec,Is_Frozen_Bit_Index_Vec,Estimated_Is_Frozen_Bit_Index_Vec,...
    Estimated_LLR_Avg_Vec,Estimated_LLR_Var_Vec,Estimated_SER_Vec,Simulated_BLER,Simulated_BER,...
    Simulated_First_Symbol_Vec,Simulated_LLR_Avg_Vec,Simulated_LLR_Var_Vec] = listDecoderSim_GeneralChannel(NumOfWorkers,Is_Genie_Aided,g,g0,Decode_Type,L,GA_CRC_Length,...
    Constellation_Type,m,m_BIPCM,m_SCL,Constellation_Mapping_Array,N,R,...
    BLER,capacity,Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec,channelType)

Num_Of_Encoded_Symbols_In_Constellation = m;

Max_Num_Of_Iterations = floor((Min_Num_Of_Errors_Sufficient_Statistics/BLER));
Num_Of_Iterations_Print_State = 100;
Min_Num_Of_Iterations = 99;

U_N = N;

Bit_Level_U_N = N/Num_Of_Encoded_Symbols_In_Constellation;

Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
Multilevel_Bit_Level_U_N = Constellation_2D_Dimension*Bit_Level_U_N;

U = zeros(NumOfWorkers,N);
X_after_encoding = zeros(NumOfWorkers,N);
LLR = zeros(NumOfWorkers,N);
LL0 = zeros(NumOfWorkers,N);
LL1 = zeros(NumOfWorkers,N);
X_after_constellation = zeros(NumOfWorkers,Multilevel_Bit_Level_U_N);
Y = zeros(NumOfWorkers,Multilevel_Bit_Level_U_N);
Estimated_U = zeros(NumOfWorkers,N);
Estimated_L = zeros(NumOfWorkers,U_N);
Iter_Vec = zeros(1,NumOfWorkers);
temp_Errors_Vec = zeros(NumOfWorkers,N);
Errors_Vec = zeros(NumOfWorkers,N);
Num_Of_Sufficient_Statistics_Indexes = zeros(1,NumOfWorkers);
Num_Of_Sufficient_Statistics_Indexes_To_Continue = zeros(1,NumOfWorkers);
BLER_Errors_Vec = zeros(NumOfWorkers,1);
First_Bit = zeros(1,NumOfWorkers);
First_Bit_Vec = zeros(NumOfWorkers,N);
First_Symbol = zeros(1,NumOfWorkers);
First_Symbol_Vec = zeros(NumOfWorkers,U_N);
temp_LLR_Vec = zeros(NumOfWorkers,U_N);
LLR_Divide_Vec = zeros(NumOfWorkers,U_N);
LLR_Avg_Vec = zeros(NumOfWorkers,U_N);
LLR_Var_Vec = zeros(NumOfWorkers,U_N);
Bhattacharyya_Vec = zeros(NumOfWorkers,U_N);

% for ParIter=1:NumOfWorkers
parfor ParIter=1:NumOfWorkers
    
    for Iter=1:1:floor(Max_Num_Of_Iterations/NumOfWorkers)
        
        Iter_Vec(1,ParIter) = Iter_Vec(1,ParIter) + 1;
        
        %-----------------------------generate input-----------------------------%
        %Generate random input
        
        U(ParIter,:) = randi([0,1],1,N);
        
        %---------------------------------encode---------------------------------%
        %Encode input into output
        X_after_encoding(ParIter,:) = mod(U(ParIter,:)*g,2);
        
        
        
        %--------------------constellation Channel and Decode--------------------%
        % channel model
        [currLLR,~,~] = channelModel((-1).^(X_after_encoding(ParIter,:)),capacity,channelType);
%         [LLR(ParIter,:),LL0(ParIter,:),LL1(ParIter,:)] = channelModel((-1).^(X_after_encoding(ParIter,:)),capacity,channelType);
        LLR(ParIter,:) = currLLR;
        
        
        % Decode
        if(L==1)
            [Estimated_U(ParIter,:),~,Estimated_L(ParIter,:)] = gArikan_SC_Decoder(LLR(ParIter,:),U(ParIter,:),N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
        else
            if(Is_Genie_Aided)
                [Estimated_U(ParIter,:),~,Estimated_L(ParIter,:)] = gArikan_SC_Decoder(LLR(ParIter,:),U(ParIter,:),N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
            else
                [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_SCL_Decoder_TopNew(L,GA_CRC_Length,LL0(ParIter,:),LL1(ParIter,:),...
                    U(ParIter,:),Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
            end
        end
        
        if(Is_Genie_Aided)
            
            temp_Errors_Vec(ParIter,:) = (Estimated_U(ParIter,:)~=U(ParIter,:));
            
            Errors_Vec(ParIter,:) = Errors_Vec(ParIter,:) + temp_Errors_Vec(ParIter,:);
            
            Num_Of_Sufficient_Statistics_Indexes(1,ParIter) = sum(Errors_Vec(ParIter,:)>Min_Num_Of_Errors_Sufficient_Statistics/NumOfWorkers);
            
            if(mod(Iter_Vec(1,ParIter),Num_Of_Iterations_Print_State)==0)
                sprintf('worker %d finish iteration %d (out of maximal %d iterations), %d indexes (out of %d needed) has already sufficient statistics', ParIter,Iter_Vec(1,ParIter),floor(Max_Num_Of_Iterations/NumOfWorkers),Num_Of_Sufficient_Statistics_Indexes(1,ParIter),N*(1-R))
            end
            
            if(Num_Of_Sufficient_Statistics_Indexes(1,ParIter)>=N*(1-R)+(GA_CRC_Length>0)*GA_CRC_Length)
                break;
            end
            
            if(Iter_Vec(1,ParIter)>max(floor((Max_Num_Of_Iterations/NumOfWorkers)/10),Min_Num_Of_Iterations))
                Num_Of_Sufficient_Statistics_Indexes_To_Continue(1,ParIter) = sum(Errors_Vec(ParIter,:)>floor((Min_Num_Of_Errors_Sufficient_Statistics/NumOfWorkers)/((Max_Num_Of_Iterations/Iter_Vec(1,ParIter))/10)));
                if(Num_Of_Sufficient_Statistics_Indexes_To_Continue(1,ParIter)<N*(1-R)+(GA_CRC_Length>0)*GA_CRC_Length)
                    break;
                end
            end
            
            temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:).*(-1).^U(ParIter,:);
            
            
            LLR_Avg_Vec(ParIter,:) = LLR_Avg_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
            LLR_Var_Vec(ParIter,:) =  LLR_Var_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:).^2;
            
            Bhattacharyya_Vec(ParIter,:) = Bhattacharyya_Vec(ParIter,:) + sqrt(exp(-temp_LLR_Vec(ParIter,:))); %???? ???? ????? ?? ??????
            
        else %simulation - not genie aided
            
            temp_Errors_Vec(ParIter,:) = Estimated_U(ParIter,:)~=U(ParIter,:);
            
            BLER_Errors_Vec(ParIter,:) = BLER_Errors_Vec(ParIter,:) + (sum((temp_Errors_Vec(ParIter,:))>0)>0);
            
            Errors_Vec(ParIter,:) = Errors_Vec(ParIter,:) + (temp_Errors_Vec(ParIter,:));
            
            try
                
                First_Bit(1,ParIter) = find(temp_Errors_Vec(ParIter,:),1,'first');
                First_Bit_Vec(ParIter,:) = First_Bit_Vec(ParIter,:) + temp_Errors_Vec(ParIter,:).*[zeros(1,First_Bit(1,ParIter)-1),1,zeros(1,N-First_Bit(1,ParIter))];
                
                First_Symbol(1,ParIter) = First_Bit(1,ParIter);
                First_Symbol_Vec(ParIter,:) = First_Bit_Vec(ParIter,:);
                
                temp_LLR_Vec(ParIter,:) = [ones(1,First_Symbol(1,ParIter)),zeros(1,U_N-First_Symbol(1,ParIter))];
                LLR_Divide_Vec(ParIter,:) = LLR_Divide_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
                
                temp_LLR_Vec(ParIter,:) = (Estimated_L(ParIter,:).*(-1).^U(ParIter,:)).*temp_LLR_Vec(ParIter,:);
                
                
                LLR_Avg_Vec(ParIter,:) = LLR_Avg_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
                LLR_Var_Vec(ParIter,:) =  LLR_Var_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:).^2;
                
            catch
                
                LLR_Divide_Vec(ParIter,:) = LLR_Divide_Vec(ParIter,:) + 1;
                
                temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:).*(-1).^U(ParIter,:);
                
                
                LLR_Avg_Vec(ParIter,:) = LLR_Avg_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
                LLR_Var_Vec(ParIter,:) =  LLR_Var_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:).^2;
                
            end
            
            if(mod(Iter,Num_Of_Iterations_Print_State)==0)
                sprintf('worker %d finish iteration %d (out of maximal %d iterations), %d errors (out of %d needed) has already occurred', ParIter,Iter,floor(Max_Num_Of_Iterations/NumOfWorkers),BLER_Errors_Vec(ParIter,:),floor(Min_Num_Of_Errors_Sufficient_Statistics/NumOfWorkers))
            end
            
            if(BLER_Errors_Vec(ParIter,:)>=Min_Num_Of_Errors_Sufficient_Statistics/NumOfWorkers)
                break;
            end
            
            if(Iter_Vec(1,ParIter)>max(floor((Max_Num_Of_Iterations/NumOfWorkers)/10),Min_Num_Of_Iterations))
                if(BLER_Errors_Vec(ParIter,:)<floor((Min_Num_Of_Errors_Sufficient_Statistics/NumOfWorkers)/((Max_Num_Of_Iterations/Iter_Vec(1,ParIter))/10)))
                    break;
                end
            end
            
        end
        
    end
    
end

if(Is_Genie_Aided)
    
    if(NumOfWorkers==1)
        Num_Of_Errors = Errors_Vec;
        SER_Vec = Errors_Vec./Iter_Vec;
        Estimated_Bhattacharyya = Bhattacharyya_Vec./Iter_Vec;
    else
        for i=1:1:NumOfWorkers
            SER_Vec(i,:) = Errors_Vec(i,:)./sum(Iter_Vec);
            Num_Of_Errors(i,:) = Errors_Vec(i,:);
            Estimated_Bhattacharyya(i,:) = Bhattacharyya_Vec(i,:)./sum(Iter_Vec);
        end
        SER_Vec = sum(SER_Vec);
        Num_Of_Errors = sum(Num_Of_Errors);
        Estimated_Bhattacharyya = sum(Estimated_Bhattacharyya);
    end
    
    Estimated_Bhattacharyya(Estimated_Bhattacharyya<0)=0;
    
    SER_Vec = min(0.5,SER_Vec);
    [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec));
    Sorted_SER_Index_Vec = size(SER_Vec,2) - Sorted_SER_Index_Vec + 1;
    Is_Frozen_Bit_Index_Vec = zeros(1,N);
    for Index=Sorted_SER_Index_Vec(1,round(R*N)+1+(GA_CRC_Length>0)*GA_CRC_Length:N)
        Is_Frozen_Bit_Index_Vec(1,Index) = 1;
    end
    Is_Frozen_Bit_Index_Vec(1,:) = Is_Frozen_Bit_Index_Vec(1,:) == ones(1,N);
    
    Total_Number_Of_Iter = sum(Iter_Vec);
    
    %         Estimated_LLR_Avg_Vec = sum(LLR_Avg_Vec,1)./(Total_Number_Of_Iter-Num_Of_Errors);
    %         Estimated_LLR_Var_Vec = ((Total_Number_Of_Iter-Num_Of_Errors)/(Total_Number_Of_Iter-Num_Of_Errors-1))*(sum(LLR_Var_Vec,1)./(Total_Number_Of_Iter-Num_Of_Errors)-Estimated_LLR_Avg_Vec.^2);
    Estimated_LLR_Avg_Vec = sum(LLR_Avg_Vec,1)./(Total_Number_Of_Iter);
    Estimated_LLR_Var_Vec = ((Total_Number_Of_Iter)/(Total_Number_Of_Iter-1))*(sum(LLR_Var_Vec,1)./(Total_Number_Of_Iter)-Estimated_LLR_Avg_Vec.^2);
    temp_Index_LLR = Estimated_LLR_Avg_Vec>0;
    Estimated_SER_Vec = zeros(1,U_N);
    Estimated_SER_Vec(temp_Index_LLR) = 0.5.*(1+erf((-Estimated_LLR_Avg_Vec(temp_Index_LLR)./sqrt(Estimated_LLR_Var_Vec(temp_Index_LLR))./sqrt(2))));
    
    Estimated_SER_Vec(Estimated_LLR_Avg_Vec<=0) = 0.5;
    [Sorted_Estimated_SER_Vec,Sorted_Estimated_SER_Index_Vec] = sort(wrev(Estimated_SER_Vec));
    Sorted_Estimated_SER_Index_Vec = size(Estimated_SER_Vec,2) - Sorted_Estimated_SER_Index_Vec + 1;
    Estimated_Is_Frozen_Bit_Index_Vec = zeros(1,U_N);
    for Index=Sorted_Estimated_SER_Index_Vec(1,round(R*U_N)+1+(GA_CRC_Length>0)*GA_CRC_Length:U_N)
        Estimated_Is_Frozen_Bit_Index_Vec(1,Index) = 1;
    end
    Estimated_Is_Frozen_Bit_Index_Vec(1,:) = Estimated_Is_Frozen_Bit_Index_Vec(1,:) == ones(1,U_N);
    Estimated_Is_Frozen_Bit_Index_Vec = reshape([Estimated_Is_Frozen_Bit_Index_Vec;Estimated_Is_Frozen_Bit_Index_Vec],1,[]);
    
    BLER_lower_bound = max(Sorted_SER_Vec(1:round(R*U_N)));
    BLER_upper_bound = min(sum(Sorted_SER_Vec(1:round(R*U_N))),0.5);
    
    Simulated_BLER = zeros(0);
    Simulated_BER = zeros(0);
    %         Simulated_First_Bit_Vec = zeros(0);
    Simulated_First_Symbol_Vec = zeros(0);
    Simulated_LLR_Avg_Vec = zeros(0);
    Simulated_LLR_Var_Vec = zeros(0);
    
else
    
    if(NumOfWorkers==1)
        Simulated_BLER = BLER_Errors_Vec./Iter_Vec;
        Simulated_BER = sum(Errors_Vec)./(Iter_Vec*round(N*R));
        %             Simulated_First_Bit_Vec = First_Bit_Vec./Iter_Vec;
        Simulated_First_Symbol_Vec = First_Symbol_Vec./Iter_Vec;
        Simulated_LLR_Avg_Vec = LLR_Avg_Vec./LLR_Divide_Vec;
        Simulated_LLR_Var_Vec = LLR_Var_Vec./LLR_Divide_Vec-Simulated_LLR_Avg_Vec.^2;
    else
        BLER_Vec = zeros(1,NumOfWorkers);
        BER_Vec = zeros(1,NumOfWorkers);
        for i=1:1:NumOfWorkers
            BLER_Vec(i) = BLER_Errors_Vec(i)./sum(Iter_Vec);
            BER_Vec(i) = sum(Errors_Vec(i,:))./(sum(Iter_Vec)*round(N*R));
        end
        Simulated_BLER = sum(BLER_Vec);
        Simulated_BER = sum(BER_Vec);
        %             Simulated_First_Bit_Vec = sum(First_Bit_Vec)./sum(Iter_Vec);
        Simulated_First_Symbol_Vec = sum(First_Symbol_Vec)./sum(Iter_Vec);
        Simulated_LLR_Avg_Vec = sum(LLR_Avg_Vec)./sum(LLR_Divide_Vec);
        Simulated_LLR_Var_Vec = sum(LLR_Var_Vec)./sum(LLR_Divide_Vec)-Simulated_LLR_Avg_Vec.^2;
    end
    
    SER_Vec = zeros(0);
    Is_Frozen_Bit_Index_Vec = zeros(0);
    Estimated_Is_Frozen_Bit_Index_Vec = zeros(0);
    Estimated_LLR_Avg_Vec = zeros(0);
    Estimated_LLR_Var_Vec = zeros(0);
    Estimated_SER_Vec = zeros(0);
    Estimated_Bhattacharyya = zeros(0);
    
end

end

