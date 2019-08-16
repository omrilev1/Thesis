 function [Estimated_Bhattacharyya,SER_Vec,Is_Frozen_Bit_Index_Vec,Estimated_Is_Frozen_Bit_Index_Vec,Estimated_LLR_Avg_Vec,Estimated_LLR_Var_Vec,Estimated_SER_Vec,Simulated_BLER,Simulated_BER,Simulated_First_Symbol_Vec,Simulated_LLR_Avg_Vec,Simulated_LLR_Var_Vec] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,Is_Genie_Aided,g,g0,Decode_Type,L,GA_CRC_Length,Constellation_Type,m,m_BIPCM,m_SCL,Constellation_Mapping_Array,N,R,BLER,SNR,Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec)
                                                                                                                                                                                                                    
    M = 2^m;
    
    Num_Of_Encoded_Symbols_In_Constellation = m;
    
    if(M==2)
                
        future_matrix = zeros(0,0);
                
    elseif(mod(log(M)/log(4),1)==0)
            
        if(strcmp(Decode_Type,'Multi_Dimensional_MLPC'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            
        elseif(strcmp(Decode_Type,'Multi_Dimensional_BIPCM'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
           
        elseif(strcmp(Decode_Type,'Multi_Dimensional_BIMLPCM'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            
        elseif(strcmp(Decode_Type,'MLPC'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
         
        elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
           
        elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            
        elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            
        elseif(strcmp(Decode_Type,'Compound_MLPC')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));            
            
        elseif(strcmp(Decode_Type,'Separated_BIPCM'))
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            
        elseif(strcmp(Decode_Type,'BIPCM'))
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
        
        elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
          
        elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            
        elseif(strcmp(Decode_Type,'Compound_BIPCM'))
        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
        
        elseif(strcmp(Decode_Type,'BIMLPCM')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
           
        elseif(strcmp(Decode_Type,'BIMLPCM_BIPCM_MLPC')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            
        elseif(strcmp(Decode_Type,'BIMLPCM_gRS4')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
        
        elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
        
        end
        
    end
        
    Max_Num_Of_Iterations = floor((Min_Num_Of_Errors_Sufficient_Statistics/BLER));
    Num_Of_Iterations_Print_State = 100;
    Min_Num_Of_Iterations = 99;

    Pavg = 1;
    N0 = Pavg/SNR;
    Sigma = sqrt(N0/2);

    if(strcmp(Decode_Type,'gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
        U_N = N/2;
        Is_GF4 = true;
    else %gArikan
        U_N = N;
        Is_GF4 = false;
    end
    
    Bit_Level_U_N = N/Num_Of_Encoded_Symbols_In_Constellation;
    
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
    Multilevel_Bit_Level_U_N = Constellation_2D_Dimension*Bit_Level_U_N;
            
    U = zeros(NumOfWorkers,N);
    X_after_encoding = zeros(NumOfWorkers,N);
    X_before_compound = zeros(NumOfWorkers,N);
    X_after_encoding_gf = zeros(NumOfWorkers,Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
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
    
%     for ParIter=1:NumOfWorkers
    parfor ParIter=1:NumOfWorkers
        
        for Iter=1:1:floor(Max_Num_Of_Iterations/NumOfWorkers)
            
            Iter_Vec(1,ParIter) = Iter_Vec(1,ParIter) + 1;

        %-----------------------------generate input-----------------------------%
        %Generate random input
        
            U(ParIter,:) = randi([0,1],1,N);
%             U(ParIter,:) = zeros(1,N);
            
        %---------------------------------encode---------------------------------%
            %Encode input into output
            if(strcmp(Decode_Type,'Separated_BIPCM'))
                
                X_after_encoding(ParIter,:) = gArikan_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation);
            
            elseif(strcmp(Decode_Type,'Compound_BIPCM'))
                
                [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
               
            elseif(strcmp(Decode_Type,'Multi_Dimensional_MLPC'))
                
               X_after_encoding(ParIter,:) = gGF2_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation);
               
            elseif(strcmp(Decode_Type,'Multi_Dimensional_BIMLPCM'))
               
               [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
                
            elseif(strcmp(Decode_Type,'MLPC'))
                
               X_after_encoding(ParIter,:) = gGF2_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation);
            
            elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                
               [X_after_encoding(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation/2);
               
            elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                
               [X_after_encoding(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation/2);
               
            elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
                
               [X_after_encoding(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = Mixed_Kernels_gGF4_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation/2);
               
            elseif(strcmp(Decode_Type,'Compound_MLPC'))
                
                [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
            
            elseif(strcmp(Decode_Type,'BIMLPCM'))
                
               [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
           
            elseif(strcmp(Decode_Type,'BIMLPCM_BIPCM_MLPC'))
                
               [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
               
            elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
                
               [X_after_encoding(ParIter,:),X_before_compound(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m/2);
            
            elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
                
               [X_after_encoding(ParIter,:),X_before_compound(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m/2);
                    
            else %BPSK,BIPCM,Multi_Dimensional_BIPCM
                
                %X_after_encoding(ParIter,:) = gArikan_Polar_Code_Encoder(U(ParIter,:));
                
                X_after_encoding(ParIter,:) = mod(U(ParIter,:)*g,2);

            end

        %--------------------constellation Channel and Decode--------------------%
            
            if(strcmp(Constellation_Type,'BPSK'))
                
                %Transform the output bits into output symbols
                X_after_constellation(ParIter,:) = BPSK(X_after_encoding(ParIter,:));
                
                %Get the channel output
                Y(ParIter,:) = X_after_constellation(ParIter,:) + Sigma.*randn(1,Bit_Level_U_N);
                
                %Decode
                if(strcmp(Decode_Type,'gArikan'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BPSK_SCL_Decoder(L,GA_CRC_Length,Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                    
                elseif(strcmp(Decode_Type,'gRS4'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BPSK_SCL_Decoder(L,GA_CRC_Length,Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                 
                elseif(strcmp(Decode_Type,'gArikan_GF4'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BPSK_SCL_Decoder(L,GA_CRC_Length,Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                    
                end
                
            else % PAM/PSK/QAM/APSK
                
                %Transform the output bits into output symbols
                X_after_constellation(ParIter,:) = Constellation(m,Constellation_Mapping_Array,X_after_encoding(ParIter,:));
               
                %Get the channel output
                if(strcmp(Constellation_Type,'PAM'))
                    Y(ParIter,:) = X_after_constellation(ParIter,:) + Sigma.*randn(1,Multilevel_Bit_Level_U_N);
                else % PSK/QAM/APSK/Multilevel-Constellations
                    Y(ParIter,:) = X_after_constellation(ParIter,:) + Sigma.*randn(1,Multilevel_Bit_Level_U_N)+1i*Sigma.*randn(1,Multilevel_Bit_Level_U_N);
                end
                
                %Decode
                
                if(strcmp(Decode_Type,'Multi_Dimensional_MLPC'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_after_encoding(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_MLPC_SCL_Decoder(L,GA_CRC_Length,m,m_SCL,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                    
                elseif(strcmp(Decode_Type,'Multi_Dimensional_BIMLPCM'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    
                elseif(strcmp(Decode_Type,'Multi_Dimensional_BIPCM'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Multi_Dimensional_BIPCM_SCL_Decoder(L,GA_CRC_Length,m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                    
                elseif(strcmp(Decode_Type,'MLPC'))
                                            
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_after_encoding(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    
                    if(L==1)
                        if(N==m)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end
                    else
                        if(Is_Genie_Aided)
                            if(N==m)
                                [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                            else
                                [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                            end
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_MLPC_SCL_Decoder(L,GA_CRC_Length,m,m_SCL,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                
                elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                                                                
                    if(L==1)
                        if(N==m)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end
                    else
                        
                    end
                  
                elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                                                                
                    if(L==1)
                        if(N==m)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end
                    else
                        
                    end
                    
                elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
                                                                
                    if(L==1)
                        if(N==m)
                            
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_gArikan_MLPC_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end
                    else
                        
                    end
                    
                elseif(strcmp(Decode_Type,'Compound_MLPC'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Compound_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    
                elseif(strcmp(Decode_Type,'Separated_BIPCM'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Separated_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                
                elseif(strcmp(Decode_Type,'BIPCM'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                        if(Is_Genie_Aided)
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        else
                            [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIPCM_SCL_Decoder(L,GA_CRC_Length,m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                        end                                     
                    end
                   
                elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                                                            
                    end
                   
                elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
                    
                    if(L==1)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    else
                                                            
                    end
                    
                elseif(strcmp(Decode_Type,'Compound_BIPCM'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Compound_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                
                elseif(strcmp(Decode_Type,'BIMLPCM'))
                    
                    if(N==m)
                        X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);    
                    else
                        X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    end
                    
                elseif(strcmp(Decode_Type,'BIMLPCM_BIPCM_MLPC'))
                    
                    if(N==m)
                        
                    else
                        X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIMLPCM_BIPCM_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    end
                    
                elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
                    
                    if(N==m)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);    
                    else
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    end
                    
                elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
                    
                    if(N==m)
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);    
                    else
                        [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    end
                    
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
                
                if(Is_GF4)
                    temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:);
                else %gArikan
                    temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:).*(-1).^U(ParIter,:);
                end
                
%                 temp_LLR_Vec(ParIter,:) = (~temp_Errors_Vec(ParIter,:)).*temp_LLR_Vec(ParIter,:); % ???? ?????? ?? ???????
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
                    
                    if(Is_GF4)
                        First_Symbol(1,ParIter) = floor((First_Bit(1,ParIter)+1)/(N/U_N));
                        First_Symbol_Vec(ParIter,:) = sum(reshape(First_Bit_Vec(ParIter,:),2,[]));
                    else %gArikan
                        First_Symbol(1,ParIter) = First_Bit(1,ParIter);
                        First_Symbol_Vec(ParIter,:) = First_Bit_Vec(ParIter,:);
                    end
                    
                    temp_LLR_Vec(ParIter,:) = [ones(1,First_Symbol(1,ParIter)),zeros(1,U_N-First_Symbol(1,ParIter))];
                    LLR_Divide_Vec(ParIter,:) = LLR_Divide_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
                    
                    if(Is_GF4)
                        temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:).*temp_LLR_Vec(ParIter,:);
                    else %gArikan
                        temp_LLR_Vec(ParIter,:) = (Estimated_L(ParIter,:).*(-1).^U(ParIter,:)).*temp_LLR_Vec(ParIter,:);
                    end
                    
                    LLR_Avg_Vec(ParIter,:) = LLR_Avg_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:);
                    LLR_Var_Vec(ParIter,:) =  LLR_Var_Vec(ParIter,:) + temp_LLR_Vec(ParIter,:).^2;
                    
                catch
                                        
                    LLR_Divide_Vec(ParIter,:) = LLR_Divide_Vec(ParIter,:) + 1;
                    
                    if(Is_GF4)
                        temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:);
                    else %gArikan
                        temp_LLR_Vec(ParIter,:) = Estimated_L(ParIter,:).*(-1).^U(ParIter,:);
                    end
                    
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
        
        if(Is_GF4)
            SER_Vec_temp = reshape(1 - SER_Vec,2,[]);
            SER_Vec = 1 - SER_Vec_temp(1,:).*SER_Vec_temp(2,:);
            SER_Vec = min(0.75,SER_Vec);
            [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec));
            Sorted_SER_Index_Vec = size(SER_Vec,2) - Sorted_SER_Index_Vec + 1;
            Is_Frozen_Bit_Index_Vec = zeros(1,U_N);
            for Index=Sorted_SER_Index_Vec(1,R*U_N+1+(GA_CRC_Length>0)*GA_CRC_Length/2:U_N)
                Is_Frozen_Bit_Index_Vec(1,Index) = 1;
            end
            Is_Frozen_Bit_Index_Vec(1,:) = Is_Frozen_Bit_Index_Vec(1,:) == ones(1,U_N);
            Is_Frozen_Bit_Index_Vec = reshape([Is_Frozen_Bit_Index_Vec;Is_Frozen_Bit_Index_Vec],1,[]);
        else %gArikan
            SER_Vec = min(0.5,SER_Vec);
            [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec));
            Sorted_SER_Index_Vec = size(SER_Vec,2) - Sorted_SER_Index_Vec + 1;
            Is_Frozen_Bit_Index_Vec = zeros(1,N);
            for Index=Sorted_SER_Index_Vec(1,R*N+1+(GA_CRC_Length>0)*GA_CRC_Length:N)
                Is_Frozen_Bit_Index_Vec(1,Index) = 1;
            end
            Is_Frozen_Bit_Index_Vec(1,:) = Is_Frozen_Bit_Index_Vec(1,:) == ones(1,N);
        end
        
        Total_Number_Of_Iter = sum(Iter_Vec);
                    
%         Estimated_LLR_Avg_Vec = sum(LLR_Avg_Vec,1)./(Total_Number_Of_Iter-Num_Of_Errors);
%         Estimated_LLR_Var_Vec = ((Total_Number_Of_Iter-Num_Of_Errors)/(Total_Number_Of_Iter-Num_Of_Errors-1))*(sum(LLR_Var_Vec,1)./(Total_Number_Of_Iter-Num_Of_Errors)-Estimated_LLR_Avg_Vec.^2);
        Estimated_LLR_Avg_Vec = sum(LLR_Avg_Vec,1)./(Total_Number_Of_Iter);
        Estimated_LLR_Var_Vec = ((Total_Number_Of_Iter)/(Total_Number_Of_Iter-1))*(sum(LLR_Var_Vec,1)./(Total_Number_Of_Iter)-Estimated_LLR_Avg_Vec.^2);
        temp_Index_LLR = Estimated_LLR_Avg_Vec>0;
        Estimated_SER_Vec = zeros(1,U_N);
        Estimated_SER_Vec(temp_Index_LLR) = 0.5.*(1+erf((-Estimated_LLR_Avg_Vec(temp_Index_LLR)./sqrt(Estimated_LLR_Var_Vec(temp_Index_LLR))./sqrt(2))));

        if(Is_GF4)
            Estimated_SER_Vec(Estimated_LLR_Avg_Vec<=0) = 0.75;
            [Sorted_Estimated_SER_Vec,Sorted_Estimated_SER_Index_Vec] = sort(wrev(Estimated_SER_Vec));
            Sorted_Estimated_SER_Index_Vec = size(Estimated_SER_Vec,2) - Sorted_Estimated_SER_Index_Vec + 1;
            Estimated_Is_Frozen_Bit_Index_Vec = zeros(1,U_N);
            for Index=Sorted_Estimated_SER_Index_Vec(1,R*U_N+1+(GA_CRC_Length>0)*GA_CRC_Length/2:U_N)
                Estimated_Is_Frozen_Bit_Index_Vec(1,Index) = 1;
            end
            Estimated_Is_Frozen_Bit_Index_Vec(1,:) = Estimated_Is_Frozen_Bit_Index_Vec(1,:) == ones(1,U_N);
            Estimated_Is_Frozen_Bit_Index_Vec = reshape([Estimated_Is_Frozen_Bit_Index_Vec;Estimated_Is_Frozen_Bit_Index_Vec],1,[]);
        else %gArikan
            Estimated_SER_Vec(Estimated_LLR_Avg_Vec<=0) = 0.5;
            [Sorted_Estimated_SER_Vec,Sorted_Estimated_SER_Index_Vec] = sort(wrev(Estimated_SER_Vec));
            Sorted_Estimated_SER_Index_Vec = size(Estimated_SER_Vec,2) - Sorted_Estimated_SER_Index_Vec + 1;
            Estimated_Is_Frozen_Bit_Index_Vec = zeros(1,U_N);
            for Index=Sorted_Estimated_SER_Index_Vec(1,R*U_N+1+(GA_CRC_Length>0)*GA_CRC_Length:U_N)
                Estimated_Is_Frozen_Bit_Index_Vec(1,Index) = 1;
            end
            Estimated_Is_Frozen_Bit_Index_Vec(1,:) = Estimated_Is_Frozen_Bit_Index_Vec(1,:) == ones(1,U_N);
            Estimated_Is_Frozen_Bit_Index_Vec = reshape([Estimated_Is_Frozen_Bit_Index_Vec;Estimated_Is_Frozen_Bit_Index_Vec],1,[]);
        end
        
        BLER_lower_bound = max(Sorted_SER_Vec(1:U_N*R));
        BLER_upper_bound = min(sum(Sorted_SER_Vec(1:U_N*R)),0.5);
        
        Simulated_BLER = zeros(0);
        Simulated_BER = zeros(0);
%         Simulated_First_Bit_Vec = zeros(0);
        Simulated_First_Symbol_Vec = zeros(0);
        Simulated_LLR_Avg_Vec = zeros(0);
        Simulated_LLR_Var_Vec = zeros(0);
        
    else
        
        if(NumOfWorkers==1) 
            Simulated_BLER = BLER_Errors_Vec./Iter_Vec; 
            Simulated_BER = sum(Errors_Vec)./(Iter_Vec*N*R);
%             Simulated_First_Bit_Vec = First_Bit_Vec./Iter_Vec;
            Simulated_First_Symbol_Vec = First_Symbol_Vec./Iter_Vec;
            Simulated_LLR_Avg_Vec = LLR_Avg_Vec./LLR_Divide_Vec;
            Simulated_LLR_Var_Vec = LLR_Var_Vec./LLR_Divide_Vec-Simulated_LLR_Avg_Vec.^2;
        else 
            BLER_Vec = zeros(1,NumOfWorkers);
            BER_Vec = zeros(1,NumOfWorkers);
            for i=1:1:NumOfWorkers     
                BLER_Vec(i) = BLER_Errors_Vec(i)./sum(Iter_Vec);
                BER_Vec(i) = sum(Errors_Vec(i,:))./(sum(Iter_Vec)*N*R);
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

