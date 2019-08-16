 function [LLR_PDF_Vec] = Simulate_Constellation_Mapping_PDF(NumOfWorkers,g,g0,Decode_Type,Constellation_Type,m,m_BIPCM,Constellation_Mapping_Array,Num_Of_Iterations,SNR,LLR_Bin_Border_Vec)
                                                                                                                                                                                                            
    %----------------------------input parameters----------------------------%
    
    N = m;
    
    M = 2^m;
    
    Num_Of_Encoded_Symbols_In_Constellation = m;
    
    Is_Genie_Aided = true;
    
    Is_Frozen_Bit_Index_Vec = zeros(1,N);
    
    if(M==2)
                
        future_matrix = zeros(0,0);
                
    elseif(mod(log(M)/log(4),1)==0)

        if(strcmp(Decode_Type,'MLPC'))
              
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            g = 1;
        
        elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            g = 1;
           
        elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                        
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation-1,1));
            g = 1;
            
        elseif(strcmp(Decode_Type,'Compound_MLPC')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            g = 1;
            
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
            g = 1;
        
        elseif(strcmp(Decode_Type,'BIMLPCM')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            g = 1;
        
        elseif(strcmp(Decode_Type,'BIMLPCM_gRS4')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            g = 1;
        
        elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4')) 
            
            future_matrix = transpose(Create_All_Options(Num_Of_Encoded_Symbols_In_Constellation,1));
            g = 1;
            
        end
        
    end
    
    Num_Of_Iterations_Print_State = 10000;
    Num_Of_Bins = size(LLR_Bin_Border_Vec,2);

    %-----------------------calc and check parameters-----------------------%

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
    else %gArikan
        U_N = N;
        Is_GF4 = false;
    end

    Bit_Level_U_N = N/Num_Of_Encoded_Symbols_In_Constellation;
            
    U = zeros(NumOfWorkers,N);
    X_after_encoding = zeros(NumOfWorkers,N);
    X_before_compound = zeros(NumOfWorkers,N);
    X_after_encoding_gf = zeros(NumOfWorkers,Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
    X_after_constellation = zeros(NumOfWorkers,Bit_Level_U_N);
    Y = zeros(NumOfWorkers,Bit_Level_U_N);
    Estimated_U = zeros(NumOfWorkers,N);
    Estimated_L = zeros(NumOfWorkers,U_N);
    Estimated_LLR = zeros(NumOfWorkers,U_N);
    Iter_Vec = zeros(1,NumOfWorkers);
    
    if(Is_GF4)
        Bin_Index = zeros(NumOfWorkers,m/2);
        LLR_Bin_Vec = zeros(NumOfWorkers,m/2,Num_Of_Bins);
        temp_LLR_Bin_Vec = zeros(NumOfWorkers,m/2,Num_Of_Bins);
        LLR_PDF_Vec = zeros(m/2,Num_Of_Bins);
    else
        Bin_Index = zeros(NumOfWorkers,m);
        LLR_Bin_Vec = zeros(NumOfWorkers,m,Num_Of_Bins);
        temp_LLR_Bin_Vec = zeros(NumOfWorkers,m,Num_Of_Bins);
        LLR_PDF_Vec = zeros(m,Num_Of_Bins);
    end
    
    temp_Bin_Index = 1:1:Num_Of_Bins;   
    
%     for ParIter=1:NumOfWorkers
    parfor ParIter=1:NumOfWorkers,
        
        for Iter=1:1:floor(Num_Of_Iterations/NumOfWorkers)
            
            Iter_Vec(1,ParIter) = Iter_Vec(1,ParIter) + 1;

        %-----------------------------generate input-----------------------------%
        %Generate random input
        
            U(ParIter,:) = randi([0,1],1,N);
            
        %---------------------------------encode---------------------------------%
            %Encode input into output
            if(strcmp(Decode_Type,'Separated_BIPCM'))
                
                X_after_encoding(ParIter,:) = gGF2_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation);
            
            elseif(strcmp(Decode_Type,'Compound_BIPCM'))
                
                [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
               
            elseif(strcmp(Decode_Type,'MLPC'))
                
               X_after_encoding(ParIter,:) = gGF2_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation);
            
            elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                
               [X_after_encoding(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation/2);
               
            elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                
               [X_after_encoding(ParIter,:),X_after_encoding_gf(ParIter,:,:)] = gGF4_Separated_Polar_Code_Encoder(g,U(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation/2);
               
            elseif(strcmp(Decode_Type,'Compound_MLPC'))
                
                [X_after_encoding(ParIter,:),X_before_compound(ParIter,:)] = gGF2_Compound_Polar_Code_Encoder(g,g0,U(ParIter,:),m);
            
            elseif(strcmp(Decode_Type,'BIMLPCM'))
                
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
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                
                elseif(strcmp(Decode_Type,'gRS4'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                
                elseif(strcmp(Decode_Type,'gArikan_GF4'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BPSK_SC_Decoder(Y(ParIter,:),U(ParIter,:),Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    
                end
                    
            else % PAM/PSK/QAM/APSK
                
                %Transform the output bits into output symbols
                X_after_constellation(ParIter,:) = Constellation(m,Constellation_Mapping_Array,X_after_encoding(ParIter,:));
               
                %Get the channel output
                if(strcmp(Constellation_Type,'PAM'))
                    Y(ParIter,:) = X_after_constellation(ParIter,:) + Sigma.*randn(1,Bit_Level_U_N);
                else % PSK/QAM/APSK
                    Y(ParIter,:) = X_after_constellation(ParIter,:) + Sigma.*randn(1,Bit_Level_U_N)+1i*Sigma.*randn(1,Bit_Level_U_N);
                end
                
                %Decode
                
                if(strcmp(Decode_Type,'MLPC'))
                                            
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_after_encoding(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    
                elseif(strcmp(Decode_Type,'MLPC_gRS4'))
                                                                
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    
                elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
                                                                
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    
                elseif(strcmp(Decode_Type,'Compound_MLPC'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Compound_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                    
                elseif(strcmp(Decode_Type,'Separated_BIPCM'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Separated_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                
                elseif(strcmp(Decode_Type,'BIPCM'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                
                elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                  
                elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
                    
                elseif(strcmp(Decode_Type,'Compound_BIPCM'))
                    
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_Compound_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding(ParIter,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                
                elseif(strcmp(Decode_Type,'BIMLPCM'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
                
                elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gRS4_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
    
                elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
                    
                    X_after_encoding_gf(ParIter,:,:) = reshape(X_before_compound(ParIter,:),Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);
                    [Estimated_U(ParIter,:),Estimated_L(ParIter,:)] = gArikan_GF4_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y(ParIter,:),U(ParIter,:),X_after_encoding_gf(ParIter,:,:),Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0);
    
                end
                
                
            end
            
            if(Is_GF4)
                Estimated_LLR(ParIter,:) = Estimated_L(ParIter,:);
            else
                Estimated_LLR(ParIter,:) = Estimated_L(ParIter,:).*(-1).^U(ParIter,:);
            end
            
            [~,Bin_Index(ParIter,:)] = min(abs(bsxfun(@minus,Estimated_LLR(ParIter,:),LLR_Bin_Border_Vec.')));
            
            temp_LLR_Bin_Vec(ParIter,:,:) = bsxfun(@eq,Bin_Index(ParIter,:).',temp_Bin_Index);
            
            LLR_Bin_Vec(ParIter,:,:) = LLR_Bin_Vec(ParIter,:,:) + temp_LLR_Bin_Vec(ParIter,:,:);

            if(mod(Iter_Vec(1,ParIter),Num_Of_Iterations_Print_State)==0)
                sprintf('worker %d finish iteration %d (out of %d iterations)', ParIter,Iter_Vec(1,ParIter),floor(Num_Of_Iterations/NumOfWorkers))
            end
            
        end
        
    end
    
        LLR_PDF_Vec(:,:) = sum(LLR_Bin_Vec,1)./Num_Of_Iterations;
        
end

