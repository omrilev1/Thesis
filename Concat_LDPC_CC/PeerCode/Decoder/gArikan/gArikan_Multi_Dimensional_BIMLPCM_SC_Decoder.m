 function [Estimated_U,Estimated_L] = gArikan_Multi_Dimensional_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
 
    N = length(Y);
    
    U_N_div_m = N/Constellation_2D_Dimension;
    
    U_N = m*U_N_div_m;

    Estimated_U = zeros(1,U_N);
    
%     Estimated_X = zeros(1,m,U_N_div_m);
    Estimated_X = X; %***Same BLER but not same BER***
    
    Estimated_L = zeros(1,U_N);
    
    if(m_BIPCM>0)
        [Estimated_lambda] = AWGN_Multi_Dimensional_BIMLPCM_LLR(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,1,zeros(0),g0);
        [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
%         if(~Is_Genie_Aided)
%             [~,Estimated_X_temp,~] = SC_Decoder_Func(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
%             Estimated_X(1,1:m_BIPCM,:) = reshape(Estimated_X_temp(1:m_BIPCM*N),[],m_BIPCM).';
%         end
    end
    
    for i=m_BIPCM+1:1:m
        
        if(Is_Genie_Aided)
            [Estimated_lambda] = AWGN_Multi_Dimensional_BIMLPCM_LLR(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X,g0);
        else
            [Estimated_lambda] = AWGN_Multi_Dimensional_BIMLPCM_LLR(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X,g0);
        end
                
        [Estimated_U((i-1)*U_N_div_m+1:i*U_N_div_m),Estimated_X(1,i,:),Estimated_L((i-1)*U_N_div_m+1:i*U_N_div_m)] = gArikan_SC_Decoder(Estimated_lambda,U((i-1)*U_N_div_m+1:i*U_N_div_m),U_N_div_m,Is_Frozen_Bit_Index_Vec((i-1)*U_N_div_m+1:i*U_N_div_m),Is_Genie_Aided); 
    
    end
    
 end