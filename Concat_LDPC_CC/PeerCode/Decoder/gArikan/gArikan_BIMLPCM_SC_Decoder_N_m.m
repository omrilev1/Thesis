  function [Estimated_U,Estimated_L] = gArikan_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    N = length(Y);
    
    U_N = m*N;
        
    Estimated_U = zeros(1,U_N);
    
%     Estimated_X = zeros(1,m,N);
    Estimated_X = X; %***Same BLER but not same BER***
        
    Estimated_L = zeros(1,U_N);

    if(m_BIPCM>0)
        
        [Estimated_lambda] = AWGN_BIMLPCM_LLR_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,1,zeros(0),g0);
        
        [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
    end
    
    for i=m_BIPCM+1:1:m
        
        Estimated_L(1,i) = AWGN_BIMLPCM_LLR_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X,g0);
        
        Estimated_U(1,i) = floor((1-sign(Estimated_L(1,i)))/2);
        
    end

 end