 function [Estimated_U,Estimated_L] = gArikan_Separated_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y,U,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);
        
    Estimated_U = zeros(1,N*m);
    
    Estimated_L = zeros(1,N*m);
    
    lambda = AWGN_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);
        
    for i=1:1:m
            
        [Estimated_U((i-1)*N+1:i*N),~,Estimated_L((i-1)*N+1:i*N)] = gArikan_SC_Decoder(lambda(i:m:end),U((i-1)*N+1:i*N),N,Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N),Is_Genie_Aided);
    
    end
    
 end