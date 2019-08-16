 function [Estimated_U,Estimated_L] = gArikan_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
     
    Estimated_U = zeros(1,m);
        
    Estimated_L = zeros(1,m);

    for i=1:1:m
    
        Estimated_L(1,i) = AWGN_MLPC_LLR_N_m(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X);
        
        Estimated_U(1,i) = floor((1-sign(Estimated_L(1,i)))/2);
        
    end
    
 end