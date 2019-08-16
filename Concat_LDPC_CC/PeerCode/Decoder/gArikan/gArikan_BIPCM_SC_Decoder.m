 function [Estimated_U,Estimated_L] = gArikan_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y,U,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);
    
    U_N = m*N;

    lambda = AWGN_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);

    [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
 end