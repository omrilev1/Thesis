 function [Estimated_U,Estimated_L] = gArikan_Multi_Dimensional_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y,U,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
 
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
 
    N = length(Y);
    
    U_N_div_m = N/Constellation_2D_Dimension;
    
    U_N = m*U_N_div_m;
 
    lambda = AWGN_Multi_Dimensional_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);

    [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
 end