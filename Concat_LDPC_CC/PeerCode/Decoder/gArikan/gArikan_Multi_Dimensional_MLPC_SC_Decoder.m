 function [Estimated_U,Estimated_L] = gArikan_Multi_Dimensional_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
 
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
 
    N = length(Y);
    
    U_N_div_m = N/Constellation_2D_Dimension;
    
    U_N = m*U_N_div_m;

    Estimated_U = zeros(1,U_N);
    
    Estimated_X = zeros(1,m,U_N_div_m);
    
    Estimated_L = zeros(1,U_N);
    
    for i=1:1:m
    
        if(Is_Genie_Aided)
            lambda = AWGN_Multi_Dimensional_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X);
        else
            lambda = AWGN_Multi_Dimensional_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X);
        end
        
        [Estimated_U((i-1)*U_N_div_m+1:i*U_N_div_m),Estimated_X(1,i,:),Estimated_L((i-1)*U_N_div_m+1:i*U_N_div_m)] = gArikan_SC_Decoder(lambda,U((i-1)*U_N_div_m+1:i*U_N_div_m),U_N_div_m,Is_Frozen_Bit_Index_Vec((i-1)*U_N_div_m+1:i*U_N_div_m),Is_Genie_Aided);
    
    end
    
 end