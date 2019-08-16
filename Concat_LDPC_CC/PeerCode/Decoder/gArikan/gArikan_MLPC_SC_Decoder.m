 function [Estimated_U,Estimated_L] = gArikan_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
 
    N = length(Y);
    
    Estimated_U = zeros(1,N*m);
    
    Estimated_X = zeros(1,m,N);
    
    Estimated_L = zeros(1,N*m);

    for i=1:1:m
    
        if(Is_Genie_Aided)
            lambda = AWGN_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X);
        else
            lambda = AWGN_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X);
        end
        
        [Estimated_U((i-1)*N+1:i*N),Estimated_X(1,i,:),Estimated_L((i-1)*N+1:i*N)] = gArikan_SC_Decoder(lambda,U((i-1)*N+1:i*N),N,Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N),Is_Genie_Aided);
    
    end
    
 end