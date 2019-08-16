function [Estimated_U,Estimated_L] = gArikan_BIMLPCM_BIPCM_MLPC_SC_Decoder(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    N = length(Y);
    
    U_N = m*N;
    
%     Estimated_X = zeros(1,m,N);
    Estimated_X = X; %***Same BLER but not same BER***
        
    [Estimated_lambda] = AWGN_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);
    [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
    if(sum(Estimated_U~=U)~=0)
    
        for i=1:1:m

            if(Is_Genie_Aided)
                [Estimated_lambda] = AWGN_Compound_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X,g0);
            else
                [Estimated_lambda] = AWGN_Compound_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X,g0);
            end

            [Estimated_U((i-1)*N+1:i*N),Estimated_X(1,i,:),Estimated_L((i-1)*N+1:i*N)] = gArikan_SC_Decoder(Estimated_lambda,U((i-1)*N+1:i*N),N,Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N),Is_Genie_Aided); 

        end
        
    end

 end