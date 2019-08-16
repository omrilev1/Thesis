  function [Estimated_U,Estimated_L] = gArikan_Compound_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    N = length(Y);
    
    U_N = m*N;
    
    half_U_N = U_N/2;
    
    half_m = m/2;
                
    [lambda1] = AWGN_Compound_BIPCM_Not_Interleaved_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,false,zeros(m,N),g0);
    lambda = (max(0,lambda1(1:2:end)+lambda1(2:2:end)) + log(1+exp(-abs(lambda1(1:2:end)+lambda1(2:2:end))))) - (max(lambda1(1:2:end),lambda1(2:2:end)) + log(1+exp(-abs(lambda1(1:2:end)-lambda1(2:2:end)))));

    [~,Estimated_X,~] = gArikan_SC_Decoder(lambda,U(1:half_U_N),half_U_N,Is_Frozen_Bit_Index_Vec(1:half_U_N),Is_Genie_Aided);
    
    if(Is_Genie_Aided)
        X = reshape(X(1:half_U_N),half_m,[]);
        [lambda2] = AWGN_Compound_BIPCM_Not_Interleaved_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,true,X,g0);
    else
        Estimated_X = reshape(Estimated_X,half_m,[]);
        [lambda2] = AWGN_Compound_BIPCM_Not_Interleaved_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,true,Estimated_X,g0);
    end
    
    lambda = [lambda1(1:2:end,:);lambda2(2:2:end,:)];
    lambda = reshape(lambda,1,[]);
        
    [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided); 
    
  end