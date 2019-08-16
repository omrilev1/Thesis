  function [Estimated_U,Estimated_L] = gArikan_GF4_BIMLPCM_SC_Decoder(m,m_BIPCM,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    N = length(Y);
    
    N_gf2 = 2*N;
    
    U_N = m*N;
        
    Estimated_U = zeros(1,U_N);
    
    Estimated_X = X; %***Same BLER but not same BER***
        
    Estimated_L = zeros(4,U_N/2);

    if(m_BIPCM>0)
        
        [LL0,LL1,LL2,LL3] = AWGN_GF4_BIMLPCM_LL(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,1,zeros(0),g0);

        Estimated_lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]); 
        
        [Estimated_U,~,Estimated_L] = gArikan_GF4_SC_Decoder(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec(1:2:end),Is_Genie_Aided);

    end
    
    for i=m_BIPCM+1:2:m
        
        if(Is_Genie_Aided)
            [LL0,LL1,LL2,LL3] = AWGN_GF4_BIMLPCM_LL(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X,g0);
        else
            [LL0,LL1,LL2,LL3] = AWGN_GF4_BIMLPCM_LL(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X,g0);
        end

        Estimated_lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]);
        
        index = (i+1)/2;
        
        [Estimated_U((index-1)*N_gf2+1:index*N_gf2),Estimated_X_temp,Estimated_L(:,(index-1)*N+1:index*N)] = gArikan_GF4_SC_Decoder(Estimated_lambda,U((index-1)*N_gf2+1:index*N_gf2),N_gf2,Is_Frozen_Bit_Index_Vec((index-1)*N_gf2+1:2:index*N_gf2),Is_Genie_Aided); 
        
        Estimated_X(1,i:i+1,:) = reshape(Estimated_X_temp,2,[]);
        
    end
    
    U_gf4 = reshape(U,2,[]);
    U_gf4 = 2*U_gf4(1,:) + U_gf4(2,:);
    Index = [0:1:size(Estimated_L,2)-1];
    right_Estimated_L = Estimated_L(U_gf4+1+4.*Index);
    Estimated_L(U_gf4+1+4.*Index) = -inf;
    sorted_Estimated_L = sort(Estimated_L);
    Estimated_L = - right_Estimated_L + sorted_Estimated_L(2,:) - log(sum(exp(bsxfun(@minus,sorted_Estimated_L(2,:),sorted_Estimated_L(2:4,:)))));

 end