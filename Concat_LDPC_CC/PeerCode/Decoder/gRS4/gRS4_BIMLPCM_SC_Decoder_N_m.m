  function [Estimated_U,Estimated_L] = gRS4_BIMLPCM_SC_Decoder_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided,g0)
  
    N = length(Y);
    
    U_N = m*N;
        
    Estimated_U = zeros(1,U_N);
    
%     Estimated_X = zeros(1,m,N);
    Estimated_X = X; %***Same BLER but not same BER***
        
    Estimated_L = zeros(1,U_N);

    if(m_BIPCM>0)
        
        [LL0,LL1,LL2,LL3] = AWGN_GF4_BIMLPCM_LL_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,1,zeros(0),g0);
                
        Estimated_lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]);
        
        [Estimated_U,~,Estimated_L] = gRS4_SC_Decoder(Estimated_lambda,U,U_N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
    end
    
    for i=m_BIPCM+1:2:m
        
        [LL0,LL1,LL2,LL3] = AWGN_GF4_BIMLPCM_LL_N_m(m,m_BIPCM,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X,g0);
                
        lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]);
        
        [~,Estimated_U_gf4] = min(lambda);
        Estimated_L(:,(i+1)/2) = lambda;
        Estimated_U(1,i:i+1) = [(Estimated_U_gf4 - 1 - mod(Estimated_U_gf4 - 1,2))/2,mod(Estimated_U_gf4 - 1,2)];
        
    end
    
    U_gf4 = reshape(U,2,[]);
    U_gf4 = 2*U_gf4(1,:) + U_gf4(2,:);
    Index = [0:1:size(Estimated_L,2)-1];
    right_Estimated_L = Estimated_L(U_gf4+1+4.*Index);
    Estimated_L(U_gf4+1+4.*Index) = -inf;
    sorted_Estimated_L = sort(Estimated_L);
    Estimated_L = - right_Estimated_L + sorted_Estimated_L(2,:) - log(sum(exp(bsxfun(@minus,sorted_Estimated_L(2,:),sorted_Estimated_L(2:4,:)))));
        

 end