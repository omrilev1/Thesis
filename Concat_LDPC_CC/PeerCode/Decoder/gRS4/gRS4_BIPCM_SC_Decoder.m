 function [Estimated_U,Estimated_L] = gRS4_BIPCM_SC_Decoder(m,Constellation_Mapping_Array,Y,U,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);
    
    U_N = m*N;
    
    [LL0,LL1,LL2,LL3] = AWGN_GF4_BIPCM_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);
    
    lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]); 

    [Estimated_U,~,Estimated_L] = gRS4_SC_Decoder(lambda,U,U_N,Is_Frozen_Bit_Index_Vec(1:2:end),Is_Genie_Aided);
    
    U_gf4 = reshape(U,2,[]);
    U_gf4 = 2*U_gf4(1,:) + U_gf4(2,:);
    Index = [0:1:size(Estimated_L,2)-1];
    right_Estimated_L = Estimated_L(U_gf4+1+4.*Index);
    Estimated_L(U_gf4+1+4.*Index) = -inf;
    sorted_Estimated_L = sort(Estimated_L);
    Estimated_L = - right_Estimated_L + sorted_Estimated_L(2,:) - log(sum(exp(bsxfun(@minus,sorted_Estimated_L(2,:),sorted_Estimated_L(2:4,:)))));

 end