 function [Estimated_U,Estimated_L] = gArikan_GF4_MLPC_SC_Decoder_N_m(m,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
     
    Estimated_U = zeros(1,m);
        
    Estimated_L = zeros(1,m/2);

    for i=1:2:m
    
        [LL0,LL1,LL2,LL3] = AWGN_GF4_MLPC_LL_N_m(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,X);
                
        lambda = bsxfun(@minus,LL0,[LL0;LL1;LL2;LL3]);
        
        [~,Estimated_U_gf4] = min(lambda);
        Estimated_U(1,i:i+1) = [(Estimated_U_gf4 - 1 - mod(Estimated_U_gf4 - 1,2))/2,mod(Estimated_U_gf4 - 1,2)];
        
        U_gf4 = reshape(U,2,[]);
        U_gf4 = 2*U_gf4(1,:) + U_gf4(2,:);
        Index = [0:1:size(Estimated_L,2)-1];
        right_Estimated_L = Estimated_L(U_gf4+1+4.*Index);
        Estimated_L(U_gf4+1+4.*Index) = -inf;
        sorted_Estimated_L = sort(lambda);
        Estimated_L(1,i:i+1) = - right_Estimated_L + sorted_Estimated_L(2,:) - log(sum(exp(bsxfun(@minus,sorted_Estimated_L(2,:),sorted_Estimated_L(2:4,:)))));

    end
    
 end