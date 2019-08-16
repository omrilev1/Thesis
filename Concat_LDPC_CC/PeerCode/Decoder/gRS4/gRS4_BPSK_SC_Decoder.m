 function [Estimated_U,Estimated_L] = gRS4_BPSK_SC_Decoder(Y,U,Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);

    [LL0,LL1] = AWGN_BPSK_LL(Y,Sigma);
    
    LL0 = reshape(LL0,2,[]);
    LL1 = reshape(LL1,2,[]);

    LL00 = LL0(1,:) + LL0(2,:);
    LL01 = LL0(1,:) + LL1(2,:);
    LL10 = LL1(1,:) + LL0(2,:);
    LL11 = LL1(1,:) + LL1(2,:);
    
    lambda = bsxfun(@minus,LL00,[LL00;LL01;LL10;LL11]);   
    
    [Estimated_U,~,Estimated_L] = gRS4_SC_Decoder(lambda,U,N,Is_Frozen_Bit_Index_Vec(1:2:end),Is_Genie_Aided);
    
    U_gf4 = reshape(U,2,[]);
    U_gf4 = 2*U_gf4(1,:) + U_gf4(2,:);
    Index = [0:1:size(Estimated_L,2)-1];
    right_Estimated_L = Estimated_L(U_gf4+1+4.*Index);
    Estimated_L(U_gf4+1+4.*Index) = -inf;
    sorted_Estimated_L = sort(Estimated_L);
    Estimated_L = - right_Estimated_L + sorted_Estimated_L(2,:) - log(sum(exp(bsxfun(@minus,sorted_Estimated_L(2,:),sorted_Estimated_L(2:4,:)))));

 end