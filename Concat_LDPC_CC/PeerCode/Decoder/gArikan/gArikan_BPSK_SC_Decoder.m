 function [Estimated_U,Estimated_L] = gArikan_BPSK_SC_Decoder(Y,U,Sigma,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);

    lambda = AWGN_BPSK_LLR(Y,Sigma);

    [Estimated_U,~,Estimated_L] = gArikan_SC_Decoder(lambda,U,N,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided);
    
 end