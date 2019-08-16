 function [Estimated_U,Estimated_L] = gArikan_SCL_Decoder_TopNew(L,GA_CRC_Length,LL0,LL1,U,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(LL0);
        
    [Estimated_U_List,Estimated_X_List,Estimated_L0_List,Estimated_L1_List,~] = gArikan_SCL_Decoder(LL0,LL1,L,U,N,Is_Frozen_Bit_Index_Vec);
    
    is_correct = ismember(U,Estimated_U_List,'rows');
    
    if(and(GA_CRC_Length>=0,is_correct))
        
        Index = find(sum(bsxfun(@eq,Estimated_U_List,U),2)==N);
        
        Estimated_U = U;

        Estimated_L = Estimated_L0_List(Index,:)-Estimated_L1_List(Index,:);
        
    else

        p = (~Estimated_X_List)*(LL0.')+(Estimated_X_List)*(LL1.');

        [~,ML_Index] = max(p);
        
        Estimated_U = Estimated_U_List(ML_Index,:);

        Estimated_L = Estimated_L0_List(ML_Index,:)-Estimated_L1_List(ML_Index,:);
        
    end
    
 end