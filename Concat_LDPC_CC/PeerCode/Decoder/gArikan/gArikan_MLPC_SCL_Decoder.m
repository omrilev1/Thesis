function [Estimated_U,Estimated_L] = gArikan_MLPC_SCL_Decoder(L,GA_CRC_Length,m,m_SCL,Constellation_Mapping_Array,Y,U,X,Sigma,future_matrix,Is_Frozen_Bit_Index_Vec,Is_Genie_Aided)
    
    N = length(Y);
    
    N_SCL = m_SCL*N;
    
    N_gf2 = N*m;
    
    Estimated_U = zeros(1,N_gf2);
        
    Estimated_X = zeros(1,m,N);
    
    Estimated_L = zeros(1,N_gf2);
    
    Estimated_U_List = zeros(L,N_gf2);
        
    Estimated_X_List = zeros(L,m,N);
    
    Estimated_L0_List = zeros(L,N_gf2);
    
    Estimated_L1_List = zeros(L,N_gf2);
    
    
    temp_powers = (2.^(m-1:-1:0)).';
    
    Estimated_X_temp = zeros(1,m,N);
    
    Estimated_X_Index_temp = zeros(L,N);
    
    Estimated_X_Index_Add_temp = 1:2^(m-1):(2^(m-1))*N;
    
    Estimated_X_Matrix_temp = zeros(m,N);
    
    LL0_temp = zeros(2^(m-1),N);
    
    LL1_temp = zeros(2^(m-1),N);

    
    [LL0,LL1] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,1,Estimated_X_List);
    
    
    for i=1:1:m_SCL-1
    
        next_L = min(L,(2^sum(~Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N)))*size(LL0,1));
    
        [Estimated_U_List(1:next_L,(i-1)*N+1:i*N),Estimated_X_List(1:next_L,i,:),Estimated_L0_List(1:next_L,(i-1)*N+1:i*N),Estimated_L1_List(1:next_L,(i-1)*N+1:i*N),Indexes] = gArikan_SCL_Decoder(LL0,LL1,L,U((i-1)*N+1:i*N),N,Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N));

        Estimated_U_List(1:next_L,1:(i-1)*N) = Estimated_U_List(Indexes,1:(i-1)*N);
        Estimated_X_List(1:next_L,1:i-1,:) = Estimated_X_List(Indexes,1:i-1,:);
        Estimated_L0_List(1:next_L,1:(i-1)*N) = Estimated_L0_List(Indexes,1:(i-1)*N);
        Estimated_L1_List(1:next_L,1:(i-1)*N) = Estimated_L1_List(Indexes,1:(i-1)*N);
        
        if(2^i>L) %לממש אלגוריתם שבודק איזה X היו ולא סתם להריץ הכל
            
            LL0 = zeros(L,N);
            LL1 = zeros(L,N);
            
            for j=1:1:L

                [LL0(j,:),LL1(j,:)] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i+1,Estimated_X_List(j,:,:));

            end
            
        else
        
            past_matrix = future_matrix(m-i:m-1,1:2^i);

            for j=1:1:size(past_matrix,2)

                Estimated_X_temp(1,1:i,:) = repmat(past_matrix(:,j),1,N);

                [LL0_temp(j,:),LL1_temp(j,:)] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i+1,Estimated_X_temp);

            end

            for j=1:1:next_L

                Estimated_X_Matrix_temp(1:i,:) = Estimated_X_List(j,1:i,:);

                Estimated_X_Index_temp(j,:) = Estimated_X_Matrix_temp(1:i,:).'*temp_powers(end-i+1:end);

            end

            Estimated_X_Index_temp = bsxfun(@plus,Estimated_X_Index_temp(1:next_L,:),Estimated_X_Index_Add_temp);

            LL0 = LL0_temp(Estimated_X_Index_temp);
            LL1 = LL1_temp(Estimated_X_Index_temp);
            
        end
        
    end
    
    next_L = min(L,(2^sum(~Is_Frozen_Bit_Index_Vec((m_SCL-1)*N+1:m_SCL*N)))*size(LL0,1));
    
    [Estimated_U_List(1:next_L,(m_SCL-1)*N+1:m_SCL*N),Estimated_X_List(1:next_L,(m_SCL-1)*N+1:m_SCL*N),Estimated_L0_List(1:next_L,(m_SCL-1)*N+1:m_SCL*N),Estimated_L1_List(1:next_L,(m_SCL-1)*N+1:m_SCL*N),Indexes] = gArikan_SCL_Decoder(LL0,LL1,L,U((m_SCL-1)*N+1:m_SCL*N),N,Is_Frozen_Bit_Index_Vec((m_SCL-1)*N+1:m_SCL*N));
    
    Estimated_U_List(1:next_L,1:(m_SCL-1)*N) = Estimated_U_List(Indexes,1:(m_SCL-1)*N);
    Estimated_X_List(1:next_L,1:m_SCL-1,:) = Estimated_X_List(Indexes,1:m_SCL-1,:);
    Estimated_L0_List(1:next_L,1:(m_SCL-1)*N) = Estimated_L0_List(Indexes,1:(m_SCL-1)*N);
    Estimated_L1_List(1:next_L,1:(m_SCL-1)*N) = Estimated_L1_List(Indexes,1:(m_SCL-1)*N);
    
    is_correct = ismember(U(1:N_SCL),Estimated_U_List(:,1:N_SCL),'rows');
    
    if(and(GA_CRC_Length>=0,is_correct))
        
        Index = find(sum(bsxfun(@eq,Estimated_U_List(:,1:N_SCL),U(1:N_SCL)),2)==N_SCL);
        
        Estimated_U(1:N_SCL) = U(1:N_SCL);

        Estimated_X(1,1:m_SCL,:) = X(1,1:m_SCL,:);

        Estimated_L(1:N_SCL) = Estimated_L0_List(Index,1:N_SCL)-Estimated_L1_List(Index,1:N_SCL);
            
    else
        
        if(m_SCL==m) %חישוב מה המסלול הכי סביר
            
            temp_real = real(Y.');
            temp_imag = imag(Y.');
            
            past = zeros(m_SCL,N);
            
            p = zeros(1,L);
            
            for i=1:1:L
                
                past(:,:) = Estimated_X_List(i,:,:);
                
                t = past.'*temp_powers;

                X_Constellation = Constellation_Mapping_Array(t+1).';
                
                p(1,i) = sum(log(sum(exp(-((bsxfun(@minus,temp_real,real(X_Constellation))).^2+(bsxfun(@minus,temp_imag,imag(X_Constellation))).^2)*(1/(2*Sigma^2))),2)).');
            
            end
            
            [~,ML_Index] = max(p);
            
        else
            
            past_matrix = future_matrix(m-m_SCL:m-1,1:2^m_SCL);

            for i=1:1:size(past_matrix,2)

                Estimated_X_temp(1,1:m_SCL,:) = repmat(past_matrix(:,i),1,N);

                [LL0_temp(i,:),LL1_temp(i,:)] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,m_SCL+1,Estimated_X_temp);

            end

            for i=1:1:next_L

                Estimated_X_Matrix_temp(1:m_SCL,:) = Estimated_X_List(i,1:m_SCL,:);

                Estimated_X_Index_temp(i,:) = Estimated_X_Matrix_temp(1:m_SCL,:).'*temp_powers(end-m_SCL+1:end);

            end

            Estimated_X_Index_temp = bsxfun(@plus,Estimated_X_Index_temp(1:next_L,:),Estimated_X_Index_Add_temp);

            LL0 = LL0_temp(Estimated_X_Index_temp);
            LL1 = LL1_temp(Estimated_X_Index_temp);  

            p = prod((exp(LL0)+exp(LL1)).');

            [~,ML_Index] = max(p);
            
        end

        Estimated_U(1:N_SCL) = Estimated_U_List(ML_Index,1:N_SCL);

        Estimated_X(1,1:m_SCL,:) = Estimated_X_List(ML_Index,1:m_SCL,:);

        Estimated_L(1:N_SCL) = Estimated_L0_List(ML_Index,1:N_SCL)-Estimated_L1_List(ML_Index,1:N_SCL);

    end
    
    for i=m_SCL+1:1:m
    
        lambda = AWGN_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,Estimated_X);
        
        [Estimated_U((i-1)*N+1:i*N),Estimated_X(1,i,:),Estimated_L((i-1)*N+1:i*N)] = gArikan_SC_Decoder(lambda,U((i-1)*N+1:i*N),N,Is_Frozen_Bit_Index_Vec((i-1)*N+1:i*N),Is_Genie_Aided);
                                                                                                            
    end
    
end