function [ Mutual_Information ] = Constellation_Mutual_Information(Constellation_Mapping_Array,SNR_dB_EbN0,Decode_Type,R)
    
%     Constellation_Mutual_Information( Make_Constellation('QAM',2,'Gray',0),2,'BIPCM',0.5)
%     Constellation_Mutual_Information( Make_Constellation('APSK-8-8',4,'Gray',0),4,'BIPCM',0.5)
    close all;
  
    M = size(Constellation_Mapping_Array,2);
    m = log2(M);
    
    q = 1;
    k=m/q;
    
    N = M*10^1;
    
    SNR_dB_EsN0 = SNR_dB_EbN0 + 10*log10(m) + 10*log10(R);
    SNR = 10^(SNR_dB_EsN0/10);
    Pavg = 1;
    N0 = Pavg/(SNR);
    Sigma = sqrt(N0/2);

    Num_Of_Bins = (N/M)/1;
    Num_Of_STD = 4;
    
    Real_Bin_Border_Vec = zeros(M,Num_Of_Bins);
    Image_Bin_Border_Vec = zeros(M,Num_Of_Bins);
    for i=1:1:M
        Real_Bin_Border_Vec(i,:) = real(Constellation_Mapping_Array(i))-Num_Of_STD*Sigma:2*Num_Of_STD*Sigma/(Num_Of_Bins-1):real(Constellation_Mapping_Array(i))+Num_Of_STD*Sigma;
        Image_Bin_Border_Vec(i,:) = imag(Constellation_Mapping_Array(i))-Num_Of_STD*Sigma:2*Num_Of_STD*Sigma/(Num_Of_Bins-1):imag(Constellation_Mapping_Array(i))+Num_Of_STD*Sigma;
    end
    
    X = reshape(repmat(Constellation_Mapping_Array.',1,N/M),1,[]);
    Y = X + Sigma.*(randn(1,N)+1i*randn(1,N));
%     Y = reshape(Y,M,[]);
    
%     P_x = 1/M;

    temp_powers = (2.^(m-1:-1:0)).';
    temp_zeros = zeros(length(future_matrix(:,1)),1);
    exp_multiplier = (1/(2*Sigma^2));

    if(strcmp(Decode_Type,'BIPCM'))

        for i=0:1:m-1

            k=i+1;

            u = [future_matrix(:,1:k-1),temp_zeros,future_matrix(:,k:m-1)];
            t0 = u*temp_powers;
            t1 = t0+2^(m-k);

            X0 = Constellation_Mapping_Array(t0+1);
            X1 = Constellation_Mapping_Array(t1+1);

            p0(k,:) = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
            p1(k,:) = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);

        end
        
        Mutual_Information = 0;
        
        for i=1:1:m
            
            Mutual_Information = Mutual_Information + sum(p0(i).*log2((p0(i))./((sum(p0,1)+sum(p1,1))./(2^m))),1);
            Mutual_Information = Mutual_Information + sum(p1(i).*log2((p1(i))./((sum(p0,1)+sum(p1,1))./(2^m))),1);
            
        end
                
        Mutual_Information = Mutual_Information/2;
        
        
        
        
        
        
        
        
        
%         Real_Bin_Index = zeros(M,N/M);
%         Imag_Bin_Index = zeros(M,N/M);
%         Bin_Index = zeros(M,(N/M));
%         P_y_given_x = zeros(M,(N/M)^2);
%         
%         for i=1:1:M
%             [~,Real_Bin_Index(i,:)] = min(abs(bsxfun(@minus,real(Y(i,:)).',Real_Bin_Border_Vec(i,:))));
%             [~,Imag_Bin_Index(i,:)] = min(abs(bsxfun(@minus,imag(Y(i,:)).',Image_Bin_Border_Vec(i,:))));
%             Bin_Index(i,:) = Num_Of_Bins*(Imag_Bin_Index(i,:)-1)+Real_Bin_Index(i,:);
%             [Count,Index]=hist(Bin_Index(i,:),unique(Bin_Index(i,:)));
%             P_y_given_x(i,Index) = Count./(N/M);
%         end
%         
%         P_y = sum(P_y_given_x,1);
%         
%         Mutual_Information = P_x.*P_y_given_x*log2(P_y_given_x./repmat(P_y,M,1));
        
        
%         future_matrix = transpose(Create_All_Options(m-1,1));
%         
%         [ll0,ll1] = AWGN_BIPCM_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix);
% %         ll0 = reshape(ll0,m,[]);
% %         ll1 = reshape(ll1,m,[]);
%         p0 = exp(ll0);
%         p1 = exp(ll1);
%         p0 = p0./(sum(p0+p1,1)/m);
% %         p0 = bsxfun(@times,p0,1./(sum(p0+p1,1)/m));
%         
%         [~,Bin_Index] = min(abs(bsxfun(@minus,p0,Bin_Border_Vec.')));
%         
%         p0_Bin_Probability_Vec = squeeze(sum(reshape(bsxfun(@eq,Bin_Index.',temp_Bin_Index),m,[],Num_Of_Bins),2))./N;
%         
%         Mutual_Information = -sum(p0_Bin_Probability_Vec.*log2(p0_Bin_Probability_Vec),2);
%         
% %         p1 = 1-p0;
% 
% %         C = m + sum(sum(p0.*log2(p0)+p1.*log2(p1),2)./N,1);
        
    elseif(strcmp(Decode_Type,'MLPC'))
        
        Real_Bin_Index = zeros(M,N/M);
        Imag_Bin_Index = zeros(M,N/M);
        Bin_Index = zeros(M,(N/M));
        P_y_given_x = zeros(M,(N/M)^2);
        
        for i=1:1:M
            [~,Real_Bin_Index(i,:)] = min(abs(bsxfun(@minus,real(Y(i,:)).',Real_Bin_Border_Vec(i,:))));
            [~,Imag_Bin_Index(i,:)] = min(abs(bsxfun(@minus,imag(Y(i,:)).',Image_Bin_Border_Vec(i,:))));
            Bin_Index(i,:) = Num_Of_Bins*(Imag_Bin_Index(i,:)-1)+Real_Bin_Index(i,:);
            [Count,Index]=hist(Bin_Index(i,:),unique(Bin_Index(i,:)));
            P_y_given_x(i,Index) = Count./(N/M);
        end
        
        P_y = sum(P_y_given_x,1);
        
        Mutual_Information = P_x.*P_y_given_x*log2(P_y_given_x./repmat(P_y,M,1));
        
        
    
%         future_matrix = transpose(Create_All_Options(m-1,1));
%         
%         C = zeros(1,m);
%         
%         for i=1:1:m
%             
%             [ll0,ll1] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,i,U);
%             p0 = exp(ll0);
%             p1 = exp(ll1);
%             p0 = p0./(p0+p1);
%             p1 = p1./(p0+p1);
% 
%             C(i) = sum(p0.*log2(p0)+p1.*log2(p1),2)./N;
%             
%         end
%         
%         C = m + sum(C,2);
    
    end
    
%     C_AWGN = log2(1+SNR);
    
end

