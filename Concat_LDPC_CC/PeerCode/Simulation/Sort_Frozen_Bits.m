function [Sorted_By_Union_Bound,BLER] = Sort_Frozen_Bits(NumOfWorkers,N,R,BLER_Max,SNR,L,GA_CRC_Length,Min_Num_Of_Errors_Sufficient_Statistics,file_path)

    g_Arikan = [1,0;1,1];
    
    Num_Of_Iterations_Print_State = 100;
       
    Estimted_BLER_Sorted_By_Union_Bound = 10^0;
        
    g = g_Arikan;
    for i=1:1:log2(N)-1
        g = kron(g,g_Arikan);
    end
    g = bitrevorder(g);

    g0 = 1;

    Constellation_Mapping_Array = [-1,1];

    file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan',N,strrep(sprintf('%.3f',R),'.','_'));
    
    [~,temp_SER_Vec,~,~,~,~,~,~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,'gArikan',1,0,'BPSK',1,1,0,Constellation_Mapping_Array,N,R,Estimted_BLER_Sorted_By_Union_Bound,SNR(end),Min_Num_Of_Errors_Sufficient_Statistics,zeros(1,N));
    
    Optimal_Frozen_Bits_Search = 8;
    Optimal_Frozen_Bits_Search_Possibilities = 16;
    
    Optimal_Frozen_Bits_Search_Percentage = Optimal_Frozen_Bits_Search/N;
    Optimal_Frozen_Bits_Search_Percentage_Possibilities = Optimal_Frozen_Bits_Search_Possibilities/N;
    
    [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(temp_SER_Vec));
    Sorted_SER_Index_Vec = length(temp_SER_Vec) - Sorted_SER_Index_Vec + 1;
    
    start_R = floor((R+Optimal_Frozen_Bits_Search_Percentage)*N)/N;
    possibilities_R = floor((R-Optimal_Frozen_Bits_Search_Percentage_Possibilities)*N)/N;
                
    Is_Optimal_Frozen_Bit_Index_Vec = zeros(1,N);
    for Index=Sorted_SER_Index_Vec(1,start_R*N+1+(GA_CRC_Length>0)*GA_CRC_Length:N)
        Is_Optimal_Frozen_Bit_Index_Vec(1,Index) = 1;
    end
    Is_Optimal_Frozen_Bit_Index_Vec(1,:) = Is_Optimal_Frozen_Bit_Index_Vec(1,:) == ones(1,N);
        
    Frozen_Bits_Possibilities = Sorted_SER_Index_Vec(1,possibilities_R*N+1+(GA_CRC_Length>0)*GA_CRC_Length:R*N+(GA_CRC_Length>0)*GA_CRC_Length);
    Frozen_Bits_Possibilities = nchoosek(Frozen_Bits_Possibilities,floor(Optimal_Frozen_Bits_Search_Percentage*N));
    Num_Of_Frozen_Bits_Possibilities = size(Frozen_Bits_Possibilities,1);
    
    Is_Frozen_Bit_Index_Vecs = zeros(Num_Of_Frozen_Bits_Possibilities,N);
    LLR_Avg_Vec = zeros(Num_Of_Frozen_Bits_Possibilities,N);
    LLR_Var_Vec = zeros(Num_Of_Frozen_Bits_Possibilities,N);
    Estimated_SER_Vec = zeros(Num_Of_Frozen_Bits_Possibilities,N);  
    
    for i=1:1:Num_Of_Frozen_Bits_Possibilities
        
        Is_Frozen_Bit_Index_Vecs(i,:) = Is_Optimal_Frozen_Bit_Index_Vec;
        Is_Frozen_Bit_Index_Vecs(i,Frozen_Bits_Possibilities(i,:)) = true;
        
        [~,~,~,~,~,~,~,~,~,~,LLR_Avg_Vec(i,:),LLR_Var_Vec(i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,'gArikan',L,GA_CRC_Length,'BPSK',1,1,0,Constellation_Mapping_Array,N,R,Estimted_BLER_Sorted_By_Union_Bound,SNR(end),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vecs(i,:));
        
        temp_Index_LLR = LLR_Avg_Vec(i,:)>0;
        Estimated_SER_Vec(i,temp_Index_LLR) = 0.5.*(1+erf((-LLR_Avg_Vec(i,temp_Index_LLR)./sqrt(LLR_Var_Vec(i,temp_Index_LLR))./sqrt(2))));
        
        if(mod(i,Num_Of_Iterations_Print_State)==0)
            try
                save(fullfile(file_path,file_name),'Is_Frozen_Bit_Index_Vecs','Estimated_SER_Vec');
            end
            sprintf('finish frozen bits set %d (out of %d frozen bits sets)',i,Num_Of_Frozen_Bits_Possibilities)
        end
        
    end
    
    Union_Bound_Vec = sum(((~Is_Frozen_Bit_Index_Vecs).*Estimated_SER_Vec).');
        
    [~,Sorted_Union_Bound_Index_Vec] = sort(Union_Bound_Vec);
    Sorted_By_Union_Bound = Is_Frozen_Bit_Index_Vecs(Sorted_Union_Bound_Index_Vec,:);

    try
        save(fullfile(file_path,file_name),'Is_Frozen_Bit_Index_Vecs','Estimated_SER_Vec','Union_Bound_Vec','Sorted_By_Union_Bound');
    end
    
    BLER = inf(1,size(SNR,2));
    
    for i = 1:1:size(SNR,2)
        [~,~,~,~,~,~,~,BLER(i),~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,'gArikan',L,GA_CRC_Length,'BPSK',1,1,0,Constellation_Mapping_Array,N,R,BLER_Max,SNR(i),Min_Num_Of_Errors_Sufficient_Statistics,Sorted_By_Union_Bound(1,:));
        try
            save(fullfile(file_path,file_name),'Is_Frozen_Bit_Index_Vecs','Estimated_SER_Vec','Union_Bound_Vec','Sorted_By_Union_Bound','BLER');
        end
    end
    
end