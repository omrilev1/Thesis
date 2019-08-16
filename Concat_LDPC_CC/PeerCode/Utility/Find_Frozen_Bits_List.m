%function [ Frozen_Bits_List_Vec ] = Find_Frozen_Bits_List(Sorted_SER_Vec,Sorted_SER_Index_Vec,R,L_Vec,file_path,file_name)
function [ Frozen_Bits_List_Vec, diff ] = Find_Frozen_Bits_List(Sorted_SER_Vec,Sorted_SER_Index_Vec,R,L_Vec)

    Sorted_SER_Vec(Sorted_SER_Vec>0.5) = 0.5;

    Sorted_SER_Vec(Sorted_SER_Vec==0) = inf;
    
    %Min_SER = min(Sorted_SER_Vec./100);
    
    %Min_SER_log = log(Min_SER);
    
    Sorted_SER_Vec(Sorted_SER_Vec==inf) = min(Sorted_SER_Vec./100);
    
    Sorted_SER_Vec = log(Sorted_SER_Vec);

    N = length(Sorted_SER_Vec);
    
    Frozen_Bits_List_Vec = zeros(length(L_Vec),N);
    
    diff = zeros(1,length(L_Vec));
    
    Is_Frozen_Bit_Index_Vec = zeros(1,N);
    for Index=Sorted_SER_Index_Vec(1,R*N+1:N)
        Is_Frozen_Bit_Index_Vec(1,Index) = 1;
    end
    Is_Frozen_Bit_Index_Vec(1,:) = Is_Frozen_Bit_Index_Vec(1,:) == ones(1,N);

    temp = sortrows([Sorted_SER_Index_Vec;Sorted_SER_Vec].');

    Unsorted_SER_Vec = temp(:,2);

    for index=1:1:length(L_Vec)
        
        L = L_Vec(index);
    
        LLR_Sum_Length = log2(L) + 1;

        for i=1:1:R*N

            BER_Sum_Vec = zeros(1,N-LLR_Sum_Length);
            
            for j=1:1:N-LLR_Sum_Length

                if(Frozen_Bits_List_Vec(index,j))
                    BER_Sum_Vec(j) = -inf;
                    continue;
                end

                k=0;
                Sum_Size=0;
                while(Sum_Size<LLR_Sum_Length)

                    if(~Frozen_Bits_List_Vec(index,j+k))
                        BER_Sum_Vec(j) = BER_Sum_Vec(j) + Unsorted_SER_Vec(j+k);
                        Sum_Size = Sum_Size + 1;
                    end

                    k=k+1;

                end

            end

            [~,Sorted_Index] = sort(BER_Sum_Vec);

            temp_LLR_Sum_Length = 0;
            Max_Index = Sorted_Index(end);
            for j=Sorted_Index(end):1:N

                if(temp_LLR_Sum_Length==LLR_Sum_Length)
                    break;
                elseif(Frozen_Bits_List_Vec(index,j))
                    continue;
                else

                    if(Unsorted_SER_Vec(j)>Unsorted_SER_Vec(Max_Index))
                        Max_Index = j;
                    end

                    temp_LLR_Sum_Length = temp_LLR_Sum_Length + 1;
                end

            end

            Frozen_Bits_List_Vec(index,Max_Index) = 1;

        end

        Frozen_Bits_List_Vec(index,:) = Frozen_Bits_List_Vec(index,:) == ones(1,N);

        diff(index) = sum(Frozen_Bits_List_Vec(index,:) ~= Is_Frozen_Bit_Index_Vec);
        
    end

    %save(fullfile(file_path,file_name),'Sorted_SER_Vec','Sorted_SER_Index_Vec','Is_Frozen_Bit_Index_Vec','L_Vec','diff','Frozen_Bits_List_Vec');
    
end

