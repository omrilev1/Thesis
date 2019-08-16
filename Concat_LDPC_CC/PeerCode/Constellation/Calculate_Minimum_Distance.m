function [d_min_square,Es] = Calculate_Minimum_Distance(Constellation_Mapping_Array,m_BIPCM,m_GF)

    M = size(Constellation_Mapping_Array,2);
    
    m = log2(M);
    
    m_BIPCM = max(m_BIPCM,0);
    
    [d_min_square1,Es1] = Calculate_Minimum_Distance_Single_Stage(Constellation_Mapping_Array,m_GF);

    [d_min_square2,Es2] = Calculate_Minimum_Distance_Multi_Stage(Constellation_Mapping_Array,m_BIPCM,m_GF);
    
    d_min_square = zeros(1,m/m_GF);
    d_min_square(1:m_BIPCM/m_GF) = d_min_square1(1:m_BIPCM/m_GF);
    d_min_square(m_BIPCM/m_GF+1:m/m_GF) = d_min_square2(m_BIPCM/m_GF+1:m/m_GF);
    
    Es = zeros(1,m/m_GF);
    Es(1:m_BIPCM/m_GF) = Es1(1:m_BIPCM/m_GF);
    Es(m_BIPCM/m_GF+1:m/m_GF) = Es2(m_BIPCM/m_GF+1:m/m_GF);
end
    
function [d_min_square,Es] = Calculate_Minimum_Distance_Multi_Stage(Constellation_Mapping_Array,m_BIPCM,m_GF)
    
    M = size(Constellation_Mapping_Array,2);

    m = log2(M);

    d_min_square = inf(1,m/m_GF);
    Es = inf(1,m/m_GF);
    
    [d_min_square_temp,Es_temp] = Calculate_Minimum_Distance_Single_Stage(Constellation_Mapping_Array,m_GF);
    
    d_min_square(1) = d_min_square_temp(1);
    Es(1) = Es_temp(1);
    
    if(m==m_GF)
        return
    end
    
    for j=0:1:2^m_GF-1

        Constellation_Mapping_Array_temp = zeros(1,M/2^m_GF);
        Index = 1;
        
        for i=0:1:M-1

            i_string = dec2base(i, 2^m_GF, m/m_GF);
            j_string = dec2base(j, 2^m_GF, m/m_GF);
            if(i_string(1)==j_string(1));
                Constellation_Mapping_Array_temp(1,Index) = Constellation_Mapping_Array(1,i+1);
                Index = Index + 1;
            end
            
        end
        
        [d_min_square_temp,Es_temp] = Calculate_Minimum_Distance(Constellation_Mapping_Array_temp,m_BIPCM-m_GF,m_GF);
        
        for i=1:1:m/m_GF-1
            if(d_min_square_temp(i)<d_min_square(1+i))
                d_min_square(1+i) = d_min_square_temp(i);
            end
            if(Es_temp(i)<Es(1+i))
                Es(1+i) = Es_temp(i);
            end
        end
        
    end

end

function [d_min_square,Es] = Calculate_Minimum_Distance_Single_Stage(Constellation_Mapping_Array,m_GF)

    M = size(Constellation_Mapping_Array,2);

    m = log2(M);

    d_min_square = inf;

    for i=1:1:M
        for j=i+1:1:M
            d_min_square_temp = abs(Constellation_Mapping_Array(1,i) - Constellation_Mapping_Array(1,j))^2;
            if(d_min_square_temp<d_min_square)
                d_min_square = d_min_square_temp;
            end
        end
    end
    
    d_min_square = d_min_square.*ones(1,m/m_GF);
    
    Es = mean(abs(Constellation_Mapping_Array).^2)*ones(1,m/m_GF);

end
    
    

