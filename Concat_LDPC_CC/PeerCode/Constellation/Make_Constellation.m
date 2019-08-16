function [Constellation_Mapping_Array] = Make_Constellation(Constellation_Type,m,Lable_Type,g0)
    
    global Max_Num_Of_Labels;
    
    Max_Num_Of_Labels = inf;

    M = 2^m;

    Lable_Array = zeros(1,M);
    
    Constellation_Mapping_Array = zeros(1,M);
    
    if(strcmp(Constellation_Type,'APSK-4-12'))
            
        if(strcmp(Lable_Type,'All-Gray'))
            
            phi_list = [0,-pi/12];
            r_list = [1.5:0.1:3];
            
        elseif(strcmp(Lable_Type,'All-SP'))
            
            phi_list = [0,-pi/12];
            r_list = [1.5:0.01:3];

        else
        
            phi_list = 0;
            r_list = 2.19;
                    
        end
        
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
            
        if(strcmp(Lable_Type,'All-Gray'))
            
            phi_list = [0,-pi/8];
            r_list = [1.5:0.1:3];
            
        elseif(strcmp(Lable_Type,'All-SP'))
            
            phi_list = [0,-pi/8];
            r_list = [1.5:0.01:3];

        else
        
            phi_list = -pi/8;
            r_list = 2.19;
                    
        end
       
    elseif(strcmp(Constellation_Type,'16APSK-2-Rings'))
            
        if(strcmp(Lable_Type,'All-Gray'))
            
            phi_list = [0,0];
            r_list = [1.5:0.1:3];
            
        elseif(strcmp(Lable_Type,'All-SP'))
            
            phi_list = [0,0];
            r_list = [1.5:0.01:3];

        else
        
            phi_list = [0];
            r_list = 2.19;
                    
        end
                
    end
    
    if(strcmp(Constellation_Type,'PSK'))
        if(m==2)
            Constellation_Type = 'QAM';
        end
    end
    
    if(strcmp(Lable_Type,'Regular'))
        Lable_Array = 0:1:M-1;
    elseif(strcmp(Lable_Type,'SP'))
        Lable_Array = SP_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'SP-GF4'))
        Lable_Array = SP_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'Gray'))
        if(strcmp(Constellation_Type,'16APSK-2-Rings'))
            Lable_Array = Gray_Labeling(Constellation_Type,m);
            inner_ring_sizes = 8;
        else
            Lable_Array = Gray_Labeling(Constellation_Type,m);
        end
    elseif(strcmp(Lable_Type,'Gray-GF4'))
            Lable_Array = Gray_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'Random'))
        Lable_Array = Random_Labeling(m);
    elseif(strcmp(Lable_Type,'SP-Compound'))
        Lable_Array = Compound_Labeling(m,g0,SP_Labeling(Constellation_Type,m));
    elseif(strcmp(Lable_Type,'Gray-Compound'))
        Lable_Array = Compound_Labeling(m,g0,Gray_Labeling(Constellation_Type,m));
    elseif(strcmp(Lable_Type,'SP-Compound-Interleaved'))
        Lable_Array = Compound_Interleaved_Labeling(m,g0,SP_Labeling(Constellation_Type,m));
    elseif(strcmp(Lable_Type,'Gray-Compound-Interleaved'))
        Lable_Array = Compound_Interleaved_Labeling(m,g0,Gray_Labeling(Constellation_Type,m));
    elseif(strcmp(Lable_Type,'Gray-Interleaved'))
        Lable_Array = Gray_Interleaved_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'SP-Interleaved'))
        Lable_Array = SP_Interleaved_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'Compound-BIPCM'))
        if(m==2)
            %Lable_Array = [0,2,1,3];
        elseif(m==4)
            %Lable_Array = [0,4,12,8,3,7,15,11,1,5,13,9,2,6,14,10];
        else
            return
        end
    elseif(strcmp(Lable_Type,'Compound-BIPCM-Interleaved'))
        if(m==2)
            %Lable_Array = [0,1,2,3];
        elseif(m==4)
            %Lable_Array = [0,4,3,7,12,8,15,11,1,5,2,6,13,9,14,10];
        else
            return
        end
    elseif(strcmp(Lable_Type,'SP-Gray'))
        if(strcmp(Constellation_Type,'QAM'))
            if(m==4)
                Lable_Array = [0,1,3,2,4,5,7,6,12,13,15,14,8,9,11,10];
            elseif(m==8)
                %Lable_Array = SP_Gray_256QAM();
            else
                return
            end
        elseif(strcmp(Constellation_Type,'PAM'))
            if(m==4)
                %Lable_Array = [0,9,3,10,1,8,2,11,15,14,4,13,6,7,5,12];   1             
                %Lable_Array = [0,9,3,10,1,8,2,11,15,6,4,5,14,7,13,12];   1  
                %Lable_Array = [0,9,11,10,1,8,2,3,15,14,12,13,6,7,5,4];   1
                %Lable_Array = [0,7,1,10,15,8,2,9,13,6,12,11,14,5,3,4];   2
                %Lable_Array = [0,15,13,14,7,8,6,5,1,10,12,11,2,9,3,4];   2
                
                %Lable_Array = [0,1,9,8,15,6,14,7,3,2,10,11,4,5,13,12];
                %Lable_Array = [0,1,9,8,3,6,14,11,15,2,10,7,4,5,13,12];
                %Lable_Array = [0,9,1,8,15,6,14,7,3,10,2,11,4,5,13,12];
                %Lable_Array = [0,9,1,8,3,10,2,11,15,14,6,7,4,13,5,12];
                %Lable_Array = [0,9,1,8,15,14,6,7,3,10,2,11,4,13,5,12];
                
                %Lable_Array = [0,1,9,8,3,2,10,11,15,6,14,7,4,5,13,12];
                
                %Lable_Array = [0,1,5,4,12,13,9,11,3,2,6,7,15,14,10,8];
                
                Lable_Array = [0,1,5,4,12,13,9,11,3,2,6,7,15,14,10,8];
                
            else
                return
            end
        end
    elseif(strcmp(Lable_Type,'SP-Gray-UVV-GF4'))
        if(strcmp(Constellation_Type,'PAM'))
            if(m==4)
                Lable_Array = [0,1,3,11,9,8,10,2,14,15,13,5,7,6,4,12];  
            else
                return
            end
        else
            return 
        end
    elseif(strcmp(Lable_Type,'SP-Gray-GF4'))
        if(strcmp(Constellation_Type,'PAM'))
            if(m==4)
                Lable_Array = [0,1,3,11,9,8,10,2,14,15,13,5,7,6,4,12];  
            else
                return
            end
        else
            return 
        end    
    elseif(strcmp(Lable_Type,'Gray-Permutated'))
        
        %Permutation = [1,2,3,4];
        %Permutation = [1,3,2,4];
        %Permutation = [4,2,3,1];
        %Permutation = [4,3,2,1];

        %Permutation = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
        
        Lable_Array = Gray_Permutation_Labeling(Constellation_Type,m,Permutation);
    elseif(strcmp(Lable_Type,'All'))
        Lable_Array = All_Labeling(m);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'All-Gray'))
        if(strcmp(Constellation_Type,'16APSK-2-Rings'))
            [Lable_Array,inner_ring_sizes] = All_Gray_Labeling_16APSK_2_Rings();
        else
            Lable_Array = All_Gray_Labeling(m,Constellation_Type);
        end
        if(strcmp(Constellation_Type,'APSK-4-12'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        elseif(strcmp(Constellation_Type,'APSK-8-8'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        elseif(strcmp(Constellation_Type,'16APSK-2-Rings'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        else
            Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
        end
    elseif(strcmp(Lable_Type,'All-Gray-GF4'))
        Lable_Array = All_Gray_GF4_Labeling(m,Constellation_Type);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'All-SP'))
        Lable_Array = All_SP_Labeling(m,Constellation_Type);
        if(strcmp(Constellation_Type,'APSK-4-12'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        elseif(strcmp(Constellation_Type,'APSK-8-8'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        elseif(strcmp(Constellation_Type,'16APSK-2-Rings'))
            Constellation_Mapping_Array = zeros(size(r_list,2)*size(phi_list,2)*size(Lable_Array,1),M);
        else
            Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
        end
    elseif(strcmp(Lable_Type,'All-SP-GF4'))
        Lable_Array = All_SP_GF4_Labeling(m,Constellation_Type);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'All-SP-Compound'))
        Lable_Array = All_SP_Compound_Labeling(m,Constellation_Type,g0);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'All-SP-Compound-GF4'))
        Lable_Array = All_SP_Compound_GF4_Labeling(m,Constellation_Type,g0);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'All-SP-Gray-GF4'))
        Lable_Array = All_SP_Gray_GF4_Labeling(m,Constellation_Type,g0);
        Constellation_Mapping_Array = zeros(size(Lable_Array,1),M);
    elseif(strcmp(Lable_Type,'Multi-Dimensional-SP'))
        Constellation_Mapping_Array = Multi_Dimensional_SP_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'Multi-Dimensional-Gray'))
        Constellation_Mapping_Array = Multi_Dimensional_Gray_Labeling(Constellation_Type,m);
    elseif(strcmp(Lable_Type,'Multi-Dimensional-SP-Compound'))
        Constellation_Mapping_Array = Multi_Dimensional_SP_Compound_Labeling(Constellation_Type,m);
    elseif(isnumeric(Lable_Type))
        Constellation_Mapping_Array = Lable_Type;
        return
    end
    
    if(strcmp(Constellation_Type,'PAM'))
        
        for j=1:1:size(Lable_Array,1)
            
            if(~isnumeric(Lable_Type))
            
                for i=0:1:M-1

                    Constellation_Mapping_Array(j,Lable_Array(j,i+1)+1) = 2*i-M+1;

                end
            
            end
        
            try
                
                if(size(Lable_Array,1)==1)

                    hMod = modem.pammod('M', M);
                    symbols = hMod.Constellation;
                    scatterplot(symbols)
                    set(get(gca,'Children'),'Marker','d','MarkerFaceColor', 'auto');
                    hold on;
                    for i=1:M
                       text(symbols(i)-m/2*0.1,1,dec2base(Lable_Array(j,i),2,m));
                    end
                    set(gca,'yTick',(-2:2:2),'xTick',(-(M-2):2:(M-2)), 'XLim',[-M-1,M+1],'YLim', [-2,2],'Box','on','YGrid','on', 'XGrid','on'); 
                end
                
            end
        end
    
    elseif(strcmp(Constellation_Type,'PSK'))
    
        for j=1:1:size(Lable_Array,1)
        
            if(~isnumeric(Lable_Type))
            
                for i=0:1:M-1

                    angle = (2*pi()/M);

                    Constellation_Mapping_Array(j,Lable_Array(j,i+1)+1) = cos(angle*i)+1i*sin(angle*i);

                end

            end
            
            try

                if(size(Lable_Array,1)==1)
                
                    delta = 10^-10;

                    hMod = modem.pskmod('M', M);
                    symbols = hMod.Constellation;
                    real_symbols = real(symbols);
                    imag_symbols = imag(symbols);

                    real_symbols(abs(real_symbols)<delta) = 0;
                    imag_symbols(abs(imag_symbols)<delta) = 0;

                    scatterplot(symbols)
                    set(get(gca,'Children'),'Marker','d','MarkerFaceColor', 'auto');
                    hold on;
                    for i=1:M
                       text(real_symbols(i)+(sign(real_symbols(i))-(sign(real_symbols(i))<0)-1)*m/2*0.025,imag_symbols(i)+sign(imag_symbols(i))*m/2*0.1, dec2base(Lable_Array(j,i),2,m));
                    end
                    set(gca,'yTick',(0),'xTick',(0), 'XLim',[-2,2],'YLim', [-2,2],'Box','on','YGrid','on', 'XGrid','on'); 
                
                end
                
            end
            
        end
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        for j=1:1:size(Lable_Array,1)
        
            if(~isnumeric(Lable_Type))
            
                for i=0:1:M-1

                    real_A = mod(i,2^(m/2));
                    imag_A = (i - real_A)/(2^(m/2));

                    Constellation_Mapping_Array(j,Lable_Array(j,i+1)+1) = (2*real_A-((2^(m/2))-1))+1i*(((2^(m/2))-1)-2*imag_A);

                end

            end
                
            try

                if(size(Lable_Array,1)==1)
                
                    hMod = modem.qammod('M', M);
                    symbols = hMod.Constellation;
                    scatterplot(symbols)
                    set(get(gca,'Children'),'Marker','d','MarkerFaceColor', 'auto');
                    hold on;
                    for i=1:M
                       text((-imag(symbols(i))-m/2*0.1),(-real(symbols(i))+m/2*0.1), dec2base(Lable_Array(j,i),2,m));
                    end
                    set(gca,'yTick',(-sqrt(M)+2:2:sqrt(M)-2),'xTick',(-sqrt(M)+2:2:sqrt(M)-2), 'XLim',[-sqrt(M)-1,sqrt(M)+1],'YLim', [-sqrt(M)-1,sqrt(M)+1],'Box','on','YGrid','on', 'XGrid','on'); 
                
                end
                
            end
            
        end
    
    elseif(strcmp(Constellation_Type,'APSK-4-12'))
        
        for r_index=1:1:size(r_list,2)
            
            r = r_list(r_index);
            
            for phi_index=1:1:size(phi_list,2)
                
                phi = phi_list(phi_index);
        
                for j=1:1:size(Lable_Array,1)

                    if(~isnumeric(Lable_Type))

                        for i=0:1:M-1

                            if(i<4)
                                angle = 2*pi()/4;
                                Constellation_Mapping_Array(((r_index-1)*size(phi_list,2)+(phi_index-1))*size(Lable_Array,1)+j,Lable_Array(j,i+1)+1) = cos(pi/4+angle*i)+1i*sin(pi/4+angle*i); % QPSK inner ring
                            else
                                angle = 2*pi()/12;
                                Constellation_Mapping_Array(((r_index-1)*size(phi_list,2)+(phi_index-1))*size(Lable_Array,1)+j,Lable_Array(j,i+1)+1) = r*(cos(pi/12+phi+angle*(i-4))+1i*sin(pi/12+phi+angle*(i-4))); % 12PSK outer ring
                            end

                        end

                    end

                    try
                        
                    end

                end
                
            end
            
        end
        
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
        
        angle = 2*pi()/8;
        
        for r_index=1:1:size(r_list,2)
            
            r = r_list(r_index);
            
            for phi_index=1:1:size(phi_list,2)
                
                phi = phi_list(phi_index);
        
                for j=1:1:size(Lable_Array,1)

                    if(~isnumeric(Lable_Type))
                        
                        for i=0:1:M-1

                            if(i<8)
                                Constellation_Mapping_Array(((r_index-1)*size(phi_list,2)+(phi_index-1))*size(Lable_Array,1)+j,Lable_Array(j,i+1)+1) = cos(angle*i)+1i*sin(angle*i); % 8PSK inner ring
                            else
                                Constellation_Mapping_Array(((r_index-1)*size(phi_list,2)+(phi_index-1))*size(Lable_Array,1)+j,Lable_Array(j,i+1)+1) = r*(cos(pi/8+phi+angle*(i-8))+1i*sin(pi/8+phi+angle*(i-8))); % 8PSK outer ring
                            end

                        end

                    end

                    try
                        
                    end

                end
                
            end
            
        end
        
    elseif(strcmp(Constellation_Type,'16APSK-2-Rings'))
        
        if(~isnumeric(Lable_Type))
        
            for r_index=1:1:size(r_list,2)

                r = r_list(r_index);

                for j=1:1:size(Lable_Array,1)

                    inner_ring_size = inner_ring_sizes(j);
                    outer_ring_size = M - inner_ring_size;

                    if(strcmp(Lable_Type,'All-Gray'))

                        phi_list = [0,-pi/outer_ring_size];

                    elseif(strcmp(Lable_Type,'All-SP'))

                        phi_list = [0,-pi/outer_ring_size];

                    else

                        phi_list = 0;

                    end

                    for phi_index=1:1:size(phi_list,2) 

                        phi = phi_list(phi_index);

                        for i=0:1:inner_ring_size-1
                            angle = 2*pi/inner_ring_size;
                            Constellation_Mapping_Array(((r_index-1)*size(Lable_Array,1)+(j-1))*size(phi_list,2)+phi_index,Lable_Array(j,i+1)+1) = cos(angle*i)+1i*sin(angle*i); % inner ring
                        end
                        
                        for i=inner_ring_size:1:M-1
                            angle = 2*pi/outer_ring_size;
                            Constellation_Mapping_Array(((r_index-1)*size(Lable_Array,1)+(j-1))*size(phi_list,2)+phi_index,Lable_Array(j,i+1)+1) = r*(cos(pi/outer_ring_size+phi+angle*(i-inner_ring_size))+1i*sin(pi/outer_ring_size+phi+angle*(i-inner_ring_size))); % outer ring
                        end

                        try

                        end

                    end

                end

            end
            
        end
    
    end
    
    if(strcmp(Constellation_Type,sprintf('%d-bit-2x8PSK',m)))
        
        Constellation_Mapping_Array = Constellation_Mapping_Array./sqrt(2);
        
    elseif(strcmp(Constellation_Type,sprintf('%d-bit-3x4PSK',m)))
        
        Constellation_Mapping_Array = Constellation_Mapping_Array./sqrt(3);
        
    else
        
        for i=1:1:size(Constellation_Mapping_Array,1)

            Constellation_Mapping_Array(i,:) = Constellation_Mapping_Array(i,:)./sqrt(sum(abs(Constellation_Mapping_Array(i,:)).^2)/M);

        end
        
    end
    
end

function [Lable_Array] = Gray_Labeling(Constellation_Type,m)

    if(strcmp(Constellation_Type,'APSK-4-12'))
        Lable_Array = [5,13,9,1,3,7,6,4,12,14,15,11,10,8,0,2];
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
        Lable_Array = [0,1,3,2,6,7,5,4,8,9,11,10,14,15,13,12];
    elseif(strcmp(Constellation_Type,'16APSK-2-Rings'))
        Lable_Array = [0,1,3,2,6,7,5,4,8,9,11,10,14,15,13,12];
    else
        
        M = 2^m;

        Lable_Array = zeros(1,M);

        for i=0:1:M-1
            Lable_Array(i+1) = bin2gray(i,lower(Constellation_Type),M);
        end 
        
    end

end

function [Lable_Array] = Gray_Interleaved_Labeling(Constellation_Type,m)

    M = 2^m;
    
    Lable_Array = zeros(1,M);
    
    for i=0:1:M-1

        i_string = dec2bin(i,m);

        for j=1:1:m
            i_bin(j) = i_string(j)=='1';
        end

        old_i_bin = zeros(1,length(i_bin(1,:)));
        old_i_bin(1,1:2:end) = i_bin(1,1:m/2);
        old_i_bin(1,2:2:end) = i_bin(1,m/2+1:m);
        
        old_i_int = 0;
        for j=1:1:m
            old_i_int = old_i_int + 2^(m-j)*old_i_bin(j);
        end
    
        Lable_Array(i+1) = bin2gray(old_i_int,lower(Constellation_Type),M);
        
    end 

end

function [Lable_Array] = Random_Labeling(m)

    M = 2^m;
    
    Lable_Array =  randperm(M)-1;

end

function [Lable_Array] = SP_Labeling(Constellation_Type,m)

    M = 2^m;
    
    Lable_Array = zeros(1,M);
    
    if(strcmp(Constellation_Type,'PSK'))
        Lable_Array = PSK_SP_Labeling_Recurtion(Lable_Array,m);
    elseif(strcmp(Constellation_Type,'PAM'))
        Lable_Array = PAM_SP_Labeling_Recurtion(Lable_Array,m);
    elseif(strcmp(Constellation_Type,'QAM'))
        Lable_Array = QAM_SP_Labeling_Recurtion(Lable_Array,m);
    elseif(strcmp(Constellation_Type,'APSK-4-12'))
        Lable_Array = [4,12,1,11,3,15,0,8,2,10,6,9,5,14,7,13];
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
        Lable_Array = [9,5,13,1,11,7,15,3,0,10,6,14,2,8,4,12];
    end
    
end

function [Lable_Array] = SP_Interleaved_Labeling(Constellation_Type,m)

    M = 2^m;
    
    Lable_Array_temp = zeros(1,M);
    
    if(strcmp(Constellation_Type,'PSK'))
        Lable_Array_temp = PSK_SP_Labeling_Recurtion(Lable_Array_temp,m);
    elseif(strcmp(Constellation_Type,'PAM'))
        Lable_Array_temp = PAM_SP_Labeling_Recurtion(Lable_Array_temp,m);
    elseif(strcmp(Constellation_Type,'QAM'))
        Lable_Array_temp = QAM_SP_Labeling_Recurtion(Lable_Array_temp,m);
    end
    
    for i=0:1:M-1

        i_string = dec2bin(i,m);

        for j=1:1:m
            i_bin(j) = i_string(j)=='1';
        end

        old_i_bin = zeros(1,length(i_bin(1,:)));
        old_i_bin(1,1:2:end) = i_bin(1,1:m/2);
        old_i_bin(1,2:2:end) = i_bin(1,m/2+1:m);
        
        old_i_int = 0;
        for j=1:1:m
            old_i_int = old_i_int + 2^(m-j)*old_i_bin(j);
        end
    
        Lable_Array(i+1) = find(Lable_Array_temp==old_i_int)-1;
        
    end 
    
end

function [Lable_Array] = PSK_SP_Labeling_Recurtion(Lable_Array,m)

    if(m==0)
        return;
    end

    M = 2^m;
            
    for i=1:1:M
        
        Lable_Array(i) = Lable_Array(i) + mod(i-1,2)*2^(m-1);
                    
    end
        
    Lable_Array1 = PSK_SP_Labeling_Recurtion(Lable_Array(1:2:end),m-1);
    Lable_Array2 = PSK_SP_Labeling_Recurtion(Lable_Array(2:2:end),m-1);
    
    Lable_Array = reshape([Lable_Array1;Lable_Array2],1,[]);
    
end

function [Lable_Array] = PAM_SP_Labeling_Recurtion(Lable_Array,m)

    if(m==0)
        return;
    end

    M = 2^m;
            
    for i=1:1:M
        
        Lable_Array(i) = Lable_Array(i) + mod(i-1,2)*2^(m-1);
                    
    end
        
    Lable_Array1 = PAM_SP_Labeling_Recurtion(Lable_Array(1:2:end),m-1);
    Lable_Array2 = PAM_SP_Labeling_Recurtion(Lable_Array(2:2:end),m-1);
    
    Lable_Array = reshape([Lable_Array1;Lable_Array2],1,[]);
    
end

function [Lable_Array] = QAM_SP_Labeling_Recurtion(Lable_Array,m)

    if(m==0)
        return;
    end

    M = 2^m;
        
    Lable_Matrix = reshape(Lable_Array,sqrt(M),sqrt(M))';
    
    for i=1:1:sqrt(M)
        
        for j=1:1:sqrt(M)
            
            Lable_Matrix(i,j) = Lable_Matrix(i,j) + mod(mod(j-1,2)+mod(i-1,2),2)*2^(m-1) + mod(i-1,2)*2^(m-2);
            
        end
        
    end
    
    for i=1:1:sqrt(M/4)
        
        for j=1:1:sqrt(M/4)
        
            Lable_Array1((i-1)*sqrt(M/4)+j) = Lable_Matrix((i-1)*2+1,(j-1)*2+1);
            Lable_Array2((i-1)*sqrt(M/4)+j) = Lable_Matrix((i-1)*2+1,(j-1)*2+2);
            Lable_Array3((i-1)*sqrt(M/4)+j) = Lable_Matrix((i-1)*2+2,(j-1)*2+1);
            Lable_Array4((i-1)*sqrt(M/4)+j) = Lable_Matrix((i-1)*2+2,(j-1)*2+2);

        end
        
    end
    
    Lable_Array1 = QAM_SP_Labeling_Recurtion(Lable_Array1,m-2);
    Lable_Array2 = QAM_SP_Labeling_Recurtion(Lable_Array2,m-2);
    Lable_Array3 = QAM_SP_Labeling_Recurtion(Lable_Array3,m-2);
    Lable_Array4 = QAM_SP_Labeling_Recurtion(Lable_Array4,m-2);
    
    for i=1:1:sqrt(M)
        
        for j=1:1:sqrt(M)
            
            if(mod(i,2))
                if(mod(j,2))
                    Lable_Matrix(i,j) = Lable_Array1(floor((i-1)/2)*sqrt(M/4)+floor(j/2)+1);
                else
                    Lable_Matrix(i,j) = Lable_Array2(floor((i-1)/2)*sqrt(M/4)+floor(j/2));
                end
            else
                if(mod(j,2))
                    Lable_Matrix(i,j) = Lable_Array3(floor((i-2)/2)*sqrt(M/4)+floor(j/2)+1);
                else
                    Lable_Matrix(i,j) = Lable_Array4(floor((i-2)/2)*sqrt(M/4)+floor(j/2));
                end
            end
            
        end
        
    end
    
    Lable_Array = reshape(Lable_Matrix',1,[]);
    
end

function [new_Lable_Array] = Compound_Interleaved_Labeling(m,g0,Lable_Array)

    M = 2^m;

    for i=0:1:M-1

        i_string = dec2bin(i,m);

        for j=1:1:m
            i_bin(j) = i_string(j)=='1';
        end

        i_bin = reshape(reshape(i_bin,[],2)',1,[]);
        
        old_i_bin = mod(i_bin*g0,2);
        old_i_bin = [old_i_bin(:,1:2:end),old_i_bin(:,2:2:end)];

        old_i_int = 0;
        for j=1:1:m
            old_i_int = old_i_int + 2^(m-j)*old_i_bin(j);
        end

        new_Lable_Array(i+1) = Lable_Array(old_i_int+1);
    end
   
end

function [new_Lable_Array] = Compound_Labeling(m,g0,Lable_Array)

    M = 2^m;

    for i=0:1:M-1

        i_string = dec2bin(Lable_Array(i+1),m);

        for j=1:1:m
            i_bin(j) = i_string(j)=='1';
        end

        old_i_bin = mod(i_bin*g0,2);

        old_i_int = 0;
        for j=1:1:m
            old_i_int = old_i_int + 2^(m-j)*old_i_bin(j);
        end

        new_Lable_Array(i+1) = old_i_int;
    end
   
end

function [Lable_Array] = Gray_Permutation_Labeling(Constellation_Type,m,Permutation)

    M = 2^m;
    
    Lable_Array = zeros(1,M);
    
    for i=0:1:M-1

        i_string = dec2bin(i,m);

        for j=1:1:m
            i_bin(j) = i_string(j)=='1';
        end

        old_i_bin = zeros(1,length(i_bin(1,:)));
        
        for k=1:1:m
            
            old_i_bin(k) = i_bin(Permutation(k));
        
        end
        
        old_i_int = 0;
        for j=1:1:m
            old_i_int = old_i_int + 2^(m-j)*old_i_bin(j);
        end
    
        Lable_Array(i+1) = bin2gray(old_i_int,lower(Constellation_Type),M);
        
    end 

end

function [ Gray_Labelings ] = All_Gray_Labeling_PSK(m,prev_Gray_Labelings)

    Gray_Labelings = All_Gray_Labeling_PAM(m,prev_Gray_Labelings);
    
end

function [ Gray_Labelings ] = All_Gray_GF4_Labeling_PSK(m,prev_Gray_Labelings)

    Gray_Labelings = All_Gray_GF4_Labeling_PAM(m,prev_Gray_Labelings);
    
end

function [ Gray_Labelings ] = All_Gray_Labeling_PAM(m,prev_Gray_Labelings)

    global Max_Num_Of_Labels;

    Index = find(prev_Gray_Labelings(1,:)==0,2);
    
    Index = Index(end);
    
    if(Index==1)
        Gray_Labelings = prev_Gray_Labelings;
        return;
    end
    
    M = 2^m;
    
    Labeling = prev_Gray_Labelings(1,1:Index-1);
    
    Counter = 0;

    if(Index==2)
        Unused_Numbers = 1:1:M-1;
    else
        Unused_Numbers = find(~sum(bsxfun(@eq,1:M-1,Labeling.')));
    end
    temp_bin = dec2bin([Labeling(Index-1),Unused_Numbers],m);
    Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    
    temp_zeros_vec = zeros(1,M-Index);
    
    for j=Good_Numbers
        
        temp_next_Gray_Labelings = All_Gray_Labeling_PAM(m,[Labeling,j,temp_zeros_vec]);
        
        if(isempty(temp_next_Gray_Labelings))
            continue;
        end
        
        temp_Gray_Labelings(Counter+1:Counter+size(temp_next_Gray_Labelings,1),:) = temp_next_Gray_Labelings;
        
        Counter = size(temp_Gray_Labelings,1);
        
        if(Counter>Max_Num_Of_Labels)
            break;
        end
            
    end
    
    if(Counter)
        Gray_Labelings = temp_Gray_Labelings;
        sprintf('Gray labelings finder finished symbol %d (out of %d) and already found %d candidates so far',Index,M,size(Gray_Labelings,1))
    else
        Gray_Labelings = [];
    end

end

function [ Gray_Labelings ] = All_Gray_GF4_Labeling_PAM(m,prev_Gray_Labelings)

    global Max_Num_Of_Labels;

    Index = find(prev_Gray_Labelings(1,:)==0,2);
    
    Index = Index(end);
    
    if(Index==1)
        Gray_Labelings = prev_Gray_Labelings;
        return;
    end
    
    M = 2^m;
    
    Labeling = prev_Gray_Labelings(1,1:Index-1);
    
    Counter = 0;

    if(Index==2)
        Unused_Numbers = 1:1:M-1;
    else
        Unused_Numbers = find(~sum(bsxfun(@eq,1:M-1,Labeling.')));
    end
    temp_quad = dec2base([Labeling(Index-1),Unused_Numbers],4,m/2);
    Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_quad(2:end,:),temp_quad(1,:)),2)==1);
    
    temp_zeros_vec = zeros(1,M-Index);
    
    for j=Good_Numbers
        
        temp_next_Gray_Labelings = All_Gray_GF4_Labeling_PAM(m,[Labeling,j,temp_zeros_vec]);
        
        if(isempty(temp_next_Gray_Labelings))
            continue;
        end
        
        temp_Gray_Labelings(Counter+1:Counter+size(temp_next_Gray_Labelings,1),:) = temp_next_Gray_Labelings;
        
        Counter = size(temp_Gray_Labelings,1);
        
        if(Counter>Max_Num_Of_Labels)
            break;
        end
            
    end
    
    if(Counter)
        Gray_Labelings = temp_Gray_Labelings;
        sprintf('Gray labelings finder finished symbol %d (out of %d) and already found %d candidates so far',Index,M,size(Gray_Labelings,1))
    else
        Gray_Labelings = [];
    end

end

function [ Gray_Labelings ] = All_Gray_Labeling_QAM(m,prev_Gray_Labelings)

    global Max_Num_Of_Labels;

    Index = find(prev_Gray_Labelings(1,:)==0,2);
    
    Index = Index(end);
    
    if(Index==1)
        Gray_Labelings = prev_Gray_Labelings;
        return;
    end
    
    M = 2^m;
    
    sqrt_M = sqrt(M);
    
    Is_Not_First_Column = mod(Index,sqrt_M)~=1;
    Is_Not_First_Row = Index>sqrt_M;
        
    Labeling = prev_Gray_Labelings(1,1:Index-1);
    
    Counter = 0;
    
    if(Index==2)
        Unused_Numbers = 1:1:M-1;
    else
        Unused_Numbers = find(~sum(bsxfun(@eq,1:M-1,Labeling.')));
    end
    
    if(and(Is_Not_First_Column,Is_Not_First_Row))
        temp_bin = dec2bin([Labeling(Index-1),Labeling(Index-sqrt(M)),Unused_Numbers],m);
        temp_Good_Numbers_Index = sum(bsxfun(@ne,temp_bin(3:end,:),temp_bin(1,:)),2)==1;
        Good_Numbers = Unused_Numbers(temp_Good_Numbers_Index);
        temp_bin = [temp_bin(2,:);temp_bin([0;0;temp_Good_Numbers_Index]==1,:)];
        Good_Numbers = Good_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    elseif(Is_Not_First_Row)
        temp_bin = dec2bin([Labeling(Index-sqrt(M)),Unused_Numbers],m);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    elseif(Is_Not_First_Column)
        temp_bin = dec2bin([Labeling(Index-1),Unused_Numbers],m);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    end
    
    for j=Good_Numbers
        
        temp_next_Gray_Labelings = All_Gray_Labeling_QAM(m,[Labeling,j,zeros(1,M-Index)]);
        
        if(isempty(temp_next_Gray_Labelings))
            continue;
        end
        
        temp_Gray_Labelings(Counter+1:Counter+size(temp_next_Gray_Labelings,1),:) = temp_next_Gray_Labelings;
        
        Counter = size(temp_Gray_Labelings,1);
        
        if(Counter>Max_Num_Of_Labels)
            break;
        end
        
    end
    
    if(Counter)
        Gray_Labelings = temp_Gray_Labelings;
        sprintf('Gray labelings finder finished symbol %d (out of %d) and already found %d candidates so far',Index,M,size(Gray_Labelings,1))
    else
        Gray_Labelings = [];
    end
    
end

function [ Gray_Labelings ] = All_Gray_GF4_Labeling_QAM(m,prev_Gray_Labelings)

    global Max_Num_Of_Labels;

    Index = find(prev_Gray_Labelings(1,:)==0,2);
    
    Index = Index(end);
    
    if(Index==1)
        Gray_Labelings = prev_Gray_Labelings;
        return;
    end
    
    M = 2^m;
    
    sqrt_M = sqrt(M);
    
    Is_Not_First_Column = mod(Index,sqrt_M)~=1;
    Is_Not_First_Row = Index>sqrt_M;
        
    Labeling = prev_Gray_Labelings(1,1:Index-1);
    
    Counter = 0;
    
    if(Index==2)
        Unused_Numbers = 1:1:M-1;
    else
        Unused_Numbers = find(~sum(bsxfun(@eq,1:M-1,Labeling.')));
    end
    
    if(and(Is_Not_First_Column,Is_Not_First_Row))
        temp_quad = dec2base([Labeling(Index-1),Labeling(Index-sqrt(M)),Unused_Numbers],4,m/2);
        temp_Good_Numbers_Index = sum(bsxfun(@ne,temp_quad(3:end,:),temp_quad(1,:)),2)==1;
        Good_Numbers = Unused_Numbers(temp_Good_Numbers_Index);
        temp_quad = [temp_quad(2,:);temp_quad([0;0;temp_Good_Numbers_Index]==1,:)];
        Good_Numbers = Good_Numbers(sum(bsxfun(@ne,temp_quad(2:end,:),temp_quad(1,:)),2)==1);
    elseif(Is_Not_First_Row)
        temp_quad = dec2base([Labeling(Index-sqrt(M)),Unused_Numbers],4,m/2);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_quad(2:end,:),temp_quad(1,:)),2)==1);
    elseif(Is_Not_First_Column)
        temp_quad = dec2base([Labeling(Index-1),Unused_Numbers],4,m/2);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_quad(2:end,:),temp_quad(1,:)),2)==1);
    end
    
    for j=Good_Numbers
        
        temp_next_Gray_Labelings = All_Gray_GF4_Labeling_QAM(m,[Labeling,j,zeros(1,M-Index)]);
        
        if(isempty(temp_next_Gray_Labelings))
            continue;
        end
        
        temp_Gray_Labelings(Counter+1:Counter+size(temp_next_Gray_Labelings,1),:) = temp_next_Gray_Labelings;
        
        Counter = size(temp_Gray_Labelings,1);
        
        if(Counter>Max_Num_Of_Labels)
            break;
        end
        
    end
    
    if(Counter)
        Gray_Labelings = temp_Gray_Labelings;
        sprintf('Gray labelings finder finished symbol %d (out of %d) and already found %d candidates so far',Index,M,size(Gray_Labelings,1))
    else
        Gray_Labelings = [];
    end
    
end

function [ SP_Labelings ] = All_SP_Labeling_PSK(m)

    SP_Labelings = All_SP_Labeling_PAM(m);

end

function [ SP_Labelings ] = All_SP_GF4_Labeling_PSK(m)

    SP_Labelings = All_SP_GF4_Labeling_PAM(m);

end

function [ SP_Labelings ] = All_SP_Labeling_PAM(m)

    if(m==0)
        return;
    end

    M = 2^m;
    
%     Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];
    
    Lable_Array0 = zeros(1,M);
    Lable_Array1 = zeros(1,M);
    Lable_Array2 = zeros(1,M);
    Lable_Array3 = zeros(1,M);
    Lable_Array4 = zeros(1,M);
    Lable_Array5 = zeros(1,M);
    Lable_Array6 = zeros(1,M);
    Lable_Array7 = zeros(1,M);

    for i=1:1:M

            Lable_Array0(1,i) = xor(0,mod(i-1,2))*2^(m-1) + xor(0,mod(i-1-mod(i-1,2),4))*2^(m-2);
            Lable_Array1(1,i) = xor(0,mod(i-1,2))*2^(m-1) + xor(1,mod(i-1-mod(i-1,2),4))*2^(m-2);
            Lable_Array2(1,i) = xor(1,mod(i-1,2))*2^(m-1) + xor(0,mod(i-1-mod(i-1,2),4))*2^(m-2);
            Lable_Array3(1,i) = xor(1,mod(i-1,2))*2^(m-1) + xor(1,mod(i-1-mod(i-1,2),4))*2^(m-2);
            
            Lable_Array4(1,i) = xor(0,mod(i-1,2))*2^(m-1) + xor(xor(0,mod(i-1-mod(i-1,2),4)),mod(i,2))*2^(m-2);
            Lable_Array5(1,i) = xor(0,mod(i-1,2))*2^(m-1) + xor(xor(1,mod(i-1-mod(i-1,2),4)),mod(i,2))*2^(m-2);
            Lable_Array6(1,i) = xor(1,mod(i-1,2))*2^(m-1) + xor(xor(0,mod(i-1-mod(i-1,2),4)),mod(i,2))*2^(m-2);
            Lable_Array7(1,i) = xor(1,mod(i-1,2))*2^(m-1) + xor(xor(1,mod(i-1-mod(i-1,2),4)),mod(i,2))*2^(m-2);

    end

    SP_Labelings = [Lable_Array0;Lable_Array1;Lable_Array2;Lable_Array3;Lable_Array4;Lable_Array5;Lable_Array6;Lable_Array7];
    
    if(m==2)
        return;
    end
    
    next_SP_Labelings = All_SP_Labeling_PAM(m-2);
    
    SP_Labelings = All_Combinations_PAM(m,[Lable_Array0;Lable_Array1;Lable_Array2;Lable_Array3;Lable_Array4;Lable_Array5;Lable_Array6;Lable_Array7],next_SP_Labelings);
  
end

function [ SP_Labelings ] = All_SP_GF4_Labeling_PAM(m)

    if(m==0)
        return;
    end

    M = 2^m;
    
    Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];
    
    Lable_Arrays = zeros(24,M);
    
    temp_Option_Vec_perms = perms(Option_Vec);
        
    for i=1:size(Option_Vec,2):M

        Lable_Arrays(:,i:i+size(Option_Vec,2)-1) = temp_Option_Vec_perms;
            
    end

    SP_Labelings = Lable_Arrays;
    
    if(m==2)
        return;
    end
    
    next_SP_Labelings = All_SP_GF4_Labeling_PAM(m-2);
    
    SP_Labelings = All_Combinations_GF4_PAM(m,Lable_Arrays,next_SP_Labelings);
  
end

function [ SP_Labelings ] = All_SP_Labeling_QAM(m)

    if(m==0)
        return;
    end

    M = 2^m;
    
%     Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];
    
    Lable_Matrix0 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix1 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix2 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix3 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix4 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix5 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix6 = zeros(sqrt(M),sqrt(M));
    Lable_Matrix7 = zeros(sqrt(M),sqrt(M));

    for i=1:1:sqrt(M)

        for j=1:1:sqrt(M)

            Lable_Matrix0(i,j) = xor(0,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(0,mod(i-1,2))*2^(m-2);
            Lable_Matrix1(i,j) = xor(0,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(1,mod(i-1,2))*2^(m-2);
            Lable_Matrix2(i,j) = xor(1,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(0,mod(i-1,2))*2^(m-2);
            Lable_Matrix3(i,j) = xor(1,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(1,mod(i-1,2))*2^(m-2);
            
            Lable_Matrix4(i,j) = xor(0,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(xor(0,xor(mod(i-1,2),mod(j,2))),mod(i,2))*2^(m-2);
            Lable_Matrix5(i,j) = xor(0,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(xor(1,xor(mod(i-1,2),mod(j,2))),mod(i,2))*2^(m-2);
            Lable_Matrix6(i,j) = xor(1,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(xor(0,xor(mod(i-1,2),mod(j,2))),mod(i,2))*2^(m-2);
            Lable_Matrix7(i,j) = xor(1,mod(mod(j-1,2)+mod(i-1,2),2))*2^(m-1) + xor(xor(1,xor(mod(i-1,2),mod(j,2))),mod(i,2))*2^(m-2);

        end

    end
        
    Lable_Array0 = reshape(Lable_Matrix0.',1,[]);
    Lable_Array1 = reshape(Lable_Matrix1.',1,[]);
    Lable_Array2 = reshape(Lable_Matrix2.',1,[]);
    Lable_Array3 = reshape(Lable_Matrix3.',1,[]);
    Lable_Array4 = reshape(Lable_Matrix4.',1,[]);
    Lable_Array5 = reshape(Lable_Matrix5.',1,[]);
    Lable_Array6 = reshape(Lable_Matrix6.',1,[]);
    Lable_Array7 = reshape(Lable_Matrix7.',1,[]);

    SP_Labelings = [Lable_Array0;Lable_Array1;Lable_Array2;Lable_Array3;Lable_Array4;Lable_Array5;Lable_Array6;Lable_Array7];
    
    if(m==2)
        return;
    end
    
    next_SP_Labelings = All_SP_Labeling_QAM(m-2);
    
    SP_Labelings = All_Combinations_QAM(m,[Lable_Array0;Lable_Array1;Lable_Array2;Lable_Array3;Lable_Array4;Lable_Array5;Lable_Array6;Lable_Array7],next_SP_Labelings);
  
end

function [ SP_Labelings ] = All_SP_GF4_Labeling_QAM(m)

    if(m==0)
        return;
    end

    M = 2^m;
    
    Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];
    
    Lable_Arrays = zeros(24,M);
    
    temp_Option_Vec_perms = perms(Option_Vec);
    
    for i=1:1:sqrt(M)/2
        
        for j=1:2:sqrt(M)
            
            Lable_Arrays(:,2*(i-1)*sqrt(M)+j) = temp_Option_Vec_perms(:,1);
            Lable_Arrays(:,2*(i-1)*sqrt(M)+j+1) = temp_Option_Vec_perms(:,2);
            Lable_Arrays(:,2*(i-1)*sqrt(M)+sqrt(M)+j) = temp_Option_Vec_perms(:,3);
            Lable_Arrays(:,2*(i-1)*sqrt(M)+sqrt(M)+j+1) = temp_Option_Vec_perms(:,4);
        
        end
        
    end
    
    SP_Labelings = Lable_Arrays;
    
    if(m==2)
        return;
    end
    
    next_SP_Labelings = All_SP_GF4_Labeling_QAM(m-2);
    
    SP_Labelings = All_Combinations_GF4_QAM(m,Lable_Arrays,next_SP_Labelings);
  
end

function [ All_Labelings ] = All_Combinations_PSK(m,SP_Labelings,next_SP_Labelings)

    All_Labelings = All_Combinations_PAM(m,SP_Labelings,next_SP_Labelings);

end

function [ All_Labelings ] = All_Combinations_GF4_PSK(m,SP_Labelings,next_SP_Labelings)

    All_Labelings = All_Combinations_GF4_PAM(m,SP_Labelings,next_SP_Labelings);

end

function [ All_Labelings ] = All_Combinations_PAM(m,SP_Labelings,next_SP_Labelings)

    All_Labelings = All_Combinations_QAM(m,SP_Labelings,next_SP_Labelings);

end

function [ All_Labelings ] = All_Combinations_GF4_PAM(m,SP_Labelings,next_SP_Labelings)

    All_Labelings = All_Combinations_GF4_QAM(m,SP_Labelings,next_SP_Labelings);

end

function [ All_Labelings ] = All_Combinations_QAM(m,SP_Labelings,next_SP_Labelings)

    M = 2^m;

    Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];

    Lable_Array0 = SP_Labelings(1,:);
    Lable_Array1 = SP_Labelings(2,:);
    Lable_Array2 = SP_Labelings(3,:);
    Lable_Array3 = SP_Labelings(4,:);
    Lable_Array4 = SP_Labelings(5,:);
    Lable_Array5 = SP_Labelings(6,:);
    Lable_Array6 = SP_Labelings(7,:);
    Lable_Array7 = SP_Labelings(8,:);

    permutations = Create_All_Permutation(M/4,log2(size(next_SP_Labelings,1)))+1;
    
    All_Labelings = [repmat(Lable_Array0,size(permutations,1),1);repmat(Lable_Array1,size(permutations,1),1);repmat(Lable_Array2,size(permutations,1),1);repmat(Lable_Array3,size(permutations,1),1);repmat(Lable_Array4,size(permutations,1),1);repmat(Lable_Array5,size(permutations,1),1);repmat(Lable_Array6,size(permutations,1),1);repmat(Lable_Array7,size(permutations,1),1)];
    
    for i=1:1:size(permutations,1)
        
        for j=1:1:size(Option_Vec,2)
            
            All_Labelings(size(permutations,1)*0+i,All_Labelings(size(permutations,1)*0+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*0+i,All_Labelings(size(permutations,1)*0+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*1+i,All_Labelings(size(permutations,1)*1+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*1+i,All_Labelings(size(permutations,1)*1+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*2+i,All_Labelings(size(permutations,1)*2+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*2+i,All_Labelings(size(permutations,1)*2+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*3+i,All_Labelings(size(permutations,1)*3+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*3+i,All_Labelings(size(permutations,1)*3+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*4+i,All_Labelings(size(permutations,1)*4+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*4+i,All_Labelings(size(permutations,1)*4+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*5+i,All_Labelings(size(permutations,1)*5+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*5+i,All_Labelings(size(permutations,1)*5+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*6+i,All_Labelings(size(permutations,1)*6+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*6+i,All_Labelings(size(permutations,1)*6+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            All_Labelings(size(permutations,1)*7+i,All_Labelings(size(permutations,1)*7+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*7+i,All_Labelings(size(permutations,1)*7+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            
        end
        
    end

end

function [ All_Labelings ] = All_Combinations_GF4_QAM(m,SP_Labelings,next_SP_Labelings)

    M = 2^m;

    Option_Vec = [0,2^(m-2),2^(m-1),2^(m-1)+2^(m-2)];

    Lable_Arrays = SP_Labelings;

    permutations = Create_All_Permutation(M/4,log2(size(next_SP_Labelings,1)))+1;
    
    All_Labelings = repmat(Lable_Arrays,size(permutations,1),1);
    
    for i=1:1:size(permutations,1)
        
        for j=1:1:size(Option_Vec,2)
            
            for k=0:1:size(Lable_Arrays,1)-1
            
                All_Labelings(size(permutations,1)*k+i,All_Labelings(size(permutations,1)*k+i,:)==Option_Vec(j)) = All_Labelings(size(permutations,1)*k+i,All_Labelings(size(permutations,1)*k+i,:)==Option_Vec(j)) + next_SP_Labelings(permutations(i,j),:);
            
            end
            
        end
        
    end

end

function [ SP_Compound_Labelings ] = All_SP_Compound_Labeling(m,Constellation_Type,g0)
    
    SP_Labelings = All_SP_Labeling(m,Constellation_Type);
    SP_Compound_Labelings = zeros(size(SP_Labelings));
    
    for i=1:1:size(SP_Labelings,1)
        SP_Compound_Labelings(i,:) = Compound_Labeling(m,g0,SP_Labelings(i,:));
    end

end


function [ SP_Gray_Labelings ] = All_SP_Gray_GF4_Labeling(m,Constellation_Type,g0)
    
    SP_Labelings = All_SP_GF4_Labeling(m,Constellation_Type);
    SP_Compound_Labelings = zeros(size(SP_Labelings));
    SP_Gray_Labelings = zeros(1,2^m);
    Index = 1;
    for i=1:1:size(SP_Labelings,1)
        SP_Compound_Labelings(i,:) = Compound_Labeling(m,g0,SP_Labelings(i,:));
        if(Is_Gray_Labeling_GF4(SP_Compound_Labelings(i,:),Constellation_Type))
            SP_Gray_Labelings(Index,:) = SP_Compound_Labelings(i,:);
            Index = Index +1;
        end
    end

end

function [ SP_Compound_Labelings ] = All_SP_Compound_GF4_Labeling(m,Constellation_Type,g0)
    
    SP_Labelings = All_SP_GF4_Labeling(m,Constellation_Type);
    SP_Compound_Labelings = zeros(size(SP_Labelings));
    
    for i=1:1:size(SP_Labelings,1)
        SP_Compound_Labelings(i,:) = Compound_Labeling(m,g0,SP_Labelings(i,:));
    end

end

function [ SP_Labelings ] = All_SP_Labeling_APSK_4_12()

    SP_Labelings = zeros(11,16);
    
%     SP_Labelings(0,:) = [0   ,1   ,2   ,3   ,4   ,5   ,6   ,7   ,8   ,9   ,10  ,11  ,12  ,13  ,14  ,15  ];
    SP_Labelings(1,:) = [7   ,15   ,4   ,12   ,1   ,14   ,3   ,8   ,6   ,10   ,0  ,13  ,2  ,9  ,5  ,11  ];
    SP_Labelings(2,:) = [5   ,13   ,6   ,14   ,7   ,9   ,2   ,15   ,0   ,10   ,4  ,8  ,3  ,12  ,1  ,11  ];
    SP_Labelings(3,:) = [7   ,15   ,0   ,8   ,5   ,10   ,1   ,13   ,4   ,9   ,3  ,12  ,6  ,11  ,2  ,14  ];
    SP_Labelings(4,:) = [5   ,13   ,4   ,12   ,1   ,14   ,2   ,9   ,7   ,10   ,0  ,15  ,3  ,8  ,6  ,11  ];
    SP_Labelings(5,:) = [5   ,13   ,4   ,12   ,7   ,9   ,6   ,15   ,0   ,11   ,6  ,8  ,2  ,14  ,1  ,10  ];
    SP_Labelings(6,:) = [6   ,14   ,3   ,11   ,4   ,9   ,2   ,12   ,7   ,10   ,0  ,15  ,5  ,8  ,1  ,13  ];
    SP_Labelings(7,:) = [5   ,7   ,4   ,6   ,0   ,15   ,11   ,3   ,12   ,9   ,1  ,14  ,10  ,2  ,13  ,8  ];
    SP_Labelings(8,:) = [5   ,13   ,7   ,15   ,6   ,11   ,1   ,14   ,2   ,9   ,4  ,10  ,0  ,12  ,3  ,8  ];
    SP_Labelings(9,:) = [2   ,10   ,1   ,9   ,6   ,11   ,4   ,14   ,0   ,12   ,7  ,8  ,5  ,15  ,3  ,13  ];
    SP_Labelings(10,:) = [9   ,11   ,8   ,10   ,5   ,1   ,13   ,6   ,2   ,15   ,4  ,0  ,12  ,7  ,3  ,14  ];
    SP_Labelings(11,:) = [4   ,12   ,1   ,11   ,3   ,15   ,0   ,8   ,2   ,10   ,6  ,9  ,5  ,14  ,7  ,13  ];
    
end

function [ SP_Labelings ] = All_SP_Labeling_APSK_8_8()

    SP_Labelings = zeros(9,16);
    
%     SP_Labelings(0,:) = [0   ,1   ,2   ,3   ,4   ,5   ,6   ,7   ,8   ,9   ,10  ,11  ,12  ,13  ,14  ,15  ];
    SP_Labelings(1,:) = [1   ,9   ,6   ,14   ,2   ,10   ,4   ,12   ,5   ,13   ,0  ,8  ,7  ,15  ,3  ,11  ];
    SP_Labelings(2,:) = [4   ,12   ,6   ,14   ,5   ,13   ,7   ,15   ,2   ,10   ,0  ,8  ,3  ,11  ,1  ,9  ];
    SP_Labelings(3,:) = [7   ,15   ,5   ,13   ,2   ,10   ,1   ,9   ,3   ,11   ,0  ,8  ,6  ,14  ,4  ,12  ];
    SP_Labelings(4,:) = [5   ,2   ,7   ,0   ,4   ,3   ,6   ,1   ,13   ,10   ,15  ,8  ,12  ,11  ,14  ,9  ];
    SP_Labelings(5,:) = [7   ,15   ,3   ,11   ,4   ,12   ,1   ,9   ,10   ,5   ,13  ,0  ,8  ,6  ,14  ,2  ];
    SP_Labelings(6,:) = [7   ,15   ,5   ,13   ,6   ,14   ,4   ,12   ,9   ,3   ,11  ,0  ,8  ,2  ,10  ,1  ];
    SP_Labelings(7,:) = [4   ,12   ,6   ,14   ,2   ,10   ,1   ,9   ,3   ,11   ,0  ,8  ,5  ,13  ,7  ,15  ];
    SP_Labelings(8,:) = [3   ,15   ,4   ,12   ,6   ,9   ,1   ,11   ,13   ,0   ,8  ,2  ,10  ,5  ,14  ,7  ];
    SP_Labelings(9,:) = [9,5,13,1,11,7,15,3,0,10,6,14,2,8,4,12];
    
end

function [ SP_Labelings ] = All_SP_Labeling(m,Constellation_Type)

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('SP_%d%s',2^m,Constellation_Type)));
        return;
    end

    if(strcmp(Constellation_Type,'PAM'))
    
        SP_Labelings = All_SP_Labeling_PAM(m);
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        SP_Labelings = All_SP_Labeling_PSK(m);
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        SP_Labelings = All_SP_Labeling_QAM(m);
        
    elseif(strcmp(Constellation_Type,'APSK-4-12'))
        
        SP_Labelings = All_SP_Labeling_APSK_4_12();
    
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
        
        SP_Labelings = All_SP_Labeling_APSK_8_8();
        
    end
    
    try
        save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('SP_%d%s',2^m,Constellation_Type)),'SP_Labelings');
    end

end

function [ SP_Labelings ] = All_SP_GF4_Labeling(m,Constellation_Type)

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('SP_GF4_%d%s',2^m,Constellation_Type)));
        return;
    end

    if(strcmp(Constellation_Type,'PAM'))
    
        SP_Labelings = All_SP_GF4_Labeling_PAM(m);
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        SP_Labelings = All_SP_GF4_Labeling_PSK(m);
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        SP_Labelings = All_SP_GF4_Labeling_QAM(m);
        
    end

    try
        save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('SP_GF4_%d%s',2^m,Constellation_Type)),'SP_Labelings');
    end
    
end

function [ Gray_Labelings ] = All_Gray_Labeling_APSK_By_Rings(M_list,prev_Gray_Labelings)

    global Max_Num_Of_Labels;

    Index = find(prev_Gray_Labelings(1,:)==0,2);
    
    Index = Index(end);
    
    if(Index==1)
        Gray_Labelings = prev_Gray_Labelings;
        return;
    end
    
    M = sum(M_list);
    
    m = log2(M);
    
    Labeling = prev_Gray_Labelings(1,1:Index-1);
    
    Counter = 0;

    if(Index==2)
        Unused_Numbers = 1:1:M-1;
    else
        Unused_Numbers = find(~sum(bsxfun(@eq,1:M-1,Labeling.')));
    end
    
    if(sum(Index==M_list+1)) %first symbol in a ring
        last_Ring_Index_Vec = cumsum(M_list);
        Ring_Index = find(Index>last_Ring_Index_Vec,1,'last');
        temp_bin = dec2bin([Labeling(last_Ring_Index_Vec(Ring_Index)-M_list(Ring_Index)+1),Unused_Numbers],m);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    else
        temp_bin = dec2bin([Labeling(Index-1),Unused_Numbers],m);
        Good_Numbers = Unused_Numbers(sum(bsxfun(@ne,temp_bin(2:end,:),temp_bin(1,:)),2)==1);
    end
    
    temp_zeros_vec = zeros(1,M-Index);
    
    for j=Good_Numbers
        
        temp_next_Gray_Labelings = All_Gray_Labeling_APSK_By_Rings(M_list,[Labeling,j,temp_zeros_vec]);
        
        if(isempty(temp_next_Gray_Labelings))
            continue;
        end
        
        temp_Gray_Labelings(Counter+1:Counter+size(temp_next_Gray_Labelings,1),:) = temp_next_Gray_Labelings;
        
        Counter = size(temp_Gray_Labelings,1);
        
        if(Counter>Max_Num_Of_Labels)
            break;
        end
            
    end
    
    if(Counter)
        Gray_Labelings = temp_Gray_Labelings;
        sprintf('Gray labelings finder finished symbol %d (out of %d) and already found %d candidates so far',Index,M,size(Gray_Labelings,1))
    else
        Gray_Labelings = [];
    end

end

function [ Gray_Labelings, inner_ring_sizes ] = All_Gray_Labeling_16APSK_2_Rings()

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_16APSK')));
        return;
    end

    M = 16;
    
    minimal_ring_size = 4;
    
    Gray_Labelings = zeros(0);
    inner_ring_sizes = zeros(0);

    for inner_ring_size = minimal_ring_size:1:M-minimal_ring_size

        temp_Gray_Labelings = All_Gray_Labeling_APSK_By_Rings([inner_ring_size,M-inner_ring_size],zeros(1,M));
        
        Gray_Labelings(end+1:end+size(temp_Gray_Labelings,1),:) = temp_Gray_Labelings;
        
        inner_ring_sizes(end+1:end+size(temp_Gray_Labelings,1),:) = inner_ring_size;
        
    end
    
    try
        save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_16APSK')));
    end
    
end

function [ Gray_Labelings ] = All_Gray_Labeling_APSK_4_12()

    Gray_Labelings = All_Gray_Labeling_APSK_By_Rings([4,12],zeros(1,16));

    Gray_Labelings(size(Gray_Labelings,1)+1,:) = [5   ,13   ,9   ,1   ,3   ,7   ,6   ,4   ,12   ,14   ,15  ,11  ,10  ,8  ,0  ,2  ];
    
    Gray_Labelings(size(Gray_Labelings,1)+1,:) = [3   ,11   ,9   ,1   ,7   ,6   ,2   ,10   ,14   ,15   ,13  ,12  ,8  ,0  ,1  ,5  ];
    
    Gray_Labelings(size(Gray_Labelings,1)+1,:) = [0   ,2   ,3   ,1   ,8   ,12   ,4   ,6   ,14   ,10   ,11  ,15  ,7  ,5  ,13  ,9  ];
    
end

function [ Gray_Labelings ] = All_Gray_Labeling_APSK_8_8()

    Gray_Labelings = All_Gray_Labeling_APSK_By_Rings([8,8],zeros(1,16));
    
    Gray_Labelings(size(Gray_Labelings,1)+1,:) = [0   ,1   ,3   ,2   ,6   ,7   ,5   ,4   ,8   ,9   ,11  ,10  ,14  ,15  ,13  ,12  ];

end

function [ Gray_Labelings ] = All_Gray_Labeling(m,Constellation_Type)

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_%d%s',2^m,Constellation_Type)));
        return;
    end
    
    if(strcmp(Constellation_Type,'PAM'))
    
        Gray_Labelings = All_Gray_Labeling_PAM(m,zeros(1,2^m));
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        Gray_Labelings = All_Gray_Labeling_PSK(m,zeros(1,2^m));
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        Gray_Labelings = All_Gray_Labeling_QAM(m,zeros(1,2^m));
        
    elseif(strcmp(Constellation_Type,'APSK-4-12'))
        
        Gray_Labelings = All_Gray_Labeling_APSK_4_12();
        
    elseif(strcmp(Constellation_Type,'APSK-8-8'))
        
        Gray_Labelings = All_Gray_Labeling_APSK_8_8();
        
    end
        
    try
        save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_%d%s',2^m,Constellation_Type)),'Gray_Labelings');
    end

end

function [ Gray_Labelings ] = All_Gray_GF4_Labeling(m,Constellation_Type)

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_GF4_%d%s',2^m,Constellation_Type)));
        return;
    end

    if(strcmp(Constellation_Type,'PAM'))
    
        Gray_Labelings = All_Gray_GF4_Labeling_PAM(m,zeros(1,2^m));
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        Gray_Labelings = All_Gray_GF4_Labeling_PSK(m,zeros(1,2^m));
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        Gray_Labelings = All_Gray_GF4_Labeling_QAM(m,zeros(1,2^m));
        
    end
        
    try
        save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('Gray_GF4_%d%s',2^m,Constellation_Type)),'Gray_Labelings');
    end

end

% function [ Labelings ] = All_Labeling(m)
% 
%     M = 2^m;
% 
%     Labelings = zeros(factorial(M-1),M);
% 
%     Labelings(:,2:end) = perms(1:1:M-1);
% 
% end

function [ Is_Gray ] = Is_Gray_Labeling(Labeling,Constellation_Type)

    if(strcmp(Constellation_Type,'PAM'))
    
        Is_Gray = Is_Gray_Labeling_PAM(Labeling);
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        Is_Gray = Is_Gray_Labeling_PSK(Labeling);
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        Is_Gray = Is_Gray_Labeling_QAM(Labeling);
        
    end

end

function [ Is_Gray ] = Is_Gray_Labeling_GF4(Labeling,Constellation_Type)

    if(strcmp(Constellation_Type,'PAM'))
    
        Is_Gray = Is_Gray_Labeling_PAM_GF4(Labeling);
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        Is_Gray = Is_Gray_Labeling_PSK_GF4(Labeling);
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        Is_Gray = Is_Gray_Labeling_QAM_GF4(Labeling);
        
    end

end

function [ Is_Gray ] = Is_Gray_Labeling_PAM(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    
    for i=1:1:M-1
        if(sum(dec2bin(Labeling(i),m) ~= dec2bin(Labeling(i+1),m))~=1)
            Is_Gray = false;
            break;
        end
    end
        
end

function [ Is_Gray ] = Is_Gray_Labeling_PAM_GF4(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    
    for i=1:1:M-1
        if(sum(dec2base(Labeling(i),4,m/2) ~= dec2base(Labeling(i+1),4,m/2))~=1)
            Is_Gray = false;
            break;
        end
    end
        
end

function [ Is_Gray ] = Is_Gray_Labeling_PSK(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    for i=1:1:M-1
        if(sum(dec2bin(Labeling(i),m) ~= dec2bin(Labeling(i+1),m))~=1)
            Is_Gray = false;
            break;
        end
    end

    if(sum(dec2bin(Labeling(1),m) ~= dec2bin(Labeling(end),m))~=1)
        Is_Gray = false;
    end
        
end

function [ Is_Gray ] = Is_Gray_Labeling_PSK_GF4(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    for i=1:1:M-1
        if(sum(dec2base(Labeling(i),4,m/2) ~= dec2base(Labeling(i+1),4,m/2))~=1)
            Is_Gray = false;
            break;
        end
    end

    if(sum(dec2base(Labeling(1),4,m/2) ~= dec2base(Labeling(end),4,m/2))~=1)
        Is_Gray = false;
    end
        
end

function [ Is_Gray ] = Is_Gray_Labeling_QAM(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    for i=1:1:sqrt(M)
        if(~Is_Gray)
            break;
        end
        for j=1:1:sqrt(M)-1
            if(sum(dec2bin(Labeling((i-1)*sqrt(M)+j),m) ~= dec2bin(Labeling((i-1)*sqrt(M)+j+1),m))~=1)
                Is_Gray = false;
                break;
            end
        end
    end

    for i=1:1:sqrt(M)
        if(~Is_Gray)
            break;
        end
        for j=1:1:sqrt(M)-1
            if(sum(dec2bin(Labeling((j-1)*sqrt(M)+i),m) ~= dec2bin(Labeling(j*sqrt(M)+i),m))~=1)
                Is_Gray = false;
                break;
            end
        end
    end
        
end

function [ Is_Gray ] = Is_Gray_Labeling_QAM_GF4(Labeling)

    Is_Gray = true;

    M = size(Labeling,2);

    m = log2(M);
    
    for i=1:1:sqrt(M)
        if(~Is_Gray)
            break;
        end
        for j=1:1:sqrt(M)-1
            if(sum(dec2base(Labeling((i-1)*sqrt(M)+j),4,m/2) ~= dec2base(Labeling((i-1)*sqrt(M)+j+1),4,m/2))~=1)
                Is_Gray = false;
                break;
            end
        end
    end

    for i=1:1:sqrt(M)
        if(~Is_Gray)
            break;
        end
        for j=1:1:sqrt(M)-1
            if(sum(dec2base(Labeling((j-1)*sqrt(M)+i),4,m/2) ~= dec2base(Labeling(j*sqrt(M)+i),4,m/2))~=1)
                Is_Gray = false;
                break;
            end
        end
    end
        
end

function [ Is_SP ] = Is_SP_Labeling(Labeling,Constellation_Type)
   
    if(strcmp(Constellation_Type,'PAM'))
    
        Is_SP = Is_SP_Labeling_PAM(Labeling);
        
    elseif(strcmp(Constellation_Type,'PSK'))
        
        Is_SP = Is_SP_Labeling_PSK(Labeling);
        
    elseif(strcmp(Constellation_Type,'QAM'))
        
        Is_SP = Is_SP_Labeling_QAM(Labeling);
        
    end

end

function [ Is_SP ] = Is_SP_Labeling_PAM(Labeling)

    Is_SP = true;

    M = size(Labeling,2);

    m = log2(M);
    
    if(m==1)
        return;
    end
    
    temp = dec2bin(Labeling(1),m) == dec2bin(M-1,m);

    first_bit = temp(1);

    for i=1:2:M
        lable = dec2bin(Labeling(i),m) == dec2bin(M-1,m);
        if(lable(1)~=first_bit)
            Is_SP = false;
            break;
        end
    end

    for i=2:2:M
        if(~Is_SP)
            break;
        end
        lable = dec2bin(Labeling(i),m) == dec2bin(M-1,m);
        if(lable(1)==first_bit)
            Is_SP = false;
            break;
        end
    end

    if(Is_SP)
        Is_SP = and(Is_SP_Labeling_PAM(Labeling(1:2:M)),Is_SP_Labeling_PAM(Labeling(2:2:M)));
    end

end

function [ Is_SP ] = Is_SP_Labeling_PSK(Labeling)

    Is_SP = true;

    M = size(Labeling,2);

    m = log2(M);
    
    if(m==1)
        return;
    end
        
    temp = dec2bin(Labeling(1),m) == dec2bin(M-1,m);

    first_bit = temp(1);

    for i=1:2:M
        lable = dec2bin(Labeling(i),m) == dec2bin(M-1,m);
        if(lable(1)~=first_bit)
            Is_SP = false;
            break;
        end
    end

    for i=2:2:M
        if(~Is_SP)
            break;
        end
        lable = dec2bin(Labeling(i),m) == dec2bin(M-1,m);
        if(lable(1)==first_bit)
            Is_SP = false;
            break;
        end
    end

    if(Is_SP)
        Is_SP = and(Is_SP_Labeling_PSK(Labeling(1:2:M)),Is_SP_Labeling_PSK(Labeling(2:2:M)));
    end

end

function [ Is_SP ] = Is_SP_Labeling_QAM(Labeling)

    Is_SP = true;

    M = size(Labeling,2);

    m = log2(M);
    
    if(m==1)
        return;
    end

    temp = dec2bin(Labeling(1),m) == dec2bin(M-1,m);

    first_bit = temp(1);

    Labeling_Matrix = reshape(Labeling,sqrt(M),sqrt(M)).';

    for i=1:1:sqrt(M)
        if(~Is_SP)
            break;
        end
        for j=1:1:sqrt(M)
            lable = dec2bin(Labeling_Matrix(i,j),m) == dec2bin(M-1,m);
            if(lable(1)~=xor(xor(first_bit,mod(i-1,2)),mod(j-1,2)))
                Is_SP = false;
                break;
            end
        end
    end

    Labeling1 = reshape([Labeling_Matrix(1:2:end,1:2:end),Labeling_Matrix(2:2:end,2:2:end)].',1,[]);
    Labeling1 = reshape(Labeling1,sqrt(M),sqrt(M)/2);
    Labeling2 = reshape([Labeling_Matrix(1:2:end,2:2:end),Labeling_Matrix(2:2:end,1:2:end)].',1,[]);
    Labeling2 = reshape(Labeling2,sqrt(M),sqrt(M)/2);

    temp1 = dec2bin(Labeling1(1),m) == dec2bin(M-1,m);
    temp2 = dec2bin(Labeling2(1),m) == dec2bin(M-1,m);

    second_bit1 = temp1(2);
    second_bit2 = temp2(2);

    for i=1:2:sqrt(M)
        if(~Is_SP)
            break;
        end
        for j=1:1:sqrt(M)/2
            lable1 = dec2bin(Labeling1(i,j),m) == dec2bin(M-1,m);
            lable2 = dec2bin(Labeling2(i,j),m) == dec2bin(M-1,m);
            if(or(lable1(2)~=second_bit1,lable2(2)~=second_bit2))
                Is_SP = false;
                break;
            end
        end
    end

    for i=2:2:sqrt(M)
        if(~Is_SP)
            break;
        end
        for j=1:1:sqrt(M)/2
            lable1 = dec2bin(Labeling1(i,j),m) == dec2bin(M-1,m);
            lable2 = dec2bin(Labeling2(i,j),m) == dec2bin(M-1,m);
            if(or(lable1(2)==second_bit1,lable2(2)==second_bit2))
                Is_SP = false;
                break;
            end
        end
    end

    if(Is_SP)
        Is_SP = and(Is_SP_Labeling_QAM(reshape(Labeling1.',1,[])),Is_SP_Labeling_QAM(reshape(Labeling2.',1,[])));
    end
        
end

function [x] = Create_All_Permutation(length,m)

    if(length==0)
        
        x = zeros(0);
        
        return
        
    end

    x = zeros((2^m)^length,length);

    if(length==1)
        
        x = transpose(0:1:(2^m)-1);
    
    else

        for i=1:1:(2^m)^length

            x(i,1) = floor((i-1)/((2^m)^(length-1)));

        end
        
        temp = Create_All_Permutation(length-1,m);
        
        for i=1:1:(2^m)
        
            x((i-1)*(2^m)^(length-1)+1:i*(2^m)^(length-1),2:end) = temp;
            
        end
    
    end

end

function [Lable_Array] = Multi_Dimensional_SP_Labeling(Constellation_Type,m)

    M = 2^m;

    if(strcmp(Constellation_Type,sprintf('%d-bit-2x8PSK',m)))
        
        Lable_Array = zeros(2,M);
        
        Constellation_8PSK = zeros(1,8);
        angle = 2*pi/8;
        for i=0:1:7
            Constellation_8PSK(1,i+1) = cos(angle*i)+1i*sin(angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_8PSK(0+1);Constellation_8PSK(2+1)];
            Lable_Array(:,2) = [Constellation_8PSK(4+1);Constellation_8PSK(6+1)];
            Lable_Array(:,3) = [Constellation_8PSK(2+1);Constellation_8PSK(4+1)];
            Lable_Array(:,4) = [Constellation_8PSK(6+1);Constellation_8PSK(0+1)];
            Lable_Array(:,5) = [Constellation_8PSK(0+1);Constellation_8PSK(6+1)];
            Lable_Array(:,6) = [Constellation_8PSK(4+1);Constellation_8PSK(2+1)];
            Lable_Array(:,7) = [Constellation_8PSK(2+1);Constellation_8PSK(0+1)];
            Lable_Array(:,8) = [Constellation_8PSK(6+1);Constellation_8PSK(4+1)];
            Lable_Array(:,9) = [Constellation_8PSK(1+1);Constellation_8PSK(3+1)];
            Lable_Array(:,10) = [Constellation_8PSK(5+1);Constellation_8PSK(7+1)];
            Lable_Array(:,11) = [Constellation_8PSK(3+1);Constellation_8PSK(5+1)];
            Lable_Array(:,12) = [Constellation_8PSK(7+1);Constellation_8PSK(1+1)];
            Lable_Array(:,13) = [Constellation_8PSK(1+1);Constellation_8PSK(7+1)];
            Lable_Array(:,14) = [Constellation_8PSK(5+1);Constellation_8PSK(3+1)];
            Lable_Array(:,15) = [Constellation_8PSK(3+1);Constellation_8PSK(1+1)];
            Lable_Array(:,16) = [Constellation_8PSK(7+1);Constellation_8PSK(5+1)];
        end
        
    elseif(strcmp(Constellation_Type,sprintf('%d-bit-3x4PSK',m)))
        
        Lable_Array = zeros(3,M);
        
        Constellation_4PSK = zeros(1,4);
        angle = 2*pi/4;
        phi = pi/4;
        for i=0:1:3
            Constellation_4PSK(1,i+1) = cos(phi+angle*i)+1i*sin(phi+angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,2) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
            Lable_Array(:,3) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,4) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
            Lable_Array(:,5) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,6) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,7) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,8) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,9) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,10) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,11) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,12) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,13) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
            Lable_Array(:,14) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,15) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,16) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
        end
        
    end
        
end

function [Lable_Array] = Multi_Dimensional_Gray_Labeling(Constellation_Type,m)

    M = 2^m;

    if(strcmp(Constellation_Type,sprintf('%d-bit-2x8PSK',m)))
        
        Lable_Array = zeros(2,M);
        
        Constellation_8PSK = zeros(1,8);
        angle = 2*pi/8;
        for i=0:1:7
            Constellation_8PSK(1,i+1) = cos(angle*i)+1i*sin(angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_8PSK(0+1);Constellation_8PSK(2+1)];
            Lable_Array(:,2) = [Constellation_8PSK(0+1);Constellation_8PSK(6+1)];
            Lable_Array(:,3) = [Constellation_8PSK(4+1);Constellation_8PSK(2+1)];
            Lable_Array(:,4) = [Constellation_8PSK(4+1);Constellation_8PSK(6+1)];
            Lable_Array(:,5) = [Constellation_8PSK(2+1);Constellation_8PSK(0+1)];
            Lable_Array(:,6) = [Constellation_8PSK(2+1);Constellation_8PSK(4+1)];
            Lable_Array(:,7) = [Constellation_8PSK(6+1);Constellation_8PSK(0+1)];
            Lable_Array(:,8) = [Constellation_8PSK(6+1);Constellation_8PSK(4+1)];
            Lable_Array(:,9) = [Constellation_8PSK(3+1);Constellation_8PSK(1+1)];
            Lable_Array(:,10) = [Constellation_8PSK(3+1);Constellation_8PSK(5+1)];
            Lable_Array(:,11) = [Constellation_8PSK(7+1);Constellation_8PSK(1+1)];
            Lable_Array(:,12) = [Constellation_8PSK(7+1);Constellation_8PSK(5+1)];
            Lable_Array(:,13) = [Constellation_8PSK(1+1);Constellation_8PSK(3+1)];
            Lable_Array(:,14) = [Constellation_8PSK(1+1);Constellation_8PSK(7+1)];
            Lable_Array(:,15) = [Constellation_8PSK(5+1);Constellation_8PSK(3+1)];
            Lable_Array(:,16) = [Constellation_8PSK(5+1);Constellation_8PSK(7+1)];
        end
        
    elseif(strcmp(Constellation_Type,sprintf('%d-bit-3x4PSK',m)))
        
        Lable_Array = zeros(3,M);
        
        Constellation_4PSK = zeros(1,4);
        angle = 2*pi/4;
        phi = pi/4;
        for i=0:1:3
            Constellation_4PSK(1,i+1) = cos(phi+angle*i)+1i*sin(phi+angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,2) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,3) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,4) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
            Lable_Array(:,5) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,6) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,7) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,8) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
            Lable_Array(:,9) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,10) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,11) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,12) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
            Lable_Array(:,13) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,14) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,15) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,16) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
        end
        
    end
        
end

function [Lable_Array] = Multi_Dimensional_SP_Compound_Labeling(Constellation_Type,m)

    M = 2^m;

    if(strcmp(Constellation_Type,sprintf('%d-bit-2x8PSK',m)))
        
        Lable_Array = zeros(2,M);
        
        Constellation_8PSK = zeros(1,8);
        angle = 2*pi/8;
        for i=0:1:7
            Constellation_8PSK(1,i+1) = cos(angle*i)+1i*sin(angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_8PSK(0+1);Constellation_8PSK(2+1)];
            Lable_Array(:,2) = [Constellation_8PSK(7+1);Constellation_8PSK(5+1)];
            Lable_Array(:,3) = [Constellation_8PSK(3+1);Constellation_8PSK(5+1)];
            Lable_Array(:,4) = [Constellation_8PSK(4+1);Constellation_8PSK(2+1)];
            Lable_Array(:,5) = [Constellation_8PSK(1+1);Constellation_8PSK(7+1)];
            Lable_Array(:,6) = [Constellation_8PSK(6+1);Constellation_8PSK(0+1)];
            Lable_Array(:,7) = [Constellation_8PSK(2+1);Constellation_8PSK(0+1)];
            Lable_Array(:,8) = [Constellation_8PSK(5+1);Constellation_8PSK(7+1)];
            Lable_Array(:,9) = [Constellation_8PSK(1+1);Constellation_8PSK(3+1)];
            Lable_Array(:,10) = [Constellation_8PSK(6+1);Constellation_8PSK(4+1)];
            Lable_Array(:,11) = [Constellation_8PSK(2+1);Constellation_8PSK(4+1)];
            Lable_Array(:,12) = [Constellation_8PSK(5+1);Constellation_8PSK(3+1)];
            Lable_Array(:,13) = [Constellation_8PSK(0+1);Constellation_8PSK(6+1)];
            Lable_Array(:,14) = [Constellation_8PSK(7+1);Constellation_8PSK(1+1)];
            Lable_Array(:,15) = [Constellation_8PSK(3+1);Constellation_8PSK(1+1)];
            Lable_Array(:,16) = [Constellation_8PSK(4+1);Constellation_8PSK(6+1)];
        end
        
    elseif(strcmp(Constellation_Type,sprintf('%d-bit-3x4PSK',m)))
        
        Lable_Array = zeros(3,M);
        
        Constellation_4PSK = zeros(1,4);
        angle = 2*pi/4;
        phi = pi/4;
        for i=0:1:3
            Constellation_4PSK(1,i+1) = cos(phi+angle*i)+1i*sin(phi+angle*i);
        end
        
        if(m==4)
            Lable_Array(:,1) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,2) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
            Lable_Array(:,3) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,4) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,5) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
            Lable_Array(:,6) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(2+1)];
            Lable_Array(:,7) = [Constellation_4PSK(0+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,8) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,9) = [Constellation_4PSK(1+1);Constellation_4PSK(0+1);Constellation_4PSK(3+1)];
            Lable_Array(:,10) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(2+1)];
            Lable_Array(:,11) = [Constellation_4PSK(0+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,12) = [Constellation_4PSK(3+1);Constellation_4PSK(0+1);Constellation_4PSK(1+1)];
            Lable_Array(:,13) = [Constellation_4PSK(1+1);Constellation_4PSK(2+1);Constellation_4PSK(1+1)];
            Lable_Array(:,14) = [Constellation_4PSK(2+1);Constellation_4PSK(3+1);Constellation_4PSK(0+1)];
            Lable_Array(:,15) = [Constellation_4PSK(2+1);Constellation_4PSK(1+1);Constellation_4PSK(0+1)];
            Lable_Array(:,16) = [Constellation_4PSK(3+1);Constellation_4PSK(2+1);Constellation_4PSK(3+1)];
        end
        
    end
        
end