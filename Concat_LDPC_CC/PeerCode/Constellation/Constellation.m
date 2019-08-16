function [ X_after_constellation ] = Constellation(m,Constellation_Mapping_Array,X_after_encoding)

    N = size(X_after_encoding,2);
    N_div_m = N/m;
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);

    X_before_constellation = zeros(1,N_div_m);
    X_after_constellation = zeros(Constellation_2D_Dimension,N_div_m);

    %Change alphabet into constellation alphabet
    for i=1:N_div_m
        X_before_constellation(i) = 0;
        for j=0:1:m-1
            X_before_constellation(i) = X_before_constellation(i)+X_after_encoding(i*m-j)*(2^j);
        end
    end

%------------------------------constellation------------------------------%
    %Transform the output bits into output symbols
    for i=1:N_div_m
         X_after_constellation(:,i) = Constellation_Mapping_Array(:,X_before_constellation(i)+1);
    end

    X_after_constellation = reshape(X_after_constellation,1,[]);
    
end

