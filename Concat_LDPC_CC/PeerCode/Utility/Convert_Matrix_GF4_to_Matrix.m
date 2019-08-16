function [g] = Convert_Matrix_GF4_to_Matrix(g_gf4,kernel_name)

    try
        load(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('g_%s_N=%d[bits]',kernel_name,2*size(g_gf4,1))));
        return;
    end

    zeroMat = zeros(2,2);
    oneMat = [1 0 ; 0 1];
    alphaMat = [0 1 ; 1 1 ];
    alpha2Mat = mod(alphaMat^2,2);
    
    zero = 0;
    one = 1;
    alpha = 2;
    alpha2 = 3;
    
    g = zeros(2*size(g_gf4));
    
    for row=1:1:size(g_gf4,1)
         for column=1:1:size(g_gf4,2)
             if(g_gf4(row,column) == zero)
                 continue;
             elseif(g_gf4(row,column) == one)
                 g(2*row-1:2*row,2*column-1:2*column) = oneMat;
             elseif(g_gf4(row,column) == alpha)
                 g(2*row-1:2*row,2*column-1:2*column) = alphaMat;
             elseif(g_gf4(row,column) == alpha2)
                 g(2*row-1:2*row,2*column-1:2*column) = alpha2Mat;
             end
        end
    end
    
    save(fullfile(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Memory'),sprintf('g_%s_N=%d[bits]',kernel_name,2*size(g_gf4,1))),'g');

end

