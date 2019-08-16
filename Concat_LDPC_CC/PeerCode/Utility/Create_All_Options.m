function [x] = Create_All_Options(length,m)

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
        
        temp = Create_All_Options(length-1,m);
        
        for i=1:1:(2^m)
        
            x((i-1)*(2^m)^(length-1)+1:i*(2^m)^(length-1),2:end) = temp;
            
        end
    
    end

end