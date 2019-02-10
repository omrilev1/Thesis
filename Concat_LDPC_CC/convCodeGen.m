function [trellis,tblLen] = convCodeGen(v,r)
% convCode is a struct, where input contain the constraint length
% and rate
% the output include the trellis, and viterbi decoder memory length
tblLen = 4*v;

switch v
    
    case 4
        
        switch r
            
            case 1/4
                trellis = [37,33,27,25];
            case 1/3
                trellis = [37,33,25];
            case 2/5
                trellis = [21,27,0;37,33,25];
            case 1/2
                trellis = [31,27];
        end
        
    case 5
        
        switch r
            
            case 1/4
                trellis = [77,67,55,51];
            case 1/3
                trellis = [75,53,47];
            case 2/5
                trellis = [55,73,0;73,75,51];
            case 1/2
                trellis = [65,57];
        end
    case 6
        
        switch r
            
            case 1/4
                trellis = [171,155,127,117];
            case 1/3
                trellis = [155,127,117];
            case 2/5
                trellis = [157,165,0;157,127,111];
            case 1/2
                trellis = [155,117];
                
        end
    case 7
        
        switch r
            
            case 1/4
                trellis = [353,335,277,231];
            case 1/3
                trellis = [357,251,233];
            case 2/5
                trellis = [333,275,0;357,245,235];
            case 1/2
                trellis = [345,237];
        end
        
    case 8
        
        switch r
            
            case 1/4
                trellis = [765,671,513,473];
            case 1/3
                trellis = [637,567,515];
            case 2/5
                trellis = [633,565,0;576,613,571];
            case 1/2
                trellis = [657,435];
        end
end
end

