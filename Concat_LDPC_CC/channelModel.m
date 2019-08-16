function [LLR,LL0,LL1] = channelModel(x,cap,channelType)

switch channelType
    
    case 2
        SNR = 10*log10(2^(2*cap) - 1);
        sigma = (1/sqrt(2)) * 10^(-SNR/20);
        y = x + sigma * randn(size(x));
        LLR = 2*real(y)/(sigma^2);
        
        temp_2var = 2*(10^(-SNR/20))^2;
        temp_add = - 0.5*log(pi()*temp_2var);
        LL0 = temp_add - ((y-1).^2)./temp_2var;
        LL1 = temp_add - ((y+1).^2)./temp_2var;  
    
    case 1
        p = entropyInverse(1 - cap);
        channelPattern = (rand(1,length(x)) > (1-p));
        y = mod((x(:)+1)/2 + channelPattern(:),2);
        
        LLR = zeros(size(y));LL0 = zeros(size(y));LL1 = zeros(size(y));
        LLR(y == 1) = log((1-p)/p);
        LLR(y == 0) = -1*log((1-p)/p);
        
        LL0(y == 0) = log(1-p);
        LL0(y == 1) = log(p);
        
        LL1(y == 1) = log(1-p);
        LL1(y == 0) = log(p);

    case 0
        p = 1 - cap;
        channelPattern = (rand(1,length(x)) > (1-p));
        LLR = 100*sign(x);
        LLR(channelPattern == 1) = 0;

        LL0(channelPattern == 1) = log(1/2);
        LL0((channelPattern == 0) & (x == 1)) = -100;
        LL0((channelPattern == 0) & (x == 0)) = 0;
        
        LL1(channelPattern == 1) = log(1/2);
        LL1((channelPattern == 0) & (x == 1)) = 0;
        LL1((channelPattern == 0) & (x == 0)) = -100;
        
end
LLR = reshape(LLR,1,[]);
LL0 = reshape(LL0,1,[]);
LL1 = reshape(LL1,1,[]);

function [y] = entropyInverse(x)

arg = 0.001:0.001:(1-0.001);
H = -1*(arg.*log2(arg) + (1-arg).*log2(1-arg));

% find the nearest value to x
[~,idx_min] = min(abs(H - x));
y = arg(idx_min);

