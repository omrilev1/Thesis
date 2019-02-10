function [optAlpha_Tx,optAlpha_Rx] = solve_TwoAlpha_Optimization(snr,sir)

dt = 1e-4;
alpha = 0:0.01:1;

f_alpha = zeros(size(alpha));
g_alpha = zeros(size(alpha));
obj_func = zeros(size(alpha));

snrLin = 10.^(snr/10);
sirLin = 10.^(sir/10);

for i=1:length(alpha)
    
    t1 = alpha(i)/2 : dt : 0.5;
    t2 = (alpha(i)/8) * (1:dt:2).^2;% alpha(i)/2 * (0.5:dt:1).^2;
    f_alpha(i) = 2*sum(H_b(t1)) * dt + alpha(i) * sum(H_b(t2)) * dt;
    
    alpha_Rx = snrLin * (1 + alpha(i)/sirLin) / (1 + (1 + (1/sirLin))*snrLin);
    g_alpha(i) = 0.5*log2((alpha_Rx)^2 + (alpha_Rx-1)^2 * snrLin + (alpha_Rx - alpha(i))^2 * snrLin/sirLin);
    
    obj_func(i) = g_alpha(i) + f_alpha(i);
end
[~,opt_idx] = min(obj_func);
optAlpha_Tx = alpha(opt_idx);
optAlpha_Rx = snrLin * (1 + alpha(opt_idx) /sirLin) / (1 + (1 + (1/sirLin))*snrLin);

end


function [Hb] = H_b(x)
Hb = -x.*log2(x) - (1-x).*log2(1-x);
Hb(x == 0) = 0;
end
