% calculating energy of MLM + Analog PPM
close all;

beta = 5:1:1e5;
opt_delta = zeros(size(beta));

for i=1:length(beta)
   
    delta = linspace((1/beta(i))^3,0.3,3*1e4);
    [min_val,delta_idx] = min(abs(delta.*sqrt(4*pi*log(delta*(beta(i))^3)) - pi^4/(6^3)));
    
    if min_val > 1e-4
        disp('WTF');
    end
    
    opt_delta(i) = delta(delta_idx);
end

figure;hold all;
plot(1./(log(beta)),opt_delta,'LineWidth',2);
plot(1./(log(beta)),(0.1)./(log(beta)) + 0.02,'LineWidth',2);
grid on; grid minor;
legend('numerical','analytical');
xlabel('\beta'); ylabel('optimal delta');