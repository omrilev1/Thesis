% Simulation of robust Posterior Matching using Horstein scheme
clear all; clc;

p = 0.11;
capacity = binaryEntropy(p);

eps = p;
epsFB = [0 0.01 0.05];
k = 4;

dtheta = 2^(-(k + 5));
theta_axis = 0:dtheta:(1 - dtheta);
intervalSize = round(2^(-k) / dtheta);
pdf_init = dtheta * ones(size(theta_axis));

N = (floor((k/capacity))) : 2 : (ceil((k/capacity)) + 36);

BLER = zeros(1,length(N));
nIter = 1e4;


% statistics loop
for j=1:length(epsFB)
    for i=1:length(N)
        
        
        currBLER = 0;currBLER1 = 0;
        for iter = 1 : nIter
            
            interval = round(2^k * rand);
            theta0 = interval/2^k + 2^(-(k + 1));
            
            for n=1:N(i)
                
                if n==1
                    pdf_arr = pdf_init;
                else
                    if yFB == 1
                        pdf_arr = 2*p*pdf_arr;
                        pdf_arr(theta_axis > curr_med) = 2*(1-p)*pdf_arr(theta_axis > curr_med);
                    else
                        pdf_arr = 2*p*pdf_arr;
                        pdf_arr(theta_axis < curr_med) = 2*(1-p)*pdf_arr(theta_axis < curr_med);
                    end
                end
                pdf_arr(:) = pdf_arr(:)/sum(pdf_arr(:));
                
                % curr_med = calc_med(pdf_arr(:),theta_axis);
                
                % calc median
                f_fix = pdf_arr(:)/sum(pdf_arr(:));
                cum_f = cumsum(f_fix);
                [~,med_idx] = min(abs(cum_f - 1/2));
                curr_med = theta_axis(med_idx);
                
                if theta0 > curr_med
                    x = 1;
                else
                    x = 0;
                end
                
                noise = (rand > (1-eps));
                noiseFB = (rand > (1-epsFB(j)));
                y = mod(x + noise,2);
                yFB = mod(y + noiseFB,2);
                
            end
            % decode maximum aposteriori interval
            map = sum(reshape(pdf_arr,intervalSize,[]),1);
            [~,theta_hat] = max(map);
            theta_hat1 = theta_hat - 1;
            err1 = (theta_hat1 ~= interval);
            err0 = (theta_hat ~= interval);
            err = err0 + err1;
            
            currBLER = currBLER + (err > 1);
        end
        BLER(j,i) = currBLER/nIter;
        disp(strcat('Finished N = ',num2str(N(i)),' BLER = ',num2str(BLER(1,i)),' eps_{FB} = ',num2str(epsFB(j))));
        
    end
end
figure;
currLegend = {};
for idx = 1 : length(epsFB)
    semilogy(k./N,BLER(idx,:),'LineWidth',2); hold on;
    currLegend = [currLegend;strcat('\eps_{FB} = ',num2str(epsFB(idx)))];
end
semilogy((capacity)*ones(1,10),linspace(1e-4,1,10),'--','LineWidth',2);
grid on; grid minor;
xlabel('R'); legend(currLegend,'Capacity');
title(strcat('BLER Vs Rate, p = ',num2str(p),' \epsilon = ',num2str(eps),' k = ',num2str(k)));
xlim([min(k./N),1.01*capacity]);

function H = binaryEntropy(p)
if p > 0.5
    p = 1 - p;
end
if p == 0
    H = 1;
else
    H = -p*log2(p) - (1-p)*log2(1-p);
end
end