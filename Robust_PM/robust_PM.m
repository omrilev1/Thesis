% Simulation of robust Posterior Matching using Horstein scheme
close all; clear all; clc;

dtheta = 0.8*1e-6;
theta_axis = 0:dtheta:1;

pdf_init = dtheta * ones(size(theta_axis));

p = 0.11;
capacity = binaryEntropy(p);

eps = p;
eps_feedback = 0.01;
k = 8;
pdf_arr = zeros(length(pdf_init),N);
pdf_arr(:,1) = pdf_init;

curr_med = zeros(1,N);
med_lower_bound = zeros(1,N);
med_upper_bound = zeros(1,N);
noise = zeros(1,N);
theta0 = round(2^k * rand)/2^k + 2^(-(k + 1));

N = ceil(k*capacity) : 2 : 30 * ceil(k*capacity);

delta = 0.2;

Niter = 1e3;

BER = zeros(1,N);
nIter = 1e4;
% statistics loop
for i=1:length(N)
    
    for iter = 1 : nIter
        for n=1:N(i)
            if n==1
                pdf_arr(:,1) = pdf_init;
            else
                if y == 1
                    pdf_arr(:,n) = 2*p*pdf_arr(:,n-1);
                    pdf_arr(theta_axis > curr_med(n-1),n) = 2*(1-p)*pdf_arr(theta_axis > curr_med(n-1),n-1);
                else
                    pdf_arr(:,n) = 2*p*pdf_arr(:,n-1);
                    pdf_arr(theta_axis < curr_med(n-1),n) = 2*(1-p)*pdf_arr(theta_axis < curr_med(n-1),n-1);
                end
            end
            pdf_arr(:,n) = pdf_arr(:,n)/sum(pdf_arr(:,n));
            
            curr_med(n) = calc_med(pdf_arr(:,n),theta_axis);
            if theta0 > curr_med(n)
                x = 1;
            else
                x = 0;
            end
            
            noise(n) = (rand > (1-eps));
            y = mod(x + noise(n),2);
            
            if n == 1
                med_lower_bound(1) = 0.5;
                med_upper_bound(1) = 0.5;
            elseif n == 2
                if noise(1) == 0
                    med_lower_bound(2) = 0.5/(2*(1-p));
                    med_upper_bound(2) = 0.5/(2*(1-p));
                else
                    med_lower_bound(2) = 1 - 0.5/(2*(1-p));
                    med_upper_bound(2) = 1 - 0.5/(2*(1-p));
                end
            else
                if noise(n-1) == 0
                    med_lower_bound(n) = med_lower_bound(n-1)/(2*(1-p));
                    med_upper_bound(n) = med_upper_bound(n-1)/((2*(1-p))^(1-delta));
                else
                    med_lower_bound(n) = (med_lower_bound(n-1) + med_lower_bound(n-2))/2;
                    med_upper_bound(n) = med_upper_bound(n-2);
                end
            end
        end
    end
end
figure; hold all
plot(theta_axis,pdf_arr(:,1),'-k','LineWidth',2)
plot(theta_axis,pdf_arr(:,10),'-g','LineWidth',2)
plot(theta_axis,pdf_arr(:,20),'-c','LineWidth',2)
plot(theta_axis,pdf_arr(:,30),'-r','LineWidth',2)
plot(theta0*ones(1,100),linspace(dtheta,1,100),'--','LineWidth',2)
grid on; grid minor;
xlabel('theta'); ylabel('f_{n} (\theta | Y^{n})');
legend('n = 1','n = 10','n = 20','n = 30','\theta_{0}');
title(strcat('PM posteriors p = ',num2str(p),' \epsilon = ',num2str(eps)));
set(gca,'yscale','log')

figure;
subplot(211); hold all
plot(1:N,curr_med + 1e-10,'-o','LineWidth',1.5)
plot(1:N,med_lower_bound,'-*','LineWidth',1.5)
plot(1:N,med_upper_bound,'-*','LineWidth',1.5)
grid on; grid minor; legend('median','median lower bound','median upper bound');
set(gca,'yscale','log')

subplot(212);
plot(1:N,noise,'-o','LineWidth',1.5)
grid on; grid minor; legend('noise');

function med = calc_med(f,x)
f_fix = f/sum(f);
cum_f = cumsum(f_fix);
[~,med_idx] = min(abs(cum_f - 1/2));
med = x(med_idx);
end

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