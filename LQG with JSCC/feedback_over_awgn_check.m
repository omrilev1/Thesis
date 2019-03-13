
% channel with feedback try


snr = 5:5:20; % snr to check
P = 1; % input power constraint
N_feedback = 4:1:10;
delta_snr = 15;

dist = zeros(length(snr),length(N_feedback));
N = 3e3;
snrLin = 10.^(snr/10);
PAM_order = 2.^(1:6);
curr_PAM = zeros(size(snr));

d = sqrt(12*P);
for i=1:length(snr)
    
    for j=1:length(N_feedback)
        
        curr_capacity = 0.5*N_feedback(j)*log(1+snrLin(i));
        % calculate closest PAM order
        [~,curr_PAM_idx] = min(abs(2^(floor(log2(curr_capacity))) - log2(PAM_order)));
        curr_PAM(i) = PAM_order(curr_PAM_idx);
        curr_dist = 0;
        
        for n = 1:N
            theta = d*(rand - 0.5);
            [theta_final] = feedback_over_awgn(theta,P,curr_PAM(i),snr(i),delta_snr,N_feedback(j));
            curr_dist = curr_dist + abs(theta_final-theta).^2;
        end
        dist(i,j) = curr_dist/N;
    end
end

figure; hold all
currLegend = [];
LineStyle = {'-*k','-pg','--c','-or'};
for i=1:length(snr)
    plot(N_feedback,10*log10(dist(i,:)),char(LineStyle(i)),'LineWidth',2.5)
    curString = strcat('MSE For snr = ',num2str(snr(i)),' [dB], PAM order = ',num2str(curr_PAM(i)));
end
grid on; grid minor;
set(gca,'yscale','log')
title('MSE for different SNRs and feedback iterations')
legend(currLegend)
xlabel('N feedback iterations')
ylabel('MSE [dB]')

