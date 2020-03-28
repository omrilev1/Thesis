

snr = 10:1:50;
snrLin = 10.^(snr/10);

MSE_spiral = zeros(size(snr));
bestDelta = zeros(size(snr));

% Spiral parameters 
sigma_x = 0.25; % for the spiral calculations, later we'll compenste over this
eta = 0.16;  % curve "speed"
Delta = 1e-3:1e-3:1.75; % Deltas to optimize over 

for i=1:length(snr)
    
        % Find best Spiral SDR for that case 
    alpha = eta*sqrt(2*pi^5) ./ (Delta*sigma_x*(1 - exp(-1/(2*sigma_x^2))));
    P_r = (1 - erf(Delta/(2*sqrt(2)/sqrt(snrLin(i)))));
    eps_th = (1./(sqrt(pi)*(alpha.^2))) .* P_r .* (erf(sqrt(2)/(2*sigma_x))*(4*sqrt(pi)*sigma_x^2*(alpha.^2) + eta^2*(Delta.^2)*(pi^(4.5))) - ...
        4*sqrt(2)*sigma_x*alpha*exp(-1/(2*sigma_x^2)).*(alpha + 2*eta*Delta*(pi^2)) + 8*sqrt(2)*eta*(pi^2)*sigma_x*Delta.*alpha);
    [MSE_spiral(i),minIdx] = min(eps_th + (1/snrLin(i))./(alpha.^2));
    bestDelta(i) = Delta(minIdx); 
    
end

SDR_Spiral = sigma_x^2 ./ MSE_spiral;
