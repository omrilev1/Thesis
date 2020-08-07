
%% Simulate RDF function for Gaussian source with causal side information and modulo encoder:
%% x ~ N(0,1)
%% y = x + n, n ~ N(0,\sigma_n^2)
%% w = [x + z]_{\Delta}, z ~ N(0,\sigma_z^20


% Parameters initialization
sigma_n = 1/4;
sigma_v = 1;
sigma_v_sigma_n = (sigma_n^2)*(sigma_v^2) / ((sigma_n^2) + (sigma_v^2)) ;

% baseline - two sided SI
d_twoSided = 0.001:0.001:sigma_v_sigma_n;
d_noSI = 0.01:0.01:sigma_v;

R_twoSided = 0.5 * log(sigma_v_sigma_n./d_twoSided);
R_twoSided(R_twoSided <= 0) = 0;
R_noSI = 0.5 * log(sigma_v_sigma_n./d_noSI);
R_noSI(R_noSI <= 0) = 0;

R_causalSI = zeros(size(R_twoSided));

dd = 1e-4;
dz = 1e-2;
eps = 1e-3;

Delta = dd:dd:exp(max(R_causalSI));
H_Delta = zeros(size(Delta));
for i=1:length(Delta)
    H_Delta(i) = calcGaussianModuloEntropy(Delta(i));
end

for i=1:length(R_causalSI)
    
    % Start with find the possible values of \sigma_z and \Delta
    currR = R_twoSided(i);
    sigma_z_2 = dz : dz : 1/(exp(2*currR) - 1);
    
    % Rate table for every combination of sigma_z and Delta
    H1 = zeros(length(sigma_z_2),length(Delta));
    H2 = zeros(length(sigma_z_2),length(Delta));
    
    for k=1:length(sigma_z_2)
        
        [~,idx1] = min(abs(Delta - Delta/sqrt(sigma_z_2(k) + sigma_v^2)));
        [~,idx2] = min(abs(Delta - Delta/sqrt(sigma_z_2(k))));
        
        H1(k,idx1) = H_Delta(idx1);
        H1(k,idx2) = H_Delta(idx2);
        
        H1(Delta/sqrt(sigma_z_2(k) + sigma_v^2) > max(Delta)) = 0.5*log(2*pi*(sigma_z_2(k) + sigma_v^2));
        H1(Delta/sqrt(sigma_z_2(k) + sigma_v^2) < min(Delta)) = log(Delta(Delta/sqrt(sigma_z_2(k) + sigma_v^2) < min(Delta))/sqrt(sigma_z_2(k) + sigma_v^2));
        H2(Delta/sqrt(sigma_z_2(k)) > max(Delta)) = 0.5*log(2*pi*(sigma_z_2(k)));
        H2(Delta/sqrt(sigma_z_2(k)) < min(Delta)) = log(Delta(Delta/sqrt(sigma_z_2(k)) < min(Delta))/sqrt(sigma_z_2(k)));
        
        curr_MI = H1(k,:) - H2(k,:);
        [~,validIdx] = find(abs(curr_MI - currR) < eps);
        
    end
    
    I = H1 - H2;
    
    
end

function [H] = calcGaussianModuloEntropy(Delta)
%% Calculate the entropy of the RV [x]_\Delta, X~N(0,1)

maxVal = 6;
N = 1e2;
dx = Delta/N;
x = -1*ceil(maxVal/Delta)*Delta - Delta/2 : dx : ceil(maxVal/Delta)*Delta + Delta/2 - dx;
f = (1/sqrt(2*pi))*exp(-x.^2/2);

f_mod = dx * sum(reshape(f,N,[]),2);

H = -1*dx*sum(f_mod(f_mod > 0).*log(f_mod(f_mod > 0)));



end