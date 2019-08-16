% This script calculates and plot the exponent of Horstein scheme
p = 0:0.01:0.25; % BSC(0.11) has capacity of 1/2 [bit/channel use]
eps_perc = 1;
lambda = 0:0.001:1;
beta = 0 : 0.001 : 2;
opt_beta = zeros(size(p));
opt_beta_tilde = zeros(size(p));
lambda_max = 1;
N = 10;
figure; hold all
for i=1:length(eps_perc)
    for j=1:length(p)
        eps = p(j)*eps_perc(i);
        opt_Psi = zeros(size(beta));
        opt_Psi_Tilde = zeros(size(beta));
        
        for k=1:length(beta)
            
            curr_beta = beta(k);
            curr_p = p(j);
            
            curr_psi = @(lambda) lambda*curr_beta + log2(((2*curr_p)^(lambda) + (2*(1-curr_p))^(lambda))/2);
            lambda_opt = fminbnd(curr_psi,0,lambda_max);
            opt_Psi(k) = -log2(((2*curr_p).^(lambda_opt) + (2*(1-curr_p)).^(lambda_opt))/2) - lambda_opt*curr_beta;
            
            curr_psi_tilde = @(lambda)lambda*curr_beta + log2(max(eps*(2*curr_p)^(lambda) + (1-eps)*(2*(1-curr_p))^(lambda),...
                (1-eps)*(2*curr_p)^(lambda) + eps*(2*(1-curr_p))^(lambda)));
            lambda_opt_tilde = fminbnd(curr_psi_tilde,0,lambda_max);
            opt_Psi_Tilde(k) = -log2(max(eps.*(2*curr_p).^(lambda_opt_tilde) + (1-eps).*(2*(1-curr_p)).^(lambda_opt_tilde),...
                (1-eps).*(2*curr_p).^(lambda_opt_tilde) + eps.*(2*(1-curr_p)).^(lambda_opt_tilde))) - lambda_opt_tilde*curr_beta;
        end
        
        [~,beta_opt_idx] = min(abs(beta - opt_Psi + 1/N));
        curr_beta_opt = beta(beta_opt_idx);
        
        [~,beta_opt_tilde_idx] = min(abs(beta - opt_Psi_Tilde + 1/N));
        curr_beta_opt_tilde = beta(beta_opt_tilde_idx);
        
        opt_beta(j) = curr_beta_opt;
        opt_beta_tilde(j) = curr_beta_opt_tilde;
    end
    
    plot(p,opt_beta,'-*','LineWidth',2);
    plot(p,opt_beta_tilde,'-o','LineWidth',2);
end
grid on; grid minor;
xlabel('p'); ylabel(strcat('\beta : Exponent for n = ',num2str(N)));
legend('\eps = 0','\eps = 0.25','\eps = 0.5','\eps = 0.75','\eps = 1');

epsilon = 0:0.01:0.5;
p = 0:0.01:0.5;
lambda = 0:0.01:1;
eps = 0.1;
% [eps_mesh,p_mesh] = meshgrid(epsilon,p);
[lambda_mesh,p_mesh] = meshgrid(lambda,p);

f1 = eps.*(2*p_mesh).^(lambda_mesh) + (1 - eps).*(2*(1 - p_mesh)).^(lambda_mesh);
f2 = (1 - eps).*(2*p_mesh).^(lambda_mesh) + eps.*(2*(1 - p_mesh)).^(lambda_mesh);
f3 = ((2*p_mesh).^(lambda_mesh) + (2*(1 - p_mesh)).^(lambda_mesh))/2;
f_tot = max(f1,f2);
figure; surf(lambda_mesh,p_mesh,-1*log2(f_tot));colorbar
xlabel('\lambda'); ylabel('p');title('f optimized');
figure; surf(lambda_mesh,p_mesh,-1*log2(f3));colorbar
xlabel('\lambda'); ylabel('p');title('f Anusha');

function [psi,p_grid,lambda_grid] = calcPsi_2D(p,lambda)
% This function calculate psi(p) = -log(((2*p)^lambda +
% (2*(1-p))^lambda)/2)
[p_grid,lambda_grid] = meshgrid(p,lambda);
psi = -log2(((2*p_grid).^(lambda_grid) + (2*(1-p_grid)).^(lambda_grid))/2);
end

function [psi] = calcPsi_1D(p,lambda)
% This function calculate psi(p) = -log(((2*p)^lambda +
% (2*(1-p))^lambda)/2)
psi = -log2(((2*p).^(lambda) + (2*(1-p)).^(lambda))/2);
end

function [psi,p_grid,lambda_grid,eps_grid] = calcPsiTilde_2D(p,lambda,eps)
% This function calculate
% psi(p) = -log((max{eps*(2*p)^lambda + (1-eps)*(2*(1-p))^lambda,(1-eps)*(2*p)^lambda + eps*(2*(1-p))^lambda}))
[p_grid,lambda_grid,eps_grid] = meshgrid(p,lambda,eps);
psi = -log2(max(eps_grid.*(2*p_grid).^(lambda_grid) + (1-eps_grid).*(2*(1-p_grid)).^(lambda_grid),...
    (1-eps_grid).*(2*p_grid).^(lambda_grid) + eps_grid.*(2*(1-p_grid)).^(lambda_grid)));
end

function [psi] = calcPsiTilde_1D(p,lambda,eps)
% This function calculate
% psi(p) = -log((max{eps*(2*p)^lambda + (1-eps)*(2*(1-p))^lambda,(1-eps)*(2*p)^lambda + eps*(2*(1-p))^lambda}))
psi = -log2(max(eps.*(2*p).^(lambda) + (1-eps).*(2*(1-p)).^(lambda),...
    (1-eps).*(2*p).^(lambda) + eps.*(2*(1-p)).^(lambda)));
end