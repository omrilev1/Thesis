x_fullAccess = zeros(length(SNR),T);
x_zeroAccess = zeros(length(SNR),T);
x_partialAccess = zeros(length(SNR),T);
x_hat_partialAccess = zeros(length(SNR),T);
x_hat_partialAccess_feedback = zeros(length(SNR),T);

x_hat_estim_fullAccess = zeros(length(SNR),T); % x\hat(t|t)
x_hat_predict_fullAccess = zeros(length(SNR),T); % x\hat(t+1|t)

x_hat_estim_zeroAccess = zeros(length(SNR),T); % x\hat(t|t)
x_hat_predict_zeroAccess = zeros(length(SNR),T); % x\hat(t+1|t)

x_hat_estim_partialAccess = zeros(length(SNR),T); % x\hat(t|t)
x_hat_predict_partialAccess = zeros(length(SNR),T); % x\hat(t+1|t)

x_hat_estim_zeroAccess_modulo = zeros(length(SNR),T); % x\hat(t|t)
x_hat_predict_zeroAccess_modulo = zeros(length(SNR),T); % x\hat(t+1|t)

x_hat_estim_partialAccess_feedback = zeros(length(SNR),T); % x\hat(t|t)
x_hat_predict_partialAccess_feedback = zeros(length(SNR),T); % x\hat(t+1|t)

u_partialAccess = zeros(length(SNR),T);
u_fullAccess = zeros(length(SNR),T);
u_zeroAccess = zeros(length(SNR),T);
u_zeroAccess_modulo = zeros(length(SNR),T);
u_partialAccess_feedback = zeros(length(SNR),T);

J_partialAccess = zeros(length(SNR),T);
J_fullAccess = zeros(length(SNR),T);
J_zeroAccess = zeros(length(SNR),T);
J_zeroAccess_modulo = zeros(length(SNR),T);
J_partialAccess_feedback = zeros(length(SNR),T);


P_error_estim = zeros(length(SNR),T);
P_error_predict = zeros(length(SNR),T);

P_error_estim_zeroAccess = zeros(length(SNR),T);
P_error_predict_zeroAccess = zeros(length(SNR),T);

P_error_estim_partialAccess = zeros(length(SNR),T);
P_error_predict_partialAccess = zeros(length(SNR),T);

P_error_estim_zeroAccess_modulo = zeros(length(SNR),T);
P_error_predict_zeroAccess_modulo = zeros(length(SNR),T);

P_error_estim_partialAccess_feedback = zeros(length(SNR),T);
P_error_predict_partialAccess_feedback = zeros(length(SNR),T);

P_X_linear = zeros(1,T);
P_X_modulo = zeros(1,T);

P_X_hat_linear = zeros(1,T);
P_X_hat_modulo = zeros(1,T);