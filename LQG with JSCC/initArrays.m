x_fullAccess = zeros(length(SNR),T,N_avg);
x_zeroAccess = zeros(length(SNR),T,N_avg);
x_partialAccess = zeros(length(SNR),T,N_avg);
x_hat_partialAccess = zeros(length(SNR),T,N_avg);
x_hat_partialAccess_feedback = zeros(length(SNR),T,N_avg);
x_hat_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg);

x_hat_estim_fullAccess = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_fullAccess = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

x_hat_estim_zeroAccess = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_zeroAccess = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

x_hat_estim_partialAccess = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_partialAccess = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

x_hat_estim_zeroAccess_modulo = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_zeroAccess_modulo = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

x_hat_estim_partialAccess_feedback = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_partialAccess_feedback = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

x_hat_estim_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg); % x\hat(t|t)
x_hat_predict_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg); % x\hat(t+1|t)

u_partialAccess = zeros(length(SNR),T,N_avg);
u_fullAccess = zeros(length(SNR),T,N_avg);
u_zeroAccess = zeros(length(SNR),T,N_avg);
u_zeroAccess_modulo = zeros(length(SNR),T,N_avg);
u_partialAccess_feedback = zeros(length(SNR),T,N_avg);
u_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg);

J_partialAccess = zeros(length(SNR),T,N_avg);
J_fullAccess = zeros(length(SNR),T,N_avg);
J_zeroAccess = zeros(length(SNR),T,N_avg);
J_zeroAccess_modulo = zeros(length(SNR),T,N_avg);
J_partialAccess_feedback = zeros(length(SNR),T,N_avg);
J_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg);


P_error_estim = zeros(length(SNR),T,N_avg);
P_error_predict = zeros(length(SNR),T,N_avg);

P_error_estim_zeroAccess = zeros(length(SNR),T,N_avg);
P_error_predict_zeroAccess = zeros(length(SNR),T,N_avg);

P_error_estim_partialAccess = zeros(length(SNR),T,N_avg);
P_error_predict_partialAccess = zeros(length(SNR),T,N_avg);

P_error_estim_zeroAccess_modulo = zeros(length(SNR),T,N_avg);
P_error_predict_zeroAccess_modulo = zeros(length(SNR),T,N_avg);

P_error_estim_partialAccess_feedback = zeros(length(SNR),T,N_avg);
P_error_predict_partialAccess_feedback = zeros(length(SNR),T,N_avg);

P_error_estim_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg);
P_error_predict_partialAccess_feedback_linear = zeros(length(SNR),T,N_avg);




power_fullAccess = 0;
power_zeroAccess = 0;
power_partialAccess = 0;
power_zeroAccess_modulo = 0;
power_partialAccess_feedback = 0;