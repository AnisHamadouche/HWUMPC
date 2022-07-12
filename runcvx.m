%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Authoer: Yun Wu

%% 
x_cvx = x0;
G = F * x_cvx;
%A_hat = (1/sqrt(2)) * H^(1/2);
A_hat = (1/sqrt(2)) * sqrtm(H);
%b_hat = (-1/sqrt(2)) * (H^(-1/2) * G);
b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
L = max(eigs(A_hat'*A_hat));
lambda = 1/L;
gamma_max = norm(A_hat' * b_hat,'inf');
gamma = 0.01*gamma_max;
tic;
u_cvx = lassompc_cvx(A_hat, b_hat, n, N, gamma);
time_cvx=toc;
fprintf('CVX solution: \n');
for i=1:length(u_cvx)
fprintf(' %f \t', u_cvx(i));
end
fprintf('\n');
fprintf('\n');