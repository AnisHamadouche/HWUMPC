%%Accelerated PG solver for time-critical application with MPC length N = 2 and 30 max
%%iterations

delta = 2.2E-1; %MAX tolerable rounding error (in gradient)
epsilon_0 = 1E-4; %MAX tolerable prox error
N = 2; %MPC horizon length
ABSTOL = eps; %Tolerance
MAX_ITER = 30; % gradient descent maximum iteration number
digits(64);% for vpa(.)

[H, F, x0, m, n, l, xmax, xmin, ymax, ymin] = model_spacecraft(N);

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

objective = @(u) sumsqr(A_hat*u-b_hat) + gamma*l1_norm(u);

% Run PG
% [u_prox,F_prox_subopt,e1,e2,u0,k]=lassompc_iaccprox(A_hat, b_hat, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon_0);
% F_opt = objective(u_prox(:));

%Run Maximum Hands-off Control
[u_prox,F_prox_subopt,e1,e2,u0,k]=lassompc_iaccprox(A_hat, b_hat, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon_0);
F_opt = objective(u_prox(:));

dU = reshape(repmat(u_prox,300,1),size(u_prox,1),300*size(u_prox,2));
%Plot upper bounds
acc_bounds_v1_2