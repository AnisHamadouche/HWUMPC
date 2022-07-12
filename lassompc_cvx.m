%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Authoer: Yun Wu

%% 
% This is a unconstraint solver using CVX
function u_cvx = lassompc_cvx(A, b, n, N, gamma, umax, umin)
    cvx_begin quiet
        variable u(n*N)
        minimize(sum_square(A * u - b)+gamma*l1_norm(u))
    cvx_end
    u_cvx = reshape(u,n,N);
end