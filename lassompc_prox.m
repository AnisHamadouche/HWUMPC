%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%% 
% This is an exact Proximal-Gradient solver
function [u_prox,F_prox_subopt] = lassompc_prox(A, b, n, N, gamma,L,MAX_ITER,ABSTOL,lambda)
    % some constant parameters
    tk = 1/L;
    beta = 0.5;
    % cached computations for all methods
    AtA = A'*A;
    Atb = A'*b;
%       cvx: minimize(sum_square(A * u - b)+gamma*l1_norm(u))
    f = @(u) sumsqr(A*u-b);
    objective = @(u) sumsqr(A*u-b) + gamma*l1_norm(u);
    
    x0 = sprandn(n*N,1,2);
    x = x0;
    xprev = x;

    F_prox_subopt = zeros(MAX_ITER,1);
    for k = 1:MAX_ITER
        grad_x = AtA*x - Atb;
        while 1
            %z = prox_l1(x - tk*grad_x, tk*1/lambda);
            z = prox_l1(x - tk*grad_x, lambda);
            if f(z) <= f(x) + grad_x'*(z - x) + (1/(2*tk))*sumsqr(z - x)
                break;
            end
            tk = beta*tk;
        end
        xprev = x;
        x = z;
        F_prox_subopt(k) = objective(x);
        F_prox_subopt(1) = objective(x0);
        if k > 1 && abs(F_prox_subopt(k) - F_prox_subopt(k-1)) < ABSTOL
            break;
        end
    end
    u0 = x0;
    u_prox = reshape(x,n,N);
end