%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%Non-accelerated inexact PG error upper bounds estimation 
%based on Schmidt et al. 2010 and Hamadouche et al. 2020
MAX_ITER = k;
x0 = u0; x_opt = u_prox(:);
epsilon0 = max(e2);

for k=1:MAX_ITER
    e1_norm(k) = norm(e1(:,k),2);
end

for k=1:size(e2)
    r_norm(k) = sqrt(2*e2(k)/L);
end

bound_0 = zeros(1,MAX_ITER);
bound_prop0 = zeros(1,MAX_ITER);
bound_prop4 = zeros(1,MAX_ITER);
j = 5;
% for j = 1:0.5:100
bound_prop5 = zeros(1,MAX_ITER);
%bound_aj4 = zeros(100,1)

for k = 1:MAX_ITER
    
    i = k-1;
    alpha(k)=(i+2)/2;
    bound_0(k) = (norm(x0-x_opt,2))^2*L/2/alpha(k)^2;
    Ak = sum((1:k).*e1_norm(1:k))/L+sum((1:k).*sqrt(e2(1:k)'))/sqrt(L);
    Bk = sum((1:k).^2.*e2(1:k)')/L;

    bound_prop0(k) = (norm(x0-x_opt,2)+2*Ak+sqrt(2*Bk))^2*L/2/alpha(k)^2;

    bound_prop4(k) = ((sum(alpha(1:k).^2.*e2(1:k)')+sum(alpha(1:k).*(e1_norm(1:k)+L*r_norm(1:k)))*norm(x0-x_opt,2))+L*norm(x0-x_opt,2)^2/2)/alpha(k)^2;
    bound_prop5(k) = ((mean(e2)*sum(alpha(1:k).^2)) + j*sqrt(n*i*(i+1)*(2*i+1)/6)*delta*norm(x0-x_opt,2)+L*norm(x0-x_opt,2)^2/2+j*sqrt(2*epsilon0*i*(i+1)*(2*i+1)/L/6)*delta*norm(x0-x_opt,2))/alpha(k)^2;
    %    bound_aj4(k) = (mean(e2(1:k))+mean(e1_norm(1:k).*norm(x0-x_opt,2)))+(L*norm(x0-x_opt,2)^2/2)/k
end

loglog(abs(F_prox_subopt-F_opt),'b','linewidth',2,'DisplayName','$F(\Delta U) - F(\Delta U^\star)$');
F_subopt=F_prox_subopt-F_opt;
hold
% loglog(bound_0,'y:', 'linewidth',3,'DisplayName','Optimal');
loglog(bound_prop0, 'linewidth',3,'Color',[0 0 0]+0.75,'DisplayName','Schmidt_2');
loglog(bound_prop4,':', 'linewidth',2,'Color',[0 0 0],'DisplayName','Thrm_4');
loglog(bound_prop5,'linewidth',2,'Color',[0 0 0],'DisplayName','Thrm_5');
%loglog(bound_0,'r','linewidth',0.5,'DisplayName','Error free');
loglog(abs(bound_prop0-bound_prop4),'r--','linewidth',1,'DisplayName','Imprvm_Thrm_4')
loglog(abs(bound_prop0-bound_prop5),'g--','linewidth',1,'DisplayName','Imprvm_Thrm_5')

ylabel('$F(\Delta U) - F(\Delta U^\star)$','interpreter','latex');
xlabel('$Iterations, k$','interpreter','latex');
xlim([0 k])

hl = legend('show');
set(hl, 'Interpreter','latex')