%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%Non-accelerated inexact PG error upper bounds estimation 
%based on Schmidt et al. 2010 and Hamadouche et al. 2020

x0 = u0; x_opt = u_prox(:);
epsilon0 = max(e2);

for k=1:MAX_ITER
    e1_norm(k) = norm(e1(:,k),2);
end

for k=1:size(e2)
    r_norm(k) = sqrt(2*e2(k)/L);
end

bound_prop0 = zeros(1,MAX_ITER);
bound_0 = zeros(1,MAX_ITER);
bound_prop1 = zeros(1,MAX_ITER);

j = 5; %probability parameter
bound_prop3 = zeros(1,MAX_ITER);

for k = 1:MAX_ITER
    
    bound_0(k) = (norm(x0-x_opt,2))^2*L/(2*k);

    Ak = vpa(sum(e1_norm(1:k))/L+sum(sqrt(e2(1:k)))/sqrt(L));
    Bk = vpa(sum(e2(1:k))/L);

%     bound_prop0(k) = (norm(x0-x_opt,2)+2*Ak+sqrt(2*Bk))^2*L/(2*k);
% 
%     bound_prop1(k) = (vpa(sum(e2(1:k)))+sum(e1_norm(1:k).*norm(x0-x_opt,2))+L*norm(x0-x_opt,2)^2/2)/k;
%     bound_prop3(k) = vpa(mean(e2)) + (j*(epsilon_0/2 + sqrt(n)*delta*norm(x0-x_opt,2))*sqrt(k)+L*norm(x0-x_opt,2)^2)/(2*k);
%     %    bound_aj4(k) = (mean(e2(1:k))+mean(e1_norm(1:k).*norm(x0-x_opt,2)))+(L*norm(x0-x_opt,2)^2/2)/k

    bound_prop0(k) = (sum(e2(1:k))+(sum(e1_norm(1:k)+L*r_norm(1:k)).*norm(x0-x_opt,2))+L*norm(x0-x_opt,2)^2/2)/k;
%    bound_prop1(k) = (vpa(sum(e2(1:k)))+sum(e1_norm(1:k).*norm(x0-x_opt,2))+L*norm(x0-x_opt,2)^2/2)/k;
%    bound_prop3(k) = vpa(mean(e2)) + (j*(epsilon0/2 + sqrt(n)*delta*norm(x0-x_opt,2))*sqrt(k)+L*norm(x0-x_opt,2)^2)/(2*k);
    bound_prop1(k) = (sum(e2(1:k))+sum(e1_norm(1:k)+L*r_norm(1:k)).*norm(x0-x_opt,2)+L*norm(x0-x_opt,2)^2/2)/k- L*sum(r_norm(1:k).^2)/2*k;
    bound_prop3(k) = mean(e2) + (j*(epsilon0/2 + (sqrt(n)*delta+sqrt(2*L*epsilon0))*norm(x0-x_opt,2))*2*sqrt(k)+L*norm(x0-x_opt,2)^2)/(2*k);
end

loglog(abs(F_prox_subopt-F_opt),'b','linewidth',2,'DisplayName','$F(\Delta U) - F(\Delta U^\star)$');
F_subopt = F_prox_subopt-F_opt;
hold
loglog(bound_prop0, 'linewidth',3,'Color',[0 0 0]+0.75,'DisplayName','Schmidt_1');
loglog(bound_prop1,':', 'linewidth',2,'Color',[0 0 0],'DisplayName','Thrm_1');
loglog(bound_prop3,'linewidth',2,'Color',[0 0 0],'DisplayName','Thrm_3');
%loglog(bound_0,'r','linewidth',0.5,'DisplayName','Error free');
loglog(abs(bound_prop0-bound_prop1),'r--','linewidth',1,'DisplayName','Imprvm_Thrm_1')
Imprvm_Thrm_1 = bound_prop0-bound_prop1;
loglog(abs(bound_prop0-bound_prop3),'g--','linewidth',1,'DisplayName','Imprvm_Thrm_3')
Imprvm_Thrm_3 = bound_prop0-bound_prop3;
ylabel('$F(\Delta U) - F(\Delta U^\star)$','interpreter','latex');
xlabel('$Iterations, k$','interpreter','latex');
xlim([0 k])

hl = legend('show');
set(hl, 'Interpreter','latex')