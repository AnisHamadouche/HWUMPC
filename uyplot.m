function uyplot(y1,u1,dU1,N_sim)
    k=0:(N_sim-1);
    figure(1)
    subplot(2,1,1)
    plot(k,y1')
    xlabel('Sampling instant');
    ylabel('y');
    legend('Output');
    subplot(2,1,2)
    stairs(k,u1');
    xlabel('Sampling instant');
    ylabel('u');
    legend('control');
    figure(2)
    plot(k,dU1)
    ylabel('\Delta u')
    xlabel('Sampling Instant')
end