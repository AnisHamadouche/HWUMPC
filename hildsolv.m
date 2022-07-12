function [x,lambda]=hildsolv(E,F,M,gamma_,MAX_ITER,QTOL)
    lambda=zeros(size(M,1),1);
    %QPhild() finds the global solution of quadratic optimization problem
    %using Hildreth's quadratic programming
    H=M*(E\M');
    K=gamma_+M*(E\F);
    %unconstrined solution
    x = -E\F;
    %check constraint violations
    if all(M*x <= gamma_)
        return;
    end
    %start iterative program  
    rho=10;
    for k=1:MAX_ITER
        %w(i,k) = -(H(i,1:i-1)*lambda(1:i-1,k)+H(i,(i+1):end)*lambda((i+1):end,k-1)+K(i))/H(i,i);
        lambda_p=lambda;
        n=size(lambda,1);
        for i=1:n
            w=H(i,:)*lambda-H(i,i)*lambda(i,1);
            w=w+K(i,1);
            w=-w/H(i,i);
            lambda(i,1)=max(0,w);
        end
        rho=(lambda-lambda_p)'*(lambda-lambda_p);
        if (rho<QTOL);break;end
    end
    x=-E\F-(E\M')*lambda;%lambda((i+1):end,k-1) is the correction term
end