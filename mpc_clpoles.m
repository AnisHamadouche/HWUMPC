function [poles]=mpc_clpoles(num,denum,rw)
%num=[0.6];
%denum=[1 -0.8];
    sys=tf(num,denum);
    sysprime=tf([num zeros(1,order(sys))],flip(denum));
    f1=tf([1],[1 -1]);
    f2=tf([1 0],[-1 1]);
    sys_cl=sys*sysprime*f1*f2/rw+1;
    z=zero(sys_cl);
    poles=[];
    j=1;
    for i=1:size(z)
        if abs(z(i)) <= 1
            poles(j)=z(i);
            j=j+1;
        end
    end