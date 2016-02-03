function test8()
% DFP/BFGS Quasi-Newton Final Code
clearvars
clc

fin = 10;
mesh = 100;
deltaX = fin/mesh;
[C,D,th]=centroid;

rho = 0.9;
c = 10^(-4);
alpha = 0.00005;
limit = 101;

X = 0:deltaX:fin;
Y = fin:(-deltaX):0;
Z = Y;
W = Y;
lim = 120000;
m = 1;
Gra = ones(1,limit-2);
Lra = ones(lim-1,1);
Tra = ones(lim-1,1);
DnX=ones(limit-2,3);
u=0;
while(1)
    u=u+1;

    if u==limit-1
        break
    end
        now=X(u+1);
[row,column] = find(D>now-0.001 & D<now+0.001);
    DnX(u)=column(1);

end
DnX=DnX(:,1,1);

H_old = Y;
H_new = Y;
for u=1:size(Y,2)
    H_old(u) = 1;
    H_new(u) = 1;
end

% big loop: real number of iterations
while(1)
    
    m = m+1;
    Lra(m-1)=norm(Gra,2);
    if m==lim
        m
        break
    end
        if norm(Gra,2)<0.001
        m
        break
    end
    k = 1;
        if m==5000
        L=Y;
    elseif m==10000
        O=Y;
    elseif m==20000
        P=Y;
    elseif m==2000
        S=Y;
    end
    H_old=H_new;
    Z=Y;
    % small loop: goes from 1 to 100, i.e. each point of the line
    while(1)
    k = k+1;
    if k==limit
        break
    end
    
    dy = (Y(k+1)-Y(k-1))/(2*deltaX);
    ddy = (Y(k+1)-2*Y(k)+Y(k-1))/(deltaX^2);
    Gee = G(Y(k),dy,ddy);
    p_k = -H_old(k)*Gee;
    Gra(k-1)=Gee;
  
    deltaD = alpha*p_k;
        Y(k) = Z(k) + deltaD; 
    
    deltaG = G(Y(k),dy,ddy) - G(Z(k),dy,ddy);
 
    % DFP
    H_new(k) = H_old(k) + ((deltaD * deltaD)/(deltaD * deltaG)) - ...
       (H_old(k) * deltaG * deltaG * H_old(k) )/(deltaG * H_old(k) * deltaG);

    %  BFGS
    %H_new(k) = H_old(k) + (1 + (deltaG * H_old(k) * deltaG)/(deltaG*deltaD))*...
    %  ((deltaD*deltaD)/(deltaD*deltaG)) - (H_old(k)*deltaG*deltaD+H_old(k)*...
    % deltaG*deltaD)/(deltaG*deltaD);
   
    % Rank 1 formula:
    H_new(k)=H_old(k)+(deltaD-H_old(k)*deltaG)/deltaG;
    

    % Should p_k = -Gee?
    
end
Tra(m-1)=norm(Y(2:100)-fliplr(C(DnX)));
end

hold on
grid on
plot(-X+10,-Y+10,'c')
plot(-X+10,-W+10,'g--')
plot(-X+10,-S+10,'m')
plot(-X+10,-L+10,'y.-')
plot(-X+10,-O+10,'r')
plot(-X+10,-P+10,'black')

[C,D]=centroid;

plot(D,-C+10)
xlabel('X')
ylabel('Y')
title('Figure 4: Quasi-Newton Solution')
legend('Final value','Initial value','2000 iterations','5000 iterations','10000 iterations','20000 iterations','Cycloid')

figure
semilogy(1:m-2,Lra(1:m-2))
grid on
xlabel('Iterations')
ylabel('Norm of error')
title('Figure 5: Quasi-Newton Gradient Convergence')

figure
semilogy(1:m-2,Tra(1:m-2))
grid on
xlabel('Iterations')
ylabel('Norm of gradient')
title('Figure 6: Quasi-Newton Error Convergence')

end

function [A,B,theta] = centroid()
    maxthet = 2.4120111439135253425264820657;
    theta = 0:0.0001:maxthet;
    %theta = 0:0.0239:maxthet+0.0239;
    c2=11.45834075;
    A = 0.5*c2*(1-cos(theta));
    B = 0.5*c2*(theta-sin(theta));
end

function G1 = G(y,y_1,y_2) 

    G1 = (-1)*(1+y_1^2+2*y*y_2)/(2*((y*(1+y_1^2))^(3/2)));
    
end




