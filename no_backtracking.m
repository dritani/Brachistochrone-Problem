function test6()
% VARIATION 1: "NO" BACKTRACKING
% NOTE: COMPARE CONVERGENCE INSTEAD OF SHAPE
clearvars
clc

fin = 10;
mesh = 100;
deltaX = fin/mesh;
[C,D,th]=centroid;

rho = 0.9;
c = 10^(-3);
alpha = 0.014;
limit = mesh+1;

X = 0:deltaX:fin;
Y = 10:(-deltaX):0;

Z = Y;
W = Y;
lim = 20000;
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








% big loop: real number of iterations
while(1)
    
    

    % convergence criterion? gradient G is what dimension
    m = m+1;
        Lra(m-1)=norm(Gra,2);
    if m==lim
        break
    end
        if m==50
        alpha=0.030;

    end
 
    if norm(Gra,2)<0.001
        m
        break
    end
    k = 1;
    if m==100
        L=Y;
    elseif m==250
        O=Y;
    elseif m==500
        P=Y;
    end
    Z = Y;
    % small loop: goes from 1 to 101, i.e. each point of the line
    while(1)
    k = k+1;
    if k==limit
        break
    end
    
    %NOTE: for backtracking, instead of making sure next step is smaller,
    %make sure next step is bigger due to convex nature (remember u flipped
    %the y values)
    dy = (Y(k+1)-Y(k-1))/(2*deltaX);
    ddy = (Y(k+1)-2*Y(k)+Y(k-1))/(deltaX^2);
    Gee = G(Y(k),dy,ddy);
    p_k = -Gee;
    Gra(k-1)=Gee;


    % Should p_k = -Gee?
    
% ask: do i need Z? as in, always use old values for new, or part old part
% new?
    
    Y(k) = Z(k) + alpha * p_k; 

    % wtf is the function itself, and how do i update it?
    % i.e. how to check backtracking for y(x+alpha*p_k) < y(x) +
    % c*alpha*...)?
    
% for quasi-newton: what dimension is Hk?
    

%Plot the log of the norm of the error between your solution and the exact solution
% so doesthis mean for each iteration, or only once it's done? 
% only final
%what parameters to vary?
% try adding constraints such as has to pass through a particular point in
% space
    end
    
   

Tra(m-1)=norm(Y(2:100)-fliplr(C(DnX)));
end

Y(2:100);
fliplr(C(DnX));
norm(Y(2:100)-fliplr(C(DnX)));

% different initial condition
% 
%hold on
%plot(-X+10,-Y+10,'c')
%plot(-X+10,-W+10,'g--')
%plot(-X+10,-L+10,'y.-')
%plot(-X+10,-O+10,'r')
%plot(-X+10,-P+10,'black')

%plot(D,-C+10)


%
% TOP KEK:::: 
%
% VARIATION OF PARAMETERS: DO THE FOLLOWING:

% 1."CHANGE" rho to 0.7:
% instead, don't reduce the alpha like usual

% 2. Different starting shape

% 3. Fixed midpoint.

%figure
%semilogy(1:m-2,Tra(1:m-2))

%figure
%semilogy(1:m-2,Lra(1:m-2))

Lra09 = dlmread('rho09.txt');

semilogy(1:3136,Lra09(1:3136),'g--','LineWidth',2)
hold on
grid on
semilogy(1:m-2,Lra(1:m-2),'r:','LineWidth',2)
xlabel('Iterations')
ylabel('Norm of Gradient')
legend('Rho = 0.9','Rho = 0.4')
title('Figure 7: Comparison of different Rho')
end

% what is C^2 for exact solution?
function [A,B,theta] = centroid()
    maxthet = 2.4120111439135253425264820657;
    theta = 0:0.0001:maxthet;
    %theta = 0:0.0239:maxthet+0.0239;
    c2=11.45834075;
    A = 0.5*c2*(1-cos(theta));
    B = 0.5*c2*(theta-sin(theta));
end

% L2 norm
function G1 = G(y,y_1,y_2) 

    % is this all we have to do to incorporate derivatives? nothing else?
    G1 = (-1)*(1+y_1^2+2*y*y_2)/(2*((y*(1+y_1^2))^(3/2)));
    
end




