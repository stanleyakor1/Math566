clc;clear all;close all;
format long;
%% 
% *QUESTION 1*

m=50; n=12;
t=linspace(0,1,m)';
A=t.^(0:n-1);
f=cos(4*t);

%%% 1a %%%%%

c1=A\f;

%%% 1b %%%%%

[Q1,R1]=GS(A);
c2=R1\(Q1'*f);

%%% 1c %%%%%

[Q2,R2]=MGS(A);
c3=R2\(Q2'*f);

%%% 1d %%%%%
[V3,R3]=house(A);
Q3=house2q(V3);
c4=R3\(Q3'*f);


%%% 1e %%%%%
[Q4,R4]=qr(A);
c5=R4\(Q4'*f);


%%%% 1f %%%%%%
[U,S,V] = svd(A,0);
S(S>0)=1./S(S>0);
c6 = (V*S*U')*f;

%%
T= table(c1,c2,c3,c4,c5,c6,'VariableNames',{'Normal eqns mthd','GS','MGS','Householder','M-qr','svd'})
%%
writetable(T,"table.csv")
%% 
% *observation*
% 
% The figure above shows the values that appear wrong or far from the normal 
% equation method. The values highlighted in red correspond to those values whose 
% difference with the normal equations method is greater than  ${10}^{-5}$. 

find(abs((c1-c2))>1e-4) % locate index where the absolute difference between the normal eqns  and the method is greater than 10e-5
find(abs((c1-c3))>1e-4)
%% 
% 
% 
% 
%% 
% From the Figure above, we observe that the GS method is least stable due to 
% roundoff errors and hence has more coefficients far from the ones produced by 
% the normal equations method. Next is the MGS method which is more stable than 
% the GS. Considering the threshold of ${10}^{-5}$ used in this analysis, the 
% other methods produce coefficients which are close to those of the normal equations.

abs(c1-c5) % Absolute difference between the normal equation method and the qr method
%% 
% It can be observed that the absolute difference between the qr and the normal 
% equations method is very small. 

figure(1)
res10=f-A*c1;
plot(res10,'--om')
xlim([1,size(f,1)])
legend('Matlab routine')
xlabel('t')
ylabel('Residuals')
axis tight
%%
figure(2)
res1=f-A*c5;
plot(res1,'--or')
xlim([1,size(f,1)])
xlabel('t')
ylabel('Residuals')
legend('qr')
axis tight
%% 
% *observation*
% 
% We observe that the maximum residual in both methods is approximately ${3e}^{-09}$.
%% 
% *Question 2*

%%%%%%%%%%%% 2a %%%%%%%%%%%%%%%%%%%%%%
%% 
% ${\|r\|}_w^2$ = ${\left(\textrm{Ax}-b\right)}^T W\left(\textrm{Ax}-b\right)$
% 
% = ($x^T$$A^T$$-$ $b^T$)$W$$\left(\textrm{Ax}-b\right)$
% 
% = $x^T A^T \textrm{WAx}$$-$$x^T A^T \textrm{Wb}$ $-$$b^T \textrm{WAx}$ $+$ 
% $b^T \textrm{Wb}$
% 
% Taking the derivative with respect to x and equating to zero, we have;
% 
% $$\frac{d(||r||^2_w)}{dx}= 2A^TWA -A^TWb - b^TWA = 0$$
% 
% But, $A^T \textrm{Wb}$ = $b^T \textrm{WA}$
% 
% ${2A}^T \textrm{WAx}$ = $2A^T \textrm{Wb}$
% 
% x = ${\left(A^T \textrm{WA}\right)}^{-1}$$A^T \textrm{Wb}$

%%%%%%%%%%%% 2b %%%%%%%%%%%%%%%%%%%%%%
%%
W = zeros(m);
delta = 0.8;
for i = 1:m
    W(i,i) = exp(-(abs(1/23 - t(i))/delta)^2);
end
%%
x_w = (A'*W*A)\(A'*W*f)
%%
x_w(11)
%%
c1 - x_w
%%

figure(11)
res21=f-A*x_w;
plot(res21,'--om'), hold on
plot(res10,'--k')
xlim([1,size(f,1)])
legend('weighted least squares','least squares')
xlabel('t')
ylabel('Residuals')
axis tight
%% 
% From the Figure above, we observe that the weighted least squares out-perform 
% the least squares for t < 30. 
%% 
% *Question 3*
%% 
% 
% 
% $$\kappa = \frac{||J(X)||}{\frac{||f(x)||}{||x||}}$$
% 
% $$f(x) = \frac{log(x+1)}{x}$$
% 
% $$f^\prime(x) = \frac{x - (x+1)log(x+1)}{x^2(x+1)}$$
% 
% $$\kappa = \frac{|x -(x+1)log(x+1)|}{|(x+1)x^2|} \times \frac{|x^2|}{|log(x+1)|}$$
% 
% $$\kappa= \frac{|x - (x+1)log(x+1)|}{|(x+1)log(x+1)|}$$
% 
% 
% 
% 

%%%%%%%%%%%%% 3a  %%%%%%%%%%%%%%%%%%%%
kappa = @(x)abs(x-log(x+1)*(x+1))/abs(log(x+1)*(x+1));
x=-1e-7:1e-8:1e-7;
xk=zeros(size(x));
for i=1:length(xk)
    xk(i)= kappa(x(i));
end
xk
%% 
% *Observation*
%% 
% It can be obsereved that the condition number of the points near zero are 
% very small in the order of $10^{-7}$. Thus we can conclude that when the input 
% is changed by a small amount, the corresponding output will be changed by a 
% small amount.

%%%%%%%%%%% 3b %%%%%%%%%%%%%%%%%
f=@(x)(log(x+1)./x);
j=0:1:520;
xj=2.^(-52+j./10);

figure(5)
semilogx(xj,f(xj))
xlabel('x')
ylabel('f(x)')
legend('f(x)')
%% 
% *Observation*
%% 
% $f(x)$ appears to be unstable near x = 0.

%%%%%%%%% 3c %%%%%%%%%%%%%%%
z =1+xj;
y = @(z)(log(z)./(z-1));

figure(6)

semilogx(xj,f(xj),'r',linewidth=2), hold on
semilogx(xj,y(z),'k',linewidth=2)
hold off;
xlabel('x')
ylabel('f(x)')
legend('f(x)','y(x)');
%% 
% *Observation*
% 
% $y(x)$appear to be stable near zero while $f(x)$ is unstable.
%% 
% *Question 4*

%%%%%% 4a %%%%%%%%%%%%
%% 
% Let us consider the equation, $\frac{||\delta x||}{||x||} \leq\kappa \frac{||r||}{||b||}$ 
% 
% Looking at the right hand side, we observe that $\kappa$ is the dominant factor 
% in determining how the output will change with respect to a perturbation in 
% the input. Thus we conclude that the size of the residual does not give us enough 
% information about how close the solution is to the exact but the condition number  
% $\kappa$  does.

%%%%%% 4b %%%%%%%%%%%%

A = 1./hankel(2:6,6:10);
b = [0.882 0.744 0.618 0.521 0.447]';
x = A\b
%%
%%%%%% 4c %%%%%%%%%%%%
%% 
% $A\hat{x}$= b + $\delta b$
% 
% $\hat{x}$ = $A^{-1} b$ + $A^{-1}$$\delta b$
% 
% $\hat{x}$ = x + $A^{-1}$$\delta b$
% 
% $\hat{x}$ = x + $\delta x$
% 
% $\delta x$ = $\hat{x}$  - x


kapp = cond(A)

% small pertubation of b

db = 1e-4:2.5e-5:2*1e-4;

x_cap = x + A\db';

dx =x_cap - x;

% Upper bound
kappa_db_b = kapp*norm(db)/norm(b)
%%
dx_d=norm(dx,'fro')/norm(x,'fro')
%% 
% *Observation*
%% 
% First we observe that the condition number ($\kappa$) is really large and 
% the product of the relative error in b and $\kappa$ is a lot greater than the 
% relative error in x. This implies that a small perturbation in the input will 
% lead to a large perturbation in the output.
%% 
% *QUESTION 5*

n=30;
kapp = zeros(n,1)';
for i = 1:length(kapp)
    kapp(i) = cond(vandermonde(i,i));
end
%%

figure(7)
semilogy(1:1:n,kapp,linewidth=2)
xlim([1,30])
xlabel('n')
ylabel('cond(A)')
%% 
% *Observation*
% 
% The condition number appears to grow linearly from n = 1 until n = 23 where 
% it begins to fluctuate. This could be attributed to the intrinsic nature of 
% the matrix which is ill-conditioned.

n=30;
kapp2 = zeros(n,1)';
for i = 1:length(kapp)
    kapp2(i) = cond(vandermonde(2*i-1,i));
end
%%
figure(8)
semilogy(1:1:n,kapp2,linewidth=2)
xlim([1,30])
xlabel('n')
ylabel('cond(A)')
%% 
% *Observation*
%% 
% The condition number appears to grow linearly from n = 1 until n = 25 where 
% it reaches steady state. 
%% 
% *QUESTION 6*

x=1.920:0.001:2.080;
r=zeros(size(x));
for i=1:length(x)
    r(i)=p(x(i));
end
figure(9)
plot(x,r,'r',LineWidth=2)
xlabel('x')
ylabel('p(x)')
legend('Explicit')
axis tight
%%
r_q=(x-2).^9;
figure(10)
plot(x,r_q,'k',LineWidth=2)
xlabel('x')
ylabel('p(x)')
legend('p(x) = (x-2)^9')
axis tight
%% 
% *Observation*
%% 
% We observe that the figure produced using $p(x) = (x-2)^9$ is more smooth 
% compared to the figure obtained from the explicit expression of $p(x)$. This 
% is due to the accumulation of roundoff error in the non-explicit case.
%% 
% *QUESTION 7*
%% 
% P is a permutation matrix, it is orthogonal and symmetric. Thus $P^2$ = $P$.

P=zeros(10);
P(1,end) = 1; P(2,5) = 1; P(3,8) = 1; P(4,1) = 1;
P(5,4) = 1; P(6,9) = 1; P(7,2) = 1; P(8,3) = 1; P(9,6) = 1; P(10,7)= 1;
%%
skeel_K = norm(P)
Kappa_P = cond(P)
%%
P(:,3)=10e-11*P(:,3);
skeel_K = norm(abs(inv(P))*abs(P))
%%
kappa_P = cond(P)
%% 
% 
% 
% The skeel condition number is always less than or equal to $\kappa$ and it 
% is invariant of the row or column scaling.
%% 
% 

function  r = p(x)
    r = x^9 -18*x^8 + 144*x^7 - 672*x^6 + 2016*x^5 -4032*x^4 + 5376*x^3 - 4608*x^2 + 2304*x -512;
end

function A = vandermonde(m,n)
 t=linspace(0,1,m)';
 A=t.^(0:n-1);
end

function [V,R]=house(A)
[m, n]=size(A);
assert(m>=n,'The row count must be greater than or equal to the column count')
R=A;
V=zeros(m,n);
for k=1:n
    x=R(k:end,k);
    I=eye(m-(k-1));
    vk=x;
    vk(1)=vk(1)+sign(x(1))*norm(x);
    vk=vk/norm(vk);
    V(k:end,k)=vk;
    P=I-2*(vk*vk');
    R(k:end,k:n)=P*R(k:end,k:n);
end
end


function [Q]=house2q(V)
[m, n]=size(V);
assert(m>=n,'The row count must be greater than or equal to the column count')
Q=eye(m);
for k=n:-1:1
    Q(k:m,:) = Q(k:m,:) - 2*V(k:m,k)*(V(k:m,k)'*Q(k:m,:));
end
end




function [Q,R]=MGS(A)
[m, n]=size(A);
assert(m>=n,'The row count must be greater than or equal to the column count')
Q=zeros(m,n);
R=zeros(n,n);
for i=1:n
    vi=A(:,i);
    R(i,i)=norm(vi);
    Q(:,i)=vi/R(i,i);
    for j=i+1:n
        qi=Q(:,i);
        R(i,j)=qi'*A(:,j);
        A(:,j)=A(:,j)-R(i,j)*qi;

    end
end
end

function [Q,R]=GS(A)
[m, n]=size(A);
assert(m>=n,'The row count must be greater than or equal to the column count')
Q=zeros(m,n);
R=zeros(n,n);
for j=1:n
    aj=A(:,j);
    vj=aj;
    for i=1:j-1
        qi=Q(:,i);
        R(i,j)=qi'*aj;
        vj=vj-R(i,j)*qi;
    end
    R(j,j)=norm(vj);
    Q(:,j)=vj/R(j,j);
end
end