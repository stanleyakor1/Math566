clc;clear all;close all;


%%%%%% Question 1 %%%%%%%%%%
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

[Q2,R2]=GS(A);
c3=R2\(Q2'*f);

%%% 1d %%%%%
[V3,R3]=house(A);
Q3=house2q(V3);
c4=R3\(Q3'*f);


%%% 1e %%%%%
[Q4,R4]=qr(A);
c5=R4\(Q4'*f);



%%% 1e %%%%%
[Q5,R5]=qr(A);
c6=R5\(Q5'*f);





%%%%%% Question 6 %%%%%%%%%%

x=-1.920:0.001:2.080;
r=zeros(size(x));
r_q=(x-2).^9;
for i=1:length(x)
    r(i)=p(x(i));
end
plot(x,r,'--or',LineWidth=1), hold on
plot(x,r_q,'k',LineWidth=1), hold on
axis tight;


%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  r = p(x)
    r = x^9 -18*x^8 + 144*x^7 - 672*x^6 + 2016*x^5 -4032*x^4 + 5376*x^3 - 4608*x^2 + 2304*x -512;
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