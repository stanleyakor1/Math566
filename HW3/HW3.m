clc;clear all;




%%%%%% Set up Vandermonde matrix %%%%%%%%%%%%%%%%%%
% Question 1

m=100;
n=15;
A=zeros(m,n);
A(:,1)=ones(m,1);
for j=1:m
    tj=(j-1)/(m-1);
    for i=2:n
        A(j,i)=tj^(i-1);
    end
end

%1a
[Q,R]=GS(A);
%1b
[q1,r1]=MGS(A);
%1c
norm(A-q1*r1,'inf')
norm(eye(n)-q1'*q1)
norm(A-Q*R,'inf')
norm(eye(n)-Q'*Q)


%Question 3
b=[1 2 3;4 5 6; 7 8 7; 4 2 3; 4 2 2];

[V,R]=house(b)
[Q]=house2q(V)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%           Householder routine based on the Numerical Linear  book by    %
%                       Trefethen and Bau                                 %  
%                                                                         % 
%                                                                         %  
%                                                                         %  
%                                                                         %  
%                                                                         %  
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
% Classical and Modified Gram-Schdmidt routine based on the Numerical     %
%                   Linear  book by Trefethen and Bau                     %  
%                                                                         % 
%                                                                         %  
%                                                                         %  
%                                                                         %  
%                                                                         %  
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

