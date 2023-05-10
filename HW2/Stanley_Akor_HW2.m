clear all;close all;clc

%%%% Question 2(b) %%%

A = [14 -2; -8 19];
[U,S,V]=svd(A);
S1=S;
S=diag(S);
n = 100; z = exp(1i*(0:n)*2*pi/n);
x = [real(z);imag(z)];
% Image of the unit circle under A
y = A*x;
% Plot of the unit circle
figure(1);
plot(x(1,:),x(2,:),'b-'), axis square, title('Unit circle'), hold on
plot([0,V(1,1)],[0,V(2,1)],'k')
plot([0,V(1,2)],[0,V(2,2)],'k')
grid on


figure(2)
plot(y(1,:),y(2,:),'r-'), axis square, title('Image under $A$ circle'),hold on
plot([0,S(1)*U(1,1)],[0,S(1)*U(2,1)],'k')
plot([0,S(2)*U(1,2)],[0,S(2)*U(2,2)],'k')
grid on


%%%%%%%% Question 2(c) %%%%%%%%%%%%%%%%%

A_2=max(S)

A_F = sqrt(trace(S1'*S1))