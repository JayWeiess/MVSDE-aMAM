function [d_alpha,alpha,Ck,w]=repara(myFx,P,Psi,t,r)
%% myFx is b(x,y,u,v)
% P is \varphi, Psi is \psi
% return value of w,alpha,Ck
m=length(t);

b=feval(myFx,P(1,:),P(2,:),Psi(1,:),Psi(2,:));
b=sqrt(sum(b.^2,1)).^r;
dt=t(2:m)-t(1:m-1);
% ba=4*sqrt(sum(feval(myFx,(P(1,1:m-1)+P(1,2:m))./2,(P(2,1:m-1)+P(2,2:m))./2,(Psi(1,1:m-1)+Psi(1,2:m))./2,(Psi(2,1:m-1)+Psi(2,2:m))./2).^2,1)).^r+...
%     sqrt(sum(feval(myFx,P(1,1:m-1),P(2,1:m-1),Psi(1,1:m-1),Psi(2,1:m-1)).^2,1)).^r+sqrt(sum( feval(myFx,P(1,2:m),P(2,2:m),Psi(1,2:m),Psi(2,2:m)).^2,1)).^r;
% 
% ba=ba./6;
ba=sqrt(sum(feval(myFx,(P(1,1:m-1)+P(1,2:m))./2,(P(2,1:m-1)+P(2,2:m))./2,(Psi(1,1:m-1)+Psi(1,2:m))./2,(Psi(2,1:m-1)+Psi(2,2:m))./2).^2,1)).^r;
d_alpha=dt.*ba;
Ck=sum(d_alpha);

d_alpha=d_alpha./Ck;
w=b./Ck;

alpha=zeros(1,m);
for i=2:m
    alpha(i)=alpha(i-1)+d_alpha(i-1);
end
