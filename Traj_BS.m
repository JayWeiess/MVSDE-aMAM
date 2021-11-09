function [x,y]=Traj_BS(x0,y0,T,dt)

be=10;
eta=1e-2;
V=@(u,v)  [u-u.^3-be*u.*v.^2; -(1+u.^2).*v];
Da=@(u,v) u^2+v^2+eta;
F1=@(u,v) [-v/Da(u,v); u/Da(u,v)]./(2*pi);
xe=-1;
ye=0;

t=0:dt:T;
n=length(t);
x=zeros(1,n-1);
y=zeros(1,n-1);

x(1)=x0;
y(1)=y0;
for i=2:n-1
  a=(V(x(i-1),y(i-1))-F1(x(i-1)-xe,y(i-1)-ye))*dt+[x(i-1);y(i-1)];
  x(i)=a(1);
  y(i)=a(2);
end