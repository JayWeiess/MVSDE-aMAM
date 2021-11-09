%% Plot the equi-potential contour for MVSDE

clear
clc
be=1;
x_start=-1.5;
x_end=1.5;
y_start=1;
y_end=-0.2;

dt=1e-2;
T=10;
N=10;

% d=(y_end-y_start)/N;
%  y0=y_start:d:y_end;
%  x0=ones(1,N+1).*6;

V=@(u,v)  [u-u.^3-be*u.*v.^2; -(1+u.^2).*v];
F1 =@(u,v)  [u;-v];
Da=@(u,v) u^2+v^2+eta;
F2=@(u,v) [-v/Da(u,v); u/Da(u,v)]./(2*pi);

% t=0:dt:T;
% n=length(t);
% x1=zeros(1,n-1);
% y1=zeros(1,n-1);
% 
% x1(1)=x0;
% y1(1)=y0;
% for i=2:n-1
%   a=V(x1(i-1),y1(i-1))*dt+[x1(i-1);y1(i-1)];
%   x1(i)=a(1);
%   y1(i)=a(2);
% end

d=(x_end-x_start)/N;
x0=x_start:d:x_end;
y0=ones(1,N+1).*y_start;
for i=1:N+1
[x1,y1]=Traj_Line(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
%plot(x1,y1,'*')
hold on
x1=[];
y1=[];
end
%% ----------------------------------------------------------
d=(x_end-x_start)/N;
x0=x_start:d:x_end;
y0=ones(1,N+1).*(y_end);


for i=1:N+1
[x1,y1]=Traj_Line(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
hold on
x1=[];
y1=[];
end
%% -----------------------------------------------------------
x0=[];
y0=[];

d=(y_end-y_start)/N;
y0=y_start:d:y_end;
x0=ones(1,N+1).*x_start;

for i=1:N+1
[x1,y1]=Traj_Line(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
hold on
x1=[];
y1=[];
end
%% -----------------------------------------------------------------
x0=[];
y0=[];
d=(y_end-y_start)/N;
y0=y_start:d:y_end;
x0=ones(1,N+1).*x_end;

for i=1:N+1
[x1,y1]=Traj_Line(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
hold on
x1=[];
y1=[];
end

