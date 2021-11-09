%% Plot the equi-potential contour for MVSDE

clear
clc
be=10;
x_start=-1.5;
x_end=1.5;
y_start=0.5;
y_end=-0.5;

dt=1e-2;
T=6;
N=10;

% d=(y_end-y_start)/N;
%  y0=y_start:d:y_end;
%  x0=ones(1,N+1).*6;

V=@(u,v)  [u-u.^3-be*u.*v.^2; -(1+u.^2).*v];
F1 =@(u,v)  [u;-v];
Da=@(u,v) u^2+v^2+eta;
F2=@(u,v) [-v/Da(u,v); u/Da(u,v)]./(2*pi);

CValue_x=zeros(4*(N+1),T/dt);
CValue_y=zeros(4*(N+1),T/dt);
d=(x_end-x_start)/N;
x0=x_start:d:x_end;
y0=ones(1,N+1).*(y_start);

for i=1:N+1
[x1,y1]=Traj_BS(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
CValue_x(i,:)=x1;
CValue_y(i,:)=y1;
hold on
x1=[];
y1=[];
end


d=(x_end-x_start)/N;
x0=x_start:d:x_end;
y0=ones(1,N+1).*(y_end);

for i=1:N+1
[x1,y1]=Traj_BS(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
CValue_x(i+N+1,:)=x1;
CValue_y(i+N+1,:)=y1;
hold on
x1=[];
y1=[];
end
%% -----------------------------------------------------------
x0=[];
y0=[];

d=2*(y_end-y_start)/N;
y0=y_start:d:y_end;
x0=ones(1,N/2+1).*x_start;

for i=1:N/2+1
[x1,y1]=Traj_BS(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
CValue_x(i+2*(N+1),:)=x1;
CValue_y(i+2*(N+1),:)=y1;
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
[x1,y1]=Traj_BS(x0(i),y0(i),T,dt);
plot(x1,y1,'color',[0.65,0.65,0.65]);
CValue_x(i+3*(N+1),:)=x1;
CValue_y(i+3*(N+1),:)=y1;
hold on
x1=[];
y1=[];
end

