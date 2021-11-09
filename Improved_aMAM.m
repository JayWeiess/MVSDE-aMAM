%% Compute the most likely transition path for a SDE
%% Improved aMAM
% We first try linear interpolation
clear 
clc

 [db1,db2,b,gradx_b,grady_b]=aMAM_BS;
                                                                                                     %Define 1/\tau to be 10
r=0.5;
K=15000;                                         %The iteration times 
N=400;                                           % The partition size
T=40;                                            %The original time scale
tau=400;   
P=zeros(2,N+1);                      %The matrix that stores the function {\varphi_i^k}, where the (:,i+1)-th element of this matrix stores the value of {\varphi_i}
x1=[-1;0];  x2=[1; 0];   y0=[-1;0];   x_s=[0.009145;-0.152703];         %Define the end points
P1=zeros(2,N+1);                     % The matrix that stores  {\varphi_i^{\prime k}}
%P2=zeros(2,N+1);                    % The matrix that stores  {\tilde{\varphi}_i}
P3=zeros(2,N+1);                    % An intermediate matrix

%Psi=zeros(2,N+1);

d1=T/N;
%d=(x2(1)-x1(1)) /N;
 % generate a line with ends x1 and x2 as the initial value
 d=(x2(1,1)-x1(1,1))/N;
P(1,:)=x1(1):d:x2(1);  

P(1,:)=x1(1,1):d:x2(1,1);
P(2,:)=interp1([x1(1,1), x_s(1,1),x2(1,1)],[x1(2,1), x_s(2,1),x2(2,1)],P(1,:));
P(:,1)=x1; P(:,N+1)=x2; 

% P(:,1)=x1; P(:,N+1)=x2; 
% P(2,:)=-0.5*P(1,:).^2+0.5;

% P(2,2:1:N/2)=sqrt(0.25-(P(1,2:1:N/2)+0.5).^2);
% P(2,N/2+1:1:N)=sqrt(0.25-(P(1,N/2+1:1:N)-0.5).^2);

t=0:d1:T;                                % Define a initial time partition
[~,Psi]=ode45(@funV,t,y0);
Psi=Psi';
%P(2,2:1:N)=sqrt(1-(P(1,2:1:N)).^2);
e=[1,0;0,1];                        % Define a 2-dimensional identical matrix
da=1/N;


%% k given
 for j=1:K

[d_alpha,alpha,Ck,w]=repara(b,P,Psi,t,r);
aa=0:da:1;

P1(1,1:N)=interp1(alpha,P(1,:),aa(1:N));                            % Get the value of \tilde{\varphi}_i^k            
P1(2,1:N)=interp1(alpha,P(2,:),aa(1:N));
P1(:,N+1)=x2;

t1=interp1(alpha,t,aa(1:N));                                                % Get the value of t_i^{k+1}
t1(N+1)=T;
%[alpha,Ck,w]=repara(b,P1,Psi,t1,r);

[~,Psi]=ode45(@funV,t1,y0);
Psi=Psi';

P2(1,:)=(P1(1,3:N+1)-P1(1,1:N-1))./(2*da);
P2(2,:)=(P1(2,3:N+1)-P1(2,1:N-1))./(2*da);
dw=(w(3:N+1)-w(1:N-1))./(2*da);
t=t1;
%% Solve the linear system Ax=b

A1=zeros(1,2*(N+1));                                                                           % define the elements of the matrix
A2=zeros(1,2*N-2);
A3=zeros(1,2*N+2);

for i=1:N-1
    A1(2*i+1:2*(i+1))=[-2*w(i+1)^2/da^2-tau,-2*w(i+1)^2/da^2-tau];
    A2(2*i-1:2*i)=[w(i+1)^2/da^2,w(i+1)^2/da^2] ;
    A3(2*i+1:2*i+2)=(grady_b(P1(1,i+1),P1(2,i+1),Psi(1,i+1),Psi(2,i+1)))'*b(P1(1,i+1),P1(2,i+1),Psi(1,i+1),Psi(2,i+1))+gradx_b(P1(1,i+1),P1(2,i+1),Psi(1,i+1),Psi(2,i+1))*b(Psi(1,i+1),Psi(2,i+1),Psi(1,i+1),Psi(2,i+1))-...
        ((grady_b(P1(1,i+1),P1(2,i+1),Psi(1,i+1),Psi(2,i+1)))'-grady_b(P1(1,i+1),P1(2,i+1),Psi(1,i+1),Psi(2,i+1)) +e.*dw(i))*P2(:,i).*w(i+1)-tau*P1(:,i+1);
end
A1(1:2)=[1,1];
A1(2*N+1:2*N+2)=[1,1];
A3(1:2)=x1'; A3(2*N+1:2*N+2)=x2';

v=[A1,A2,A2];                                                                                            %Store the elements value of the sparse matrix S
row=[1:1:2*N+2, 3:1:2*N,3:1:2*N];                                                      %Store the row index of the sparse matrix S
col=[1:1:2*N+2,1:1:2*N-2,5:1:2*N+2];                                                %Store the coloumn index of the sparse matrix S
A=sparse(row,col,v);                                                                                %Define the sparse matrix A
S=A\A3';                                                                                                     %Solve the linear equation;

for i=1:N+1
    P3(:,i)=[S(2*i-1,1);S(2*i,1)];                                                                 %Store the value of \tilde{\varphi}_i^k in P2
end
P=P3;

 end

figure(1)
plot(P(1,:),P(2,:))
hold on
plot(P(1,:),P(2,:),'r*')

%  figure(1)
%  
%  plot(P(1,2:N),P1(2,:));
% xlabel('x')
% ylabel('$\varphi^\prime$','interpreter','latex','FontSize', 12)
% title(' Plot of $\varphi^\prime$','interpreter','latex', 'FontSize', 18)
% 
% 
% %   figure(2)
% %   plot(P3(1,:));
% %   
%   figure(2)
%   plot(P(1,:),P(2,:),'b*');
%   xlabel('x')
%   ylabel('$\varphi$','interpreter','latex','FontSize', 12)
%   title(' Plot of $\varphi$','interpreter','latex', 'FontSize', 18)
%   
%   
%   figure(3)
%   plot(P(1,1:N+1),lambda,'r*')
%   xlabel('x')
%    ylabel('$\lambda$','interpreter','latex','FontSize', 12)
%   title(' Plot of $\lambda$','interpreter','latex', 'FontSize', 18)
%  %axis equal
% %   figure(2)
% %   plot(P(2,:));
%save("Interaction_Line.mat",'P');
