
load('Contour_BS_dt=1e-4.mat')
N=15;

for i=1:4*(N+1)
plot(CValue_x(i,:),CValue_y(i,:),'color',[0.65,0.65,0.65]);
hold on
end