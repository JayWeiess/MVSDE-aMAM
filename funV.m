function dydt=funV(~,y)
u=y(1);
v=y(2);
dydt=[u-u.^3-u.*v.^2; -(1+u.^2).*v];