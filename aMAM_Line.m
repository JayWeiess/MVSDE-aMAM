function [db1,db2,b,gradx_b,grady_b]=aMAM_Line
%% a two dimensional case 
% b(x,y,u,v)=V(x,y)-F(x-u,y-v)
be=10;
V=@(u,v) [u-u.^3-be*u.*v.^2; -(1+u.^2).*v];
F=@(u,v) [u;-v];
dV1=@(u,v) [- 3*u.^2-be*v.^2 + 1; -2*u.*v];                                       %derviative w.r.t. u of V(u,v)
dV2=@(u,v)  [-2*be*u.*v; - u.^2 - 1];                                                      %derviative w.r.t. v of V(u,v)
dF1=@(u,v) [1;0];                                                                                         %derviative w.r.t. u of F(u,v)
dF2=@(u,v) [0;-1];                                                                                        %derviative w.r.t. v of F(u,v)

b=@(x,y,u,v) V(x,y)-F(x-u,y-v);                                                                    % b(y,x)

db1=@(x,y,u,v) dV1(x,y)-dF1(x-u,y-v);                                                      %derviative w.r.t. x of b(x,y,u,v)
db2=@(x,y,u,v) dV2(x,y)-dF2(x-u,y-v);                                                      %derviative w.r.t. y of b(x,y,u,v)
gradx_b=@(x,y,u,v) [dF1(x-u,y-v),dF2(x-u,y-v)];                                        % gradient w.r.t. x of  b(y,x)
grady_b=@(x,y,u,v) [dV1(x,y)-dF1(x-u,y-v),dV2(x,y)-dF2(x-u,y-v)];      % gradient w.r.t. y of  b(y,x)
