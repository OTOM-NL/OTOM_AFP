


% 
% z=-2:.3:2;
% % x=z.^2;
% % y=z.^1.5;
% 
% [x,y]=meshgrid(z,z)
% 
% f=x.^2+y.^2 -4
% 
% mesh(x,y,f)
% 

t = 0:pi/10:2*pi;
figure
[X,Y,Z] = cylinder(2+cos(t),80);

 X=X.* (Z.^1.5);
%  Y=Y.^1.5;
Z=Z*10;

surf(X,Y,Z,'Linestyle','--')
axis equal
colormap cool