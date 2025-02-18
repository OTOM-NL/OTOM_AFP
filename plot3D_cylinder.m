function Mandrel_plot=plot3D_cylinder(R,z_cylinder_end)




if length(R)==3
   
    c1=R(2);
     c2=R(3);
     
else
       c1=R(1);
     c2=R(1);
%       R_cyl=R_cyl(1);
    
end



h=figure (1);
hold on;
set(h,'Visible','off');
% jFrame = get(h, 'JavaFrame');
% jFrame.setMinimized(1);







%plot axes system
x0=0;
y0=0;
z0=0;
plot3(x0+[0, R(1), nan, 0, 0, nan, 0, 0], ...
       y0+[0, 0, nan, 0, R(1), nan, 0, 0], ...
       z0+[0, 0, nan, 0, 0, nan, 0, z_cylinder_end] , 'k--');
text([x0+R(1)/2, x0, x0], [y0, y0+R(1)/2, y0], [z0, z0, z0+z_cylinder_end/2], ['+X';'+Y';'+Z']);
hold on;

if R(1)~=0

a=R(1);
b=R(1);

    
    
    

%plot cylinder
number_of_section=80;
[X,Y,Z] = cylinder(R(1),number_of_section);
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',0.95,'FaceColor',[0.31 0.31 0.31]);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square


number_of_section_sph=80;
[x,y,z] = sphere(number_of_section_sph);
x1=x*a;
y1=y*b;
z1=z*c1;

% lower dome
% h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',0.8,'EdgeColor',[1 1 1]); %'FaceColor','none'
h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',0.8,'EdgeColor',[1 1 1],'FaceColor',[0.94 0.94 0.94]);



% Upper Dome
x2=x*a;
y2=y*b;
z2=(z*c2)+z_cylinder_end;

% h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',0.8,'EdgeColor',[1 1 1]);  %'FaceColor','none'
h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',0.8,'EdgeColor',[1 1 1],'FaceColor',[0.74 0.84 0.84]);


Mandrel_plot=[h1,h2,h3];


elseif R(1)==0
    
    [X,Y]=meshgrid(linspace(-z_cylinder_end,z_cylinder_end,3),linspace(-z_cylinder_end,z_cylinder_end,3));
        
    Tol_representation=1e-5;
    
    Z= 0* X- Tol_representation;
    
    h1=surf(X,Y,Z,'FaceAlpha',0.8,'EdgeColor',[1 1 1],'FaceColor',[0.74 0.84 0.84]);
    
    Mandrel_plot=h1;
    
end

