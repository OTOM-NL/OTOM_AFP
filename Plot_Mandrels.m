function Mandrel_plot=Plot_Mandrels(R,z_cylinder_end,Type,finer,c1,c2,a1,b1,RGB,opaqueness)



%%
%test
% R=1;
% z_cylinder_end=2;
% Type=4;
% finer=40;
% c1=.5;
% c2=2;
% a1=1;
% b1=2;



%%

if RGB(1) >1 || RGB(2) >1 || RGB(3) >1 || RGB(1) < 0 || RGB(2) < 0 || RGB(3) < 0
    
    
    % should be modified to give a error msg
    
    
%     msgBox('RGB should be between 0 and 1');
%     errordlg('RGB should be between 0 to 1');
%     h = msgbox('Operation Completed');
%     h = msgbox('Invalid Value', 'Error','error');
    RGB=[0.31 0.31 0.31];
end





% figure (1)

%R is the radious of the cylinderical part
% z_cylinder_end is the length of cylinder
% Type which indicates 3 types of mandrel: cylinder, cylinder + domes,
% domes
% finer number of devision to have a finer shape

%% Define and plot axis system

%plot axes system
x0=0;
y0=0;
z0=0;
plot3(x0+[0, R, nan, 0, 0, nan, 0, 0], ...
       y0+[0, 0, nan, 0, R, nan, 0, 0], ...
       z0+[0, 0, nan, 0, 0, nan, 0, z_cylinder_end] , 'k--');
text([x0+R/2, x0, x0], [y0, y0+R/2, y0], [z0, z0, z0+z_cylinder_end/2], ['+X';'+Y';'+Z']);
hold on;
  view([0 0]);


%%
switch Type
    case 1
        % Cylinder

%plot cylinder
% number_of_section=50;
number_of_section=finer;


[X,Y,Z] = cylinder(R,number_of_section);
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',opaqueness,'FaceColor',RGB,'EdgeColor',[1 1 1]);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square;

hold on;
C=Z*0;

h2=fill3(X(1,:),Y(1,:),Z(1,:),C(1,:),'FaceAlpha',opaqueness,'FaceColor',RGB);
h3=fill3(X(2,:),Y(2,:),Z(2,:),C(2,:),'FaceAlpha',opaqueness,'FaceColor',RGB);



Mandrel_plot=[h1,h2,h3];

axis equal;   

    case 2
        %cylinder + spherical domes
        
        %plot cylinder
% number_of_section=50;
number_of_section=finer;


[X,Y,Z] = cylinder(R,number_of_section);
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',opaqueness,'FaceColor',RGB,'EdgeColor',[1 1 1]);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square

        
        % number_of_section_sph=40;
number_of_section_sph=finer;


[x,y,z] = sphere(number_of_section_sph);
x1=x*R;
y1=y*R;
z1=z*R;

h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB); %'FaceColor','none'


x2=x*R;
y2=y*R;
z2=(z*R)+z_cylinder_end;

h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB);  %'FaceColor','none'


Mandrel_plot=[h1,h2,h3];
axis equal;   


   case 3  
        
%         Cylinder + different height domes


        %cylinder + spherical domes
        
        %plot cylinder
% number_of_section=50;
number_of_section=finer;


[X,Y,Z] = cylinder(R,number_of_section);
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',opaqueness,'FaceColor',RGB,'EdgeColor',[1 1 1]);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square

        
        % number_of_section_sph=40;
number_of_section_sph=finer;


[x,y,z] = sphere(number_of_section_sph);
x1=x*R;
y1=y*R;
z1=z*c1;

h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB); %'FaceColor','none'


x2=x*R;
y2=y*R;
z2=(z*c2)+z_cylinder_end;

h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB);  %'FaceColor','none'


Mandrel_plot=[h1,h2,h3];
        axis equal;   
        
    case 4
    
%         Ellipsoidical-cylinder
      
%plot cylinder
% number_of_section=50;
number_of_section=finer;


[X,Y,Z] = cylinder(1,number_of_section);
X=X*a1;
Y=Y*b1;
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square;

hold on;
C=Z*0;

h2=fill3(X(1,:),Y(1,:),Z(1,:),C(1,:),'FaceAlpha',opaqueness,'FaceColor',RGB);
h3=fill3(X(2,:),Y(2,:),Z(2,:),C(2,:),'FaceAlpha',opaqueness,'FaceColor',RGB);



Mandrel_plot=[h1,h2,h3];
 axis equal;   

    case 5
        
%         Ellipsoidical-cylinder + domes


      
%plot cylinder
% number_of_section=50;
number_of_section=finer;


[X,Y,Z] = cylinder(1,number_of_section);
X=X*a1;
Y=Y*b1;
Z=Z*z_cylinder_end;
h1=mesh(X,Y,Z,'FaceAlpha',opaqueness,'FaceColor',RGB,'EdgeColor',[1 1 1]);
% assignin('base','Graphics_reflection',Mandrel_plot);
% colormap copper
% camlight(120,150)
% shading flat
axis square;

hold on;




number_of_section_sph=finer;


[x,y,z] = sphere(number_of_section_sph);
x1=x*a1;
y1=y*b1;
z1=z*c1;

h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB); %'FaceColor','none'


x2=x*a1;
y2=y*b1;
z2=(z*c2)+z_cylinder_end;

h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB);  %'FaceColor','none'


Mandrel_plot=[h1,h2,h3];



axis equal;   
        
    case 6
        
	
        
            
        % number_of_section_sph=40;
number_of_section_sph=finer;




[x,y,z] = sphere(number_of_section_sph);
x1=x*a1;
y1=y*b1;
z1=z*c1;


h2=surf(x1(1:(number_of_section_sph/2)+1,1:end),y1(1:(number_of_section_sph/2)+1,1:end),z1(1:(number_of_section_sph/2)+1,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB); %'FaceColor','RGB'




x2=x*a1;
y2=y*b1;
z2=(z*c2);

h3=surf(x2((number_of_section_sph/2)+1:end,1:end),y2((number_of_section_sph/2)+1:end,1:end),z2((number_of_section_sph/2)+1:end,1:end),'FaceAlpha',opaqueness,'EdgeColor',[1 1 1],'FaceColor',RGB);  %'FaceColor','none'


Mandrel_plot=[h2,h3];
    axis equal;  
  
        
end

