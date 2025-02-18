
% This program is able to calculate and combine of two axis systems to
% calculate intersectioin points and 1st reflection
%Input: Laser ray point and direction 3D, geometrical paratmers

%Transformation of rays and intersection points should be into axis-1 which
%is considered as a general axis system
% For caluclation in axis-2 all the parameters should be transformed into
% axis-2
% Plot of system should be performed in axis-1 systems
%The output is maximum 2 intersection points and should indicate which
%point are the first one and which are 2nd


%% Calculate for 1st reflection
%reflection from tape (should transformed into axis 1)> check intersection with cylinder
%reflection from cylinder (should transform into axis 2) > check
%intersection with tape (results should be transformed into axis 1)

%%
function counter_ray=General_3D_tape_cylinder(th_y,Tape_Sp,...
    Rxyz,R_cyl,z_cyl_end,tv3,ID,Laser_head,L_xyz0,absorbtion,W_R,H_indentation)




   h1=figure(1);
     h2=figure(2);
     h3=figure(3);
     
     h21=figure(21);
     h100=figure(100);
     
     
     
set(h1,'Position', [0 0 100 100 ]);
set(h2,'Position', [0 0  100 100 ]);
set(h3,'Position', [0 0  100 100 ]);
set(h21,'Position', [0 0  100 100 ]);
set(h100,'Position', [0 0 100 100 ]);


% set(h100,'Position', [680 558 100 100 ]);


N_tape=Tape_Sp(1); % EVEN number of points, should not be changed !!
W_tape=Tape_Sp(2) ; % width of the tape
R_tape=Tape_Sp(3);
L_flat=Tape_Sp(4);
thick_T=Tape_Sp(5); % thickness of the tape
deg_tape=Tape_Sp(6);

h=figure(1);
set(h,'Visible','off');


fileID1 = fopen('Cylinder_ints.txt','w');   % file includes the xyz + ID
fileID2 = fopen('Tape_ints.txt','w');     % file includes the xyz + ID
fileID3 = fopen('Roller_ints.txt','w'); 





 fprintf(fileID1,' X         ,Y         ,Z       ,ID  \r\n');  % 0 means which reflection
  fprintf(fileID2,'X         ,Y         ,Z       ,ID \r\n');  % 0 means which reflection
 fprintf(fileID3,'X         ,Y         ,Z       ,ID \r\n'); 
           
%ID = 0 intersection from laser, 1 from reflected rays

% Define the input ray


Rx= Rxyz(1); %-0.8239; %input('direction of ray in x-direction =');
Ry= Rxyz(2); %-0.1972;% input('direction of ray in y-direction =');
Rz=Rxyz(3); %-0.3455; %input('direction of ray in z-direction =');
 normR=norm([Rx Ry Rz]); % normalize the Rx,Ry, Rz
Rx=Rx/normR;
Ry=Ry/normR;
Rz=Rz/normR;


 
    laser_direction_G=[Rx;Ry;Rz];
% geometrical of cylinder

% R_cyl= 15; %80; %input('The radius of cylinder =');
% z_cyl_end= 50; %input('Enetr the end point of cylinder =');
% plot cylinder
Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);
assignin('base','Mandrel_plot',Mandrel_plot);

%% plot roller

%%

% plot tape ....
% N_tape = 360; % EVEN number of points, should not be changed !!
% W_tape = 6; % width of the tape
% R_tape=5;
% L_flat=R_tape*3;
% 
% thick_T=0.4e-3; % thickness of the tape

%tv  > translation vector
% tv=[0;R_cyl;z_cyl_end/2];  % for the tape
% tv_R=[0;R_cyl+R_tape+thick_T;z_cyl_end/2];  % for center of Roller

% H_indentation=0.00;

 theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree
 
 Nip_Mov =R_tape*sind(theta_ind); 


% current position of the tape/Roller
tv=[0;R_cyl-H_indentation;tv3];  % for the tape
tv_R=[0;R_cyl+R_tape+thick_T-H_indentation;tv3];  % for center of Roller

%Y-Rotation in degrees
% th_y=-170;
% N=360;
% deg_tape=N/8; % degrees for the cylinder section

[n1,n2,n3,x0,y0,z0,h_R,h_T]=Tape_plot_3D(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind);
assignin('base','h_R',h_R);
assignin('base','h_T',h_T);




%n1,n2,n3,x0,y0,z0 are the point and normal of intersection of curve and
%flat part

axis equal;



% an ID which define the laser distribution
% ID=0 uniform, 1 for tophat, 2 for Gaussian, 3 for linear
% laser energy distribution definition
% Ray distribution is uniform

% ID=2;

% The width and long of the laser head

% Laser_head_Ax=2;
% Laser_head_Ay=4;
% 
% Laser_head_nx=10;
% Laser_head_ny=8;




Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);

Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);



counter_ray=Laser_head_nx*Laser_head_ny;

% counter_ray=0;

% L_xyz0=[50 27 39];

% Gauss_Par_X=100;
% Gauss_Par_Y=100;

% ID=[2,Gauss_Par_X,Gauss_Par_Y]



laser_point_Actual_all=Laser_ray_generator (Rx,Ry,Rz,Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,L_xyz0,ID);

h=figure(1);
set(h,'Visible','off');

%%
% absorbtion=0.8;
% First_hit_energy=absorbtion;  % should use better frensel equation or sth else
% second_hit_energy=1-First_hit_energy;

    p1 =   5.965e-11  ;
       p2 =  -1.598e-08  ;
       p3 =   1.695e-06  ;
       p4 =  -8.707e-05  ;
       p5 =    0.002224 ;
       p6 =    -0.02664 ;
       p7 =      0.1532 ;
       p8 =       10.01 ;
       
        
     Abs=@(x) 1-(p1*x^7 + p2*x^6 + p3*x^5 + p4*x^4 + p5*x^3 + ...
                    p6*x^2 + p7*x + p8)/100; % function of absorbtion




%%



for ii=1:counter_ray
    Lx=laser_point_Actual_all(ii,1);
    Ly=laser_point_Actual_all(ii,2);% input('y position of Laser =');
    Lz=laser_point_Actual_all(ii,3);%  input('z position of Laser =');
    

%% 

% This two loops define the number of rays
% for ii=-5:0.5:5
%     for jj=-2:0.6:2
%         counter_ray=counter_ray+1;
        
        % there should a loop, in each time, Lx Ly LZ will be gotten
%         from the laser_ray_generator
        
        
% Lx= 100+ii; % input('x position of Laser = ');    corrdinate of laser
% Ly=42+ii;% input('y position of Laser =');
% Lz=54+jj;%  input('z position of Laser =');
%% 
 laser_source_P_G=[Lx;Ly;Lz];

xyz_int_G=[];
part_index=0;

xyz_int_L=zeros(3,1);
xyz_int_G_C=zeros(3,1);
xyz_int_G_T=zeros(3,1);
xyz_int_G=zeros(3,1);
xyz_int_G_Roller=zeros(3,1);
xyz_int_L_Roller=zeros(3,1);

% Check intersection with cylinder
% a function should be performed to do that

 
  %G_C >> General, Cylinder
[xyz_int_G_C]=Line_cylinder_intersection(laser_source_P_G,laser_direction_G,R_cyl,z_cyl_end);
%if Yes, calculate   > OK go for reflection calculation
%if No, check for tape > transformation of input ray


% if isempty(xyz_int_G),  calculate for the tape intersection
    
    
    %input ray must be in local axis for tape calculations...
    % Working on the transformation matrix
  
    laser_source_P_L=Transformation_G2L (tv,th_y,laser_source_P_G);
    laser_direction_L=Transformation_G2L (zeros(3,1),th_y,laser_direction_G);
    
    [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(laser_source_P_L,laser_direction_L,R_tape,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
        %        intersection_points_L=([x_int_L;y_int_L;z_int_L]);
        xyz_int_G_T=Transformation_L2G (tv,th_y,xyz_int_L);
    
        
        %% for the Roller
        %assume W_R=2*W
        
        laser_source_P_L_Roller=Transformation_G2L (tv_R,th_y,laser_source_P_G);
        
        [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(laser_source_P_L_Roller,laser_direction_L,R_tape,W_R);
                xyz_int_G_Roller=Transformation_L2G (tv_R,th_y,xyz_int_L_Roller);
        %%
        
        
        %     check which xyz_int are more closer to laser source
        
        %%
          dis_C=1e5;
            dis_T=1e5;
            dis_Roller=1e5;
             nm=0; % to check that there is no intersection with all part
       
        if (xyz_int_G_T)
        dis_T=norm(laser_source_P_G-xyz_int_G_T);
        else
            nm=nm-1;
        end
            
          if (xyz_int_G_C)
       dis_C=norm(laser_source_P_G-xyz_int_G_C);
       else
            nm=nm-1;
          end
                
            if (xyz_int_G_Roller)
         dis_Roller=norm(laser_source_P_G-xyz_int_G_Roller);
         else
            nm=nm-1;
            end
                
             [m,n]=min([dis_C,dis_T,dis_Roller]);
             
             if   nm==-3 %~(xyz_int_G_T) | ~(xyz_int_G_C)| ~(xyz_int_G_Roller)  % if no intersection with all aprts
                 n=-1;
             end
             
             if n==1
                %cylinder substrate
%                 \beta is the angle between normal and laser for measuring
%                 angle of incident
                xyz_int_G=xyz_int_G_C;
                [Ref_dir_G,beta]=refl_cylin_plot3D(R_cyl,xyz_int_G,laser_direction_G);
                part_index=1;
                
                First_hit_energy= Abs(beta);
                
%                  fprintf(fileID1,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));  % 0 means which reflection
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);  % 0 means which reflection
              
                 
                 
            elseif n==2
                %with Tape 
                xyz_int_G=xyz_int_G_T;
                [Ref_dir_L,beta]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,laser_direction_L,tv,th_y); %input are in
                %     local axis system
                Ref_dir_G=Transformation_L2G (zeros(3,1),th_y,Ref_dir_L);
                part_index=2;
                
                   First_hit_energy= Abs(beta);
                
%                 fprintf(fileID2,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                 fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);
                
                
            elseif n==3
                %with Roller
                 xyz_int_G=xyz_int_G_Roller;
                [Ref_dir_G,beta]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,laser_direction_L,tv_R,th_y);
 
                part_index=3;
                
                   First_hit_energy= Abs(beta);
                
%                 fprintf(fileID3,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                 fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,laser_point_Actual_all(ii,4)*First_hit_energy);
               

                
            end
          
        
   

  
 plot3D_line_interP(R_cyl,z_cyl_end,Rx,Lx,Ry,Ly,Rz,Lz,xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
 
%% reflection

 
if isempty(xyz_int_G)
%     msgbox ('No intersection with all Parts!');
elseif  part_index==2 
    % continue the reflection calculation..
    xyz_int_G2=zeros(3,1);
        [xyz_int_G2]=Line_cylinder_intersection(xyz_int_G,Ref_dir_G,R_cyl,z_cyl_end);
        %intersection points
        if xyz_int_G2
plot3( xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3),'y.'  );
% text([xyz_int_G2(1)],[xyz_int_G2(2)],[xyz_int_G2(3)], ['P2']);

% intersection with cylinder, the reflected ray comes from tape
%                  fprintf(fileID1,' %f %f %f , 1 \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));


  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G2,Ref_dir_G);


second_hit_energy=(1-First_hit_energy)*Abs(beta2);


                 
 fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), laser_point_Actual_all(ii,4)*second_hit_energy);
 
                


        end
    
        
      elseif  part_index==3   % maybe hit the tape or mandrel
            xyz_int_L2=Transformation_G2L (tv,th_y,xyz_int_G); % checking for the tape
          Ref_dir_L2=Transformation_G2L (zeros(3,1),th_y,Ref_dir_G);
          
              %for the tape
     [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
      xyz_int_G_T=Transformation_L2G (tv,th_y,xyz_int_L);
      
        % for the mandrel

        [xyz_int_G_M]=Line_cylinder_intersection(xyz_int_G,Ref_dir_G,R_cyl,z_cyl_end);
      
      
        
        
        
                     
        if ~(isempty(xyz_int_G_M) || isempty(xyz_int_G_T) )  % intersection with 3 parts
            
            dis_T=norm(xyz_int_G-xyz_int_G_T);
            dis_M=norm(xyz_int_G-xyz_int_G_M);
            [m,n]=min([dis_T,dis_M]);
            
            
          
              if n==1
                %with Tape 
                 plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));


  [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,th_y); 
  %input are in   local axis system
 second_hit_energy=(1-First_hit_energy)*Abs(beta2);
               
          fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
    
              elseif n==2
                %with Mandrel
                    plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));

  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
second_hit_energy=(1-First_hit_energy)*Abs(beta2);

           fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3), laser_point_Actual_all(ii,4)*second_hit_energy);
              
              end
            
              
              elseif xyz_int_G_T    
      plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));

  [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,th_y); 
  %input are in   local axis system
 second_hit_energy=(1-First_hit_energy)*Abs(beta2);

               
 fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
     
     elseif xyz_int_G_M
   %with Mandrel
                    plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));


 [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
second_hit_energy=(1-First_hit_energy)*Abs(beta2);

%                
 fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),laser_point_Actual_all(ii,4)*second_hit_energy);
               
               
    
    end
        
        
        
        
        
elseif part_index==1;
        
      xyz_int_L2=Transformation_G2L (tv,th_y,xyz_int_G);
          xyz_int_L2_R=Transformation_G2L (tv_R,th_y,xyz_int_G);
    Ref_dir_L2=Transformation_G2L (zeros(3,1),th_y,Ref_dir_G);
    
    xyz_int_L=zeros(3,1);
    
    %for the tape
     [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
      xyz_int_G=Transformation_L2G (tv,th_y,xyz_int_L);
      
      % for the roller
         [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(xyz_int_L2_R,Ref_dir_L2,R_tape,W_R);
                xyz_int_G_Roller=Transformation_L2G (tv_R,th_y,xyz_int_L_Roller);
                
                
            
                
                
        if ~(isempty(xyz_int_G_Roller) || isempty(xyz_int_G) )  % intersection with 3 parts
            
            dis_T=norm(xyz_int_L2-xyz_int_L);
            dis_Roller=norm(xyz_int_L2_R-xyz_int_L_Roller);
            [m,n]=min([dis_T,dis_Roller]);
            
            
            %%
              if n==1
                %with Tape 
                 plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID2,' %f %f %f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));


  [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,th_y); 
  %input are in   local axis system
 second_hit_energy=(1-First_hit_energy)*Abs(beta2);


               
               fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), laser_point_Actual_all(ii,4)*second_hit_energy);
              
    
              elseif n==2
                %with Roller
                    plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));

  %input are in   local axis system

 
 
  [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,th_y);
   second_hit_energy=(1-First_hit_energy)*Abs(beta2);

               
 fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3), laser_point_Actual_all(ii,4)*second_hit_energy);
               
              end
            
              elseif xyz_int_L     
         plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
   
         % intersection with tape, the reflected ray comes from cylinder
%      fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));

 [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,th_y); 
 second_hit_energy=(1-First_hit_energy)*Abs(beta2);


     
        fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), laser_point_Actual_all(ii,4)*second_hit_energy);
     
     
     elseif xyz_int_G_Roller
         
    %with Roller
                    plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                 
%                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));

  [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,th_y);
   second_hit_energy=(1-First_hit_energy)*Abs(beta2);

               
 fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
    
    end
    
        
end
 
    
    
    
    
    end
% end

grid  on;
% axis([-2*R_cyl 2*R_cyl -2*R_cyl 2*R_cyl -z_cyl_end/2 1.5*z_cyl_end]);
view([180 0]);



	 fclose(fileID1);
     	 fclose(fileID2);
     	 fclose(fileID3);
%first intersection are black stars
% legend('yellow for second intersections')
% reflection from cylinder are red
% reflection from the tape are blue
% 

%% Changing the graphical object, Visibility or deleting

Graphics_laser=findobj(h,'color','g');
assignin('base','Graphics_laser',Graphics_laser);
% set(Graphics_laser,'Visible','off');
Graphics_reflection1=findobj(h,'color','r');
Graphics_reflection2=findobj(h,'color','b');
Graphics_reflection=[Graphics_reflection1;Graphics_reflection2];
assignin('base','Graphics_reflection',Graphics_reflection);

reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
assignin('base','reflection_int',reflection_int);

Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
assignin('base','Laser_int_points',Laser_int_points);

Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
assignin('base','Normal_lines',Normal_lines);

% Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
% assignin('base','Thermal_points',Thermal_points);




% set(Graphics_reflection,'Visible','off');

   
  

% Create a push button
% btn3 = uibutton(gca,'push',...
%                'Position',[250 20 90 30],...
%                'Value', true,'ButtonPushedFcn', @(btn3,event) Visible_plotButtonPushed(btn3,Graphics_laser));




rb1 = uicontrol('Style', 'radiobutton','Position',[20 50 90 40],'Value',true,... 
    'String','Rays','Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Graphics_laser));

  btn1 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser rays',...
        'Position', [20 20 90 40],...
        'Callback', 'delete(Graphics_laser)');  


rb2 = uicontrol('Style', 'radiobutton','Position',[150 50 90 40],'Value',true,... 
    'String','Reflections','Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Graphics_reflection));

     btn2 = uicontrol('Style', 'pushbutton', 'String', 'Clear Reflections',...
        'Position', [150 20 90 40],...
        'Callback', 'delete(Graphics_reflection)');    
    


rb3 = uicontrol('Style', 'radiobutton','Position',[280 50 90 40],'Value',true,... 
    'String','Mandrel','Callback', @(rb3,event) Visible_plotButtonPushed(rb3,Mandrel_plot));

     btn3 = uicontrol('Style', 'pushbutton', 'String', 'Clear Mandrel',...
        'Position', [280 20 90 40],...
        'Callback', 'delete(Mandrel_plot)');    


rb4 = uicontrol('Style', 'radiobutton','Position',[400 50 90 40],'Value',true,... 
    'String','Roller','Callback', @(rb4,event) Visible_plotButtonPushed(rb4,h_R));

 btn4 = uicontrol('Style', 'pushbutton', 'String', 'Clear Roller',...
        'Position', [400 20 90 40],...
        'Callback', 'delete(h_R)');   


% rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,... 
%     'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,Surf_T));


rb6 = uicontrol('Style', 'radiobutton','Position',[640 50 90 40],'Value',true,... 
    'String','1st reflections','Callback', @(rb6,event) Visible_plotButtonPushed(rb6,reflection_int));
 btn6 = uicontrol('Style', 'pushbutton', 'String', 'Clear reflection_int',...
        'Position', [640 20 90 40],...
        'Callback', 'delete(reflection_int)');   


rb7 = uicontrol('Style', 'radiobutton','Position',[760 50 90 40],'Value',true,... 
    'String','Laser_int','Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Laser_int_points));

 btn7 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser_int_points',...
        'Position', [760 20 90 40],...
        'Callback', 'delete(Laser_int_points)');   


% 
% rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,... 
%     'String','Thermal_P','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Thermal_points));

rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,... 
    'String','Normal lines','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Normal_lines));

 btn8 = uicontrol('Style', 'pushbutton', 'String', 'Clear Normal_lines',...
        'Position', [900 20 90 40],...
        'Callback', 'delete(Normal_lines)');   




function Visible_plotButtonPushed(rb1,Graphics_laser)
%         x = linspace(0,2*pi,100);
%         y = sin(x);
%         plot(ax,x,y)
if rb1.Value ==false
set(Graphics_laser,'Visible','off');
else 
    set(Graphics_laser,'Visible','on');
end

