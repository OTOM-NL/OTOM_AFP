
% This program is able to calculate and combine of two axis systems to
%Transformation of rays and intersection points should be into axis-1 which

% Calculate for 1st reflection


function [counter_ray,Rot_Roller_axis,tv]=General_3D_Optical_UOT(th_y,Tape_Sp,...
    Rxyz,R_cyl,z_cyl_end,Roller_Pos_TV,ID,Laser_head,absorbtion_waste,W_R,H_indentation,...
    Graphic_chekbox,...
    BRDF_mode,jobname,...
    CV_mesh,LaserH_mesh,delta_new_all,pos_new_all,nip_point_M_all,Rot_Roller_axis_all,...
    UOT_pathfile,text_status,Divergence_factor)


%%
% bump in an assumed distance >> 25 April 2019

%%

       
a=R_cyl(1);
b=R_cyl(1);
 if length(R_cyl)==3
   
    c1=R_cyl(2);
     c2=R_cyl(3);
else
       c1=R_cyl(1);
     c2=R_cyl(1);
%       R_cyl=R_cyl(1);
    
 end

 
 
%   tic

% tv3 AND absorbtion_waste are useless here


% if Graphic_chekbox(7)
%     h1=figure(1);
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%     set(h1,'Position', [0 0 100 100 ]);
% end
%
%
% if Graphic_chekbox(8)
%     h2=figure(2);
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%     set(h2,'Position', [0 0  100 100 ]);
% end

% if Graphic_chekbox(9)
%     h3=figure(3);
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%     set(h3,'Position', [0 0  100 100 ]);
%
%     h100=figure(100);
%     set(h100,'Position', [0 0 100 100 ]);
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%
% end


% if  Graphic_chekbox(3) &&  Graphic_chekbox(5)
%     h21=figure(21);
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%
%     set(h21,'Position', [0 0  100 100 ]);
% end
%
%



% an ID which define the laser distribution

Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);

Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);



counter_ray=Laser_head_nx*Laser_head_ny;


% laser_point_Actual_all(ii,4)
% laser_point_Actual_all=Laser_ray_generator (Rx,Ry,Rz,Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,L_xyz0,ID,Laser_head_Rot);


% consider as Uniform distribution
ID(1)=0;
Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);
% laser_point_Actual_all(:,4)
    %% >> Modified on 11 Aug 2018  absorption fault
%     p1 =   -6.55e-08 ;
%     p2 =    1.56e-05  ;
%     p3 =   -0.001348 ;
%     p4 =     0.05033  ;
%     p5 =       0.198 ;
%     
%     Abs=@(x) p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;

%% 25 April 2019
% General model Rat32:
%      f(x) = (p1*x^3 + p2*x^2 + p3*x + p4) /
%                (x^2 + q1*x + q2)
% Coefficients (with 95% confidence bounds):
%        p1 =      0.2602 ;
%        p2 =      -46.97 ;
%        p3 =        6542 ;
%        p4 =        5073 ;
%        q1 =        4065 ;
%        q2 =   7.148e+04 ;

% Goodness of fit:
%   SSE: 0.003576
%   R-square: 0.995
%   Adjusted R-square: 0.993
%   RMSE: 0.01726
%     
%  Abs = @(x)  (p1*x^3 + p2*x^2 + p3*x + p4) /  (x^2 + q1*x + q2);
    
    
    %% Should be modified
%     if absorbtion_waste(1) ==-1
%         Abs_mandrel=@(x) (p1*x^3 + p2*x^2 + p3*x + p4) /  (x^2 + q1*x + q2); % it is assumed only 80% is absorbed due to second time melting ... which the quality is not as first time
%         wasted_energy_portion=absorbtion_waste(2);
%     else
%         Abs_mandrel=@(x) absorbtion_waste(1);
%         wasted_energy_portion=absorbtion_waste(2);  % it means less goes to Tape
%     end
    
%% 19 Aug 2019
n_sub=absorbtion_waste(1);
n_tape=absorbtion_waste(2);

% Abs = FresnelAbs( n_tape, cosi)

    
  Abs=@(zaviye) FresnelAbs( n_tape, zaviye) ;
  
Abs_mandrel=@(zaviye) FresnelAbs( n_sub, zaviye) ; % it is assumed only 80% is absorbed due to second time melting ... which the quality is not as first time
        wasted_energy_portion=0;

        if length(absorbtion_waste)==3
    n_roller=absorbtion_waste(3);
     Abs_roller=@(zaviye) FresnelAbs( n_roller, zaviye) ;
     
else
    Abs_roller=Abs;
end


%%

% set(h100,'Position', [680 558 100 100 ]);


N_tape=Tape_Sp(1); % EVEN number of points, should not be changed !!
W_tape=Tape_Sp(2) ; % width of the tape
R_tape=Tape_Sp(3);
L_flat=Tape_Sp(4);
thick_T=Tape_Sp(5); % thickness of the tape
deg_tape=Tape_Sp(6);


% if Graphic_chekbox(7)
%     h=figure(1);
%     set(h,'Visible','off');
% end

% jobname='Example0';

% mkdir(strcat('.\Analysis_UOT\', jobname));
mkdir(UOT_pathfile);
steps=length(nip_point_M_all);
   tv_def_all=0*Roller_Pos_TV;
% path=5;
for ss=1:steps
    
    set(text_status,'String',[num2str((ss/steps)*100) '%']);
    drawnow;
    
        fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'w');
        fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'w');
        fileID3 = fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'w');
        
%     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'w');   % file includes the xyz + ID
%     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'w');     % file includes the xyz + ID
%     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'w');
    
    
    %%  instead of plotting in each loop, plot all afterwards
    
    fileID4 = fopen(strcat(UOT_pathfile, sprintf('int_xyz%d.txt',ss)),'w');
    fileID5 = fopen(strcat(UOT_pathfile, sprintf('Laser_Rays%d.txt',ss)),'w');
    fileID6 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Mandrel%d.txt',ss)),'w');
        fileID7 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Mandrel%d.txt',ss)),'w');
    fileID8 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Tape%d.txt',ss)),'w');
    fileID9 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Tape%d.txt',ss)),'w');
        fileID11 = fopen(strcat(UOT_pathfile, sprintf('Ref_int_xyz%d.txt',ss)),'w');
    
    
%     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'w');
%     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'w');
%     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'w');
%     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'w');
%     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'w');
%     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'w');
    
%     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'w');
    
    
    fprintf(fileID1,' X         ,Y         ,Z     , Energy  ,ID  \r\n');  % 0 means which reflection
    fprintf(fileID2,'X         ,Y         ,Z     , Energy   ,ID \r\n');  % 0 means which reflection
    fprintf(fileID3,'X         ,Y         ,Z     , Energy   ,ID \r\n');
    
    %ID = 0 intersection from laser, 1 from reflected rays
    
    % Define the input ray
    % Direction of laser head toward target
    Rx= delta_new_all(ss,1,1); %-0.8239; %input('direction of ray in x-direction =');
    Ry= delta_new_all(ss,1,2); %-0.1972;% input('direction of ray in y-direction =');
    Rz=delta_new_all(ss,1,3); %-0.3455; %input('direction of ray in z-direction =');
    normR=norm([Rx Ry Rz]); % normalize the Rx,Ry, Rz
    Rx=Rx/normR;
    Ry=Ry/normR;
    Rz=Rz/normR;
    
    
    
    %% plot roller
       Roller_Pos_TV=nip_point_M_all(ss,:,:);
    
    % plot tape .
%     H_indentation=0
    
    theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree
    
    Nip_Mov =R_tape*sind(theta_ind);
    
%        n_xyz=[0 0 0];

 

      % Modified on 23 April 2019 
 if Roller_Pos_TV(3) >=0 && Roller_Pos_TV(3) <=z_cyl_end
 
% gradient in cylender
n_xyz=[2*Roller_Pos_TV(1); 2*Roller_Pos_TV(2) ;0];
n_xyz=n_xyz/norm(n_xyz);

 elseif Roller_Pos_TV(3) < 0  
     % bottom dome part

     % gradient in cylender
n_xyz=[2*Roller_Pos_TV(1)/a^2; 2*Roller_Pos_TV(2)/b^2 ;2*Roller_Pos_TV(3)/c1^2];
n_xyz=n_xyz/norm(n_xyz);
 elseif Roller_Pos_TV(3) > z_cyl_end
     
     % For upper dome
     
     % gradient in cylender
n_xyz=[2*Roller_Pos_TV(1)/a^2; 2*Roller_Pos_TV(2)/b^2 ;2*(Roller_Pos_TV(3)-z_cyl_end)/c2^2];
n_xyz=n_xyz/norm(n_xyz);

        end


    tv_def=(n_xyz*H_indentation);
    %   -H_indentation+thick_T
    
    % What about thickness of the Tape ..???!!
    
 
    tv_def_all(ss,:,:)=tv_def;
 
    
    % It is already calculated in Kinematic model
    tv(1)=Roller_Pos_TV(1)+(-tv_def(1));
    tv(2)=Roller_Pos_TV(2)-tv_def(2);
    tv(3)=Roller_Pos_TV(3)-tv_def(3);
    
    
    % location of the Roller should be defined  and the gradient of the
    % surface mandrel
    %Thus:
    %  1 translation + 2 rotation should be performed
    % example
    
    
    % tv is the translation vector
    % tv=[R_cyl(1)/2;0 ; sqrt(3/4)*R_cyl(1) + z_cyl_end];  % for the Tape
    
    
%     if R_cyl (1) ~=0
%         [~,n_xyz]=Rot_Matrix_Finder_RollerTape(Roller_Pos_TV,th_y,R_cyl,z_cyl_end);
    
    
%     else
%     
%             R_cyl (1) ==0
%         Rotation around the x-axis
%         Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];
    
%         Rot_Roller_axis= [1 0 0; 0 cosd(90) -sind(90) ; 0 sind(90) cosd(90)]*Rot_y;
    
    
        % Rotation around the x-axis
        %     Rot_Roller_axis= [1 0 0; 0 cosd(90) -sind(90) ; 0 sind(90) cosd(90)]*Rot_y;
    
%         n_xyz=[0 ;0 ;-1];
%     end
    
%     [n1,n2,n3,x0,y0,z0,h_R,h_T,Center_roller]=Tape_plot_3D(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);

    
    xo_points=0+[0, R_tape, nan, 0, 0, nan, 0, 0];
yo_points=0+[0, 0, nan, 0, R_tape, nan, 0, 0];
zo_points=0+[0, 0, nan, 0, 0, nan, 0, W_tape];
   Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
   
   
    Temp=Rot_Roller_axis*[xo_points;yo_points;zo_points];
xo_points=Temp(1,:);
yo_points=Temp(2,:);
zo_points=Temp(3,:);

tv_x=tv(1);
tv_y=tv(2);
tv_z=tv(3);

% we already imposed the effect of Roller deformation
    xo_points=xo_points+tv_x;%-tv_def(1) ;   % Sth was wrong in probabely implementation that's why all shitty stupid change is done
yo_points=yo_points+tv_y;%+tv_def(2) ;
zo_points=zo_points+tv_z;%+tv_def(3) ;
    
    Center_roller=[xo_points(5), yo_points(5),zo_points(5) ];
    
    
    deg=ceil(deg_tape);

if deg_tape< theta_ind
    deg=floor(theta_ind);

end
ceil_theta_ind=floor(theta_ind);
    

R_thick=R_tape+thick_T;  

    [x,y,~] = cylinder(R_thick,N_tape); % for tape
    y=-y+R_thick;  % for tape
    x = x(:,N_tape/4+ceil_theta_ind:N_tape/4+deg);  % for continuous ends
% x(:,end) = x(:,1);
y = y(:,N_tape/4+ceil_theta_ind:N_tape/4+deg);
    
 x0=x(1,end);
y0=y(1,end);
z0=0;
    
    
    
    
    n1=- R_thick*sind(deg_tape);%   x(1,end); %input('direction of ray in x-direction =');
n2=- R_thick*cosd(deg_tape);% y(1,end)-R_thick;% input('direction of ray in y-direction =');
n3= 0.0; %input('direction of ray in z-direction =');
% normalize the Rx,Ry, Rz
normR=norm([n1 n2 n3]);
n1=n1/normR;
n2=n2/normR;
n3=n3/normR;



 
    %%
    
    
    
    % tv_R=Center_roller;
    
    tv_Tape_thickness=(n_xyz*(thick_T)); % SOME STUPID MANUPULATION IN THE NORMAL CALCULATIONS CAUSED TO HAVE A UNUSUAL NORMALS
    
    tv_R(1)=Center_roller(1)+0*tv_def(1)+tv_Tape_thickness(1);
    tv_R(2)=Center_roller(2)-0*tv_def(2)+tv_Tape_thickness(2);  % some speculation about + or - for tv_Tape_thickness
    tv_R(3)=Center_roller(3)-0*tv_def(3)+tv_Tape_thickness(3);   % some speculation about + or - for tv_Tape_thickness
    
    % plot3(tv_R(1),tv_R(2),tv_R(3),'c*','Markersize',8)
    
    
    %n1,n2,n3,x0,y0,z0 are the point and normal of intersection of curve and
    %flat part
    
    
    
    
    % counter_ray=0;
    
    % Gauss_Par_X=100;
    % Gauss_Par_Y=100;
    
    % ID=[2,Gauss_Par_X,Gauss_Par_Y]
    
    
    % if Graphic_chekbox(7)
    %     h=figure(1);
    %     set(h,'Visible','off');
    % end
    
    %
    % absorbtion_waste=0.8;
    % First_hit_energy=absorbtion_waste;  % should use better frensel equation or sth else
    % second_hit_energy=1-First_hit_energy;
    
    

    
    
    
    %% Specular MODEL
    % there should be a loop, in each time, Lx Ly LZ will be gotten
    %         from the laser_ray_generator
    
    No_graphic=1; % not to write reflections of A reflection
    
    % BRDF=BRDF_mode;
    
    
    
    X_Actual=LaserH_mesh(ss,1,:,:);
        Y_Actual=LaserH_mesh(ss,2,:,:);
        Z_Actual=LaserH_mesh(ss,3,:,:);
        
        
        LaserHead_C=reshape(pos_new_all(ss,:,:),[3,1]);
    
    if ~BRDF_mode
        for ii=1:counter_ray
          
            Lx=X_Actual(ii);
            Ly=Y_Actual(ii);% input('y position of Laser =');
            Lz=Z_Actual(ii);%  input('z position of Laser =');
            
            laser_source_P_G=[Lx;Ly;Lz];
            
            % Divergence
           Vec_C_point_actual=laser_source_P_G- LaserHead_C;
            
         b2=[Rx;Ry;Rz]+ (Vec_C_point_actual*Divergence_factor);


             b2=b2/norm(b2);
            
             
             
            % Define a direction for each point of laser head
            laser_direction_G=b2; %[Rx;Ry;Rz];
%             laser_direction_G=laser_point_Actual_all(ii,5:7)';
            
            
            xyz_int_G=[];
            part_index=0;
            
        
            
            % Check intersection with Mandrel
            % a function should be performed to do that
            
            
            %G_C >> General, Cylinder
            [xyz_int_G_C]=Line_Mandrel_intersection(laser_source_P_G,laser_direction_G,R_cyl,z_cyl_end);
            
            
            
            %input ray must be in local axis for tape calculations...
            % Working on the transformation matrix
            
            laser_source_P_L=Transformation_G2L (tv,Rot_Roller_axis,laser_source_P_G);
            laser_direction_L=Transformation_G2L (zeros(3,1),Rot_Roller_axis,laser_direction_G);
            
            [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(laser_source_P_L,laser_direction_L,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
            xyz_int_G_T=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
            
            
            %for the Roller, assume W_R=2*W
            laser_source_P_L_Roller=Transformation_G2L (tv_R,Rot_Roller_axis,laser_source_P_G);
            
            [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(laser_source_P_L_Roller,laser_direction_L,R_tape,W_R);
            xyz_int_G_Roller=Transformation_L2G (tv_R,Rot_Roller_axis,xyz_int_L_Roller);
            
            
            %%
            %     check which xyz_int are more closer to laser source
            
            dis_C=1e5;
            dis_T=1e5;
            dis_Roller=1e5;
            nm=0; % to check that there is no intersection with all part
            
            if ~isempty(xyz_int_G_T)
                dis_T=norm(laser_source_P_G-xyz_int_G_T);
            else
                nm=nm-1;
            end
            
            if ~isempty(xyz_int_G_C)
                dis_C=norm(laser_source_P_G-xyz_int_G_C);
            else
                nm=nm-1;
            end
            
            if ~isempty(xyz_int_G_Roller)
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
                [Ref_dir_G,beta]=refl_Mandrel_plot3D(R_cyl,xyz_int_G,laser_direction_G,z_cyl_end,fileID6,fileID7);
                part_index=1;
                
                First_hit_energy= Abs_mandrel(beta);
                
                %% For cheking the bump in an assumed distance >> 25 April 2019
%                 if xyz_int_G(3) >0.54 && xyz_int_G(3) <0.55 
%                      First_hit_energy= 1;
%                     wasted_energy_portion=0;
%                 end
                    
                    %%
                
                
                %                  fprintf(fileID1,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));  % 0 means which reflection
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),Power_Actual(ii)*First_hit_energy,ii);  % 0 means which reflection
                
                
                First_hit_energy=First_hit_energy+wasted_energy_portion;  % for reflection calculation recieve to the other parts
                
            elseif n==2
                %with Tape
                xyz_int_G=xyz_int_G_T;
                [Ref_dir_L,beta]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,laser_direction_L,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic); %input are in
                %     local axis system
                Ref_dir_G=Transformation_L2G (zeros(3,1),Rot_Roller_axis,Ref_dir_L);
                part_index=2;
                
                First_hit_energy= Abs(beta);
                
                %                 fprintf(fileID2,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),Power_Actual(ii)*First_hit_energy,ii);
                
                
            elseif n==3
                %with Roller
                % THE AMOUNT OF ABSORBTION FOR THE ROLLER SHOULD BE
                % DIFFERENT
                xyz_int_G=xyz_int_G_Roller;
                [Ref_dir_G,beta]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,laser_direction_L,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R);
                
                part_index=3;
                
                First_hit_energy= Abs_roller(beta);
                
                %                 fprintf(fileID3,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,Power_Actual(ii)*First_hit_energy,ii);
                
                
                
            end
            
            
            %    xyz_int_G(1),xyz_int_G(2),xyz_int_G(3)
            
            
            plot3D_line_interP(R_cyl,z_cyl_end,laser_direction_G,Lx,Ly,Lz,xyz_int_G,fileID4,fileID5);
            
            %% reflection
            
            
            if isempty(xyz_int_G)
                %     msgbox ('No intersection with all Parts!');
            elseif  part_index==2
                % continue the reflection calculation..
                xyz_int_G2=zeros(3,1);
                [xyz_int_G2]=Line_Mandrel_intersection(xyz_int_G,Ref_dir_G,R_cyl,z_cyl_end);
                %intersection points
                if xyz_int_G2
                    % plot3( xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3),'y.'  );
                    
                    fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                    
                    % text([xyz_int_G2(1)],[xyz_int_G2(2)],[xyz_int_G2(3)], ['P2']);
                    % intersection with cylinder, the reflected ray comes from tape
                    
                    [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G2,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                    
                    second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                    
                    fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), Power_Actual(ii)*second_hit_energy,ii);
                    
                end
                
                
            elseif  part_index==3   % maybe hit the tape or mandrel
                xyz_int_L2=Transformation_G2L (tv,Rot_Roller_axis,xyz_int_G); % checking for the tape
                Ref_dir_L2=Transformation_G2L (zeros(3,1),Rot_Roller_axis,Ref_dir_G);
                
                %for the tape
                [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
                xyz_int_G_T=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
                
                % for the mandrel
                [xyz_int_G_M]=Line_Mandrel_intersection(xyz_int_G,Ref_dir_G,R_cyl,z_cyl_end);
                
                
                
                
                
                
                if ~(isempty(xyz_int_G_M) || isempty(xyz_int_G_T) )  % intersection with 3 parts
                    
                    dis_T=norm(xyz_int_G-xyz_int_G_T);
                    dis_M=norm(xyz_int_G-xyz_int_G_M);
                    [m,n]=min([dis_T,dis_M]);
                    
                    
                    
                    if n==1
                        %with Tape
                        %                  plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                        
                        
                        
                        [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                        %input are in   local axis system
                        second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                        
                        fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                    elseif n==2
                        %with Mandrel
                        %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                        
                        
                        
                        [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                        
                        
                        second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                        
                        fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3), Power_Actual(ii)*second_hit_energy,ii);
                        
                    end
                    
                    
                elseif xyz_int_G_T
                    %       plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                    
                    
                    
                    %                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                    
                    [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                    %input are in   local axis system
                    second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                    
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),Power_Actual(ii)*second_hit_energy,ii);
                    
                    
                elseif xyz_int_G_M
                    %with Mandrel
                    %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                    
                    %                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                    
                    
                    %  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
                    
                    [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                    
                    second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                    
                    %
                    fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),Power_Actual(ii)*second_hit_energy,ii);
                    
                    
                    
                end
                
                
                
                
                
            elseif part_index==1
                
                xyz_int_L2=Transformation_G2L (tv,Rot_Roller_axis,xyz_int_G);
                xyz_int_L2_R=Transformation_G2L (tv_R,Rot_Roller_axis,xyz_int_G);
                Ref_dir_L2=Transformation_G2L (zeros(3,1),Rot_Roller_axis,Ref_dir_G);
                
                xyz_int_L=zeros(3,1);
                
                %for the tape
                [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
                xyz_int_G=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
                
                % for the roller
                [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(xyz_int_L2_R,Ref_dir_L2,R_tape,W_R);
                xyz_int_G_Roller=Transformation_L2G (tv_R,Rot_Roller_axis,xyz_int_L_Roller);
                
                
                
                
                
                if ~(isempty(xyz_int_G_Roller) || isempty(xyz_int_G) )  % intersection with 3 parts
                    
                    dis_T=norm(xyz_int_L2-xyz_int_L);
                    dis_Roller=norm(xyz_int_L2_R-xyz_int_L_Roller);
                    [m,n]=min([dis_T,dis_Roller]);
                    
                    
                    
                    if n==1
                        %with Tape
                        %                  plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    , it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                        
                        
                        
                        %                fprintf(fileID2,' %f %f %f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                        
                        
                        [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                        %input are in   local axis system
                        % should be modified >> how much loss to ambient !!
                        second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                        
                        
                        
                        fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                    elseif n==2
                        %with Roller
                        %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                        
                        
                        %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                        
                        %input are in   local axis system
                        
                        
                        
                        [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                        % should be modified >> how much loss to ambient !!
                        second_hit_energy=(1-First_hit_energy)*Abs_roller(beta2);
                        
                        
                        fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3), Power_Actual(ii)*second_hit_energy,ii);
                        
                    end
                    
                elseif xyz_int_L
                    %          plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                    
                    % intersection with tape, the reflected ray comes from cylinder
                    %      fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                    
                    [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                    
                    
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), Power_Actual(ii)*second_hit_energy,ii);
                    
                    
                elseif xyz_int_G_Roller
                    
                    %with Roller
                    %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    
                    % only reflections
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    
                    %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=(1-First_hit_energy)*Abs_roller(beta2);
                    
                    
                    fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),Power_Actual(ii)*second_hit_energy,ii);
                    
                    
                end
                
                
            end
            
            
            
            
            
        end
        
        
    else
        %% Non-specular model- BRDF
        
        
        
        fileID_BRDF_par = fopen('.\Supp_files\BRDF_Par_sub_tape.txt','r');
        
        out = textscan(fileID_BRDF_par,'%*s %f %f','delimiter',',');
        
        BRDF_sub=out{1};
        BRDF_Tape=out{2};
        
        fclose( fileID_BRDF_par);
        
        for ii=1:counter_ray
          Lx=X_Actual(ii);
            Ly=Y_Actual(ii);% input('y position of Laser =');
            Lz=Z_Actual(ii);%  input('z position of Laser =');
            
            laser_source_P_G=[Lx;Ly;Lz];
            
            % Define a direction for each point of laser head
            %  laser_direction_G=[Rx;Ry;Rz];
%             laser_direction_G=laser_point_Actual_all(ii,5:7)';
             laser_direction_G=[Rx;Ry;Rz];
            
            xyz_int_G=[];
            part_index=0;
            
           
            % Check intersection with Mandrel
            % a function should be performed to do that
            
            
            %G_C >> General, Cylinder
            [xyz_int_G_C]=Line_Mandrel_intersection(laser_source_P_G,laser_direction_G,R_cyl,z_cyl_end);
            
            
            
            %input ray must be in local axis for tape calculations...
            % Working on the transformation matrix
            
            laser_source_P_L=Transformation_G2L (tv,Rot_Roller_axis,laser_source_P_G);
            laser_direction_L=Transformation_G2L (zeros(3,1),Rot_Roller_axis,laser_direction_G);
            
            [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(laser_source_P_L,laser_direction_L,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
            xyz_int_G_T=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
            
            
            %for the Roller, assume W_R=2*W
            laser_source_P_L_Roller=Transformation_G2L (tv_R,Rot_Roller_axis,laser_source_P_G);
            
            [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(laser_source_P_L_Roller,laser_direction_L,R_tape,W_R);
            xyz_int_G_Roller=Transformation_L2G (tv_R,Rot_Roller_axis,xyz_int_L_Roller);
            
            
            %%
            %     check which xyz_int are more closer to laser source
            
            dis_C=1e5;
            dis_T=1e5;
            dis_Roller=1e5;
            nm=0; % to check that there is no intersection with all part
            
            if ~isempty(xyz_int_G_T)
                dis_T=norm(laser_source_P_G-xyz_int_G_T);
            else
                nm=nm-1;
            end
            
            if ~isempty(xyz_int_G_C)
                dis_C=norm(laser_source_P_G-xyz_int_G_C);
            else
                nm=nm-1;
            end
            
            if ~isempty(xyz_int_G_Roller)
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
                [Ref_dir_G,beta]=refl_Mandrel_plot3D_BRDF(R_cyl,xyz_int_G,laser_direction_G,z_cyl_end,fileID6,fileID7,BRDF_sub);
                part_index=1;
                
                First_hit_energy= Abs_mandrel(beta);
                
                %                  fprintf(fileID1,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));  % 0 means which reflection
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),Power_Actual(ii)*First_hit_energy,ii);  % 0 means which reflection
                
                
                First_hit_energy=First_hit_energy+wasted_energy_portion;  % for reflection calculation recieve to the other parts
                
            elseif n==2
                %with Tape
                xyz_int_G=xyz_int_G_T;
                [Ref_dir_G,beta]=tape_ray_reflection_3D_BRDF(R_tape,S1,S2,xyz_int_L,x0,y0,laser_direction_L,tv,Rot_Roller_axis,fileID8,fileID9,xyz_int_G,BRDF_Tape); %input are in
                %     local axis system
                %             Ref_dir_G=Transformation_L2G (zeros(3,1),Rot_Roller_axis,Ref_dir_L);
                part_index=2;
                
                % specular assumption, only one ray reflected
                %             Ref_dir_G(4,:)=1;
                
                First_hit_energy= Abs(beta);
                
                %                 fprintf(fileID2,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),Power_Actual(ii)*First_hit_energy,ii);
                
                
            elseif n==3
                %with Roller
                % THE AMOUNT OF ABSORBTION FOR THE ROLLER SHOULD BE
                % DIFFERENT
                xyz_int_G=xyz_int_G_Roller;
                [Ref_dir_G,beta]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,laser_direction_L,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R);
                
                % specular assumption, only one ray reflected
                Ref_dir_G(4,:)=1;
                
                part_index=3;
                
                First_hit_energy= Abs_roller(beta);
                
                %                 fprintf(fileID3,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,Power_Actual(ii)*First_hit_energy,ii);
                
                
                
            end
            
            
            %    xyz_int_G(1),xyz_int_G(2),xyz_int_G(3)
            
            
            plot3D_line_interP(R_cyl,z_cyl_end,laser_direction_G,Lx,Ly,Lz,xyz_int_G,fileID4,fileID5);
            
            %% reflection continues for number of Active microfacet reflection
            [Row,Col]=size(Ref_dir_G) ;
            
            %       Ref_dir_G should structure of 4*(Number of active rays)
            
            
            for CC=1:Col
                
                if isempty(xyz_int_G)
                    %     msgbox ('No intersection with all Parts!');
                elseif  part_index==2   % from Tape
                    % continue the reflection calculation..
                    xyz_int_G2=zeros(3,1);
                    [xyz_int_G2]=Line_Mandrel_intersection(xyz_int_G,Ref_dir_G(1:3,CC),R_cyl,z_cyl_end);
                    %intersection points
                    if xyz_int_G2
                        % plot3( xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3),'y.'  );
                        
                        fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                        
                        % text([xyz_int_G2(1)],[xyz_int_G2(2)],[xyz_int_G2(3)], ['P2']);
                        % intersection with cylinder, the reflected ray comes from tape
                        
                        [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G2,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                        
                        
                        second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                        
                        fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), Power_Actual(ii)*second_hit_energy,ii);
                        
                    end
                    
                    
                elseif  part_index==3   % maybe hit the tape or mandrel
                    xyz_int_L2=Transformation_G2L (tv,Rot_Roller_axis,xyz_int_G); % checking for the tape
                    Ref_dir_L2=Transformation_G2L (zeros(3,1),Rot_Roller_axis,Ref_dir_G(1:3,CC));
                    
                    %for the tape
                    [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
                    xyz_int_G_T=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
                    
                    % for the mandrel
                    [xyz_int_G_M]=Line_Mandrel_intersection(xyz_int_G,Ref_dir_G(1:3,CC),R_cyl,z_cyl_end);
                    
                    
                    
                    
                    
                    
                    if ~(isempty(xyz_int_G_M) || isempty(xyz_int_G_T) )  % intersection with 3 parts
                        
                        dis_T=norm(xyz_int_G-xyz_int_G_T);
                        dis_M=norm(xyz_int_G-xyz_int_G_M);
                        [m,n]=min([dis_T,dis_M]);
                        
                        
                        
                        if n==1
                            %with Tape
                            %                  plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                            fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                            
                            
                            
                            
                            [~,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                            %input are in   local axis system
                            second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                            
                            fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),Power_Actual(ii)*second_hit_energy,ii);
                            
                            
                        elseif n==2
                            %with Mandrel
                            %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                            fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                            
                            
                            
                            [~,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                            
                            
                            second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                            
                            fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3), Power_Actual(ii)*second_hit_energy,ii);
                            
                        end
                        
                        
                    elseif xyz_int_G_T
                        %       plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                        
                        
                        
                        %                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                        
                        [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                        %input are in   local axis system
                        second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                        
                        
                        fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                    elseif xyz_int_G_M
                        %with Mandrel
                        %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                        
                        %                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                        
                        
                        %  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
                        
                        [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                        
                        second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                        
                        %
                        fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                        
                    end
                    
                    
                    
                    
                    
                elseif part_index==1
                    
                    xyz_int_L2=Transformation_G2L (tv,Rot_Roller_axis,xyz_int_G);
                    xyz_int_L2_R=Transformation_G2L (tv_R,Rot_Roller_axis,xyz_int_G);
                    Ref_dir_L2=Transformation_G2L (zeros(3,1),Rot_Roller_axis,Ref_dir_G(1:3,CC));
                    
                    xyz_int_L=zeros(3,1);
                    
                    %for the tape
                    [xyz_int_L,S1,S2]=intersection_tape_Ray_3D(xyz_int_L2,Ref_dir_L2,R_tape+thick_T,W_tape,deg_tape,n1,n2,n3,x0,y0,z0,L_flat,W_R,theta_ind);
                    xyz_int_G2=Transformation_L2G (tv,Rot_Roller_axis,xyz_int_L);
                    
                    % for the roller
                    [xyz_int_L_Roller]=Line_cylinder_intersection_Roller(xyz_int_L2_R,Ref_dir_L2,R_tape,W_R);
                    xyz_int_G_Roller=Transformation_L2G (tv_R,Rot_Roller_axis,xyz_int_L_Roller);
                    
                    
                    
                    
                    
                    if ~(isempty(xyz_int_G_Roller) || isempty(xyz_int_G2) )  % intersection with 3 parts
                        
                        dis_T=norm(xyz_int_L2-xyz_int_L);
                        dis_Roller=norm(xyz_int_L2_R-xyz_int_L_Roller);
                        [m,n]=min([dis_T,dis_Roller]);
                        
                        
                      
                        if n==1
                            %with Tape
                            %                  plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                            fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                            
                            
                            
                            %                fprintf(fileID2,' %f %f %f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                            
                            
                            [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                            %input are in   local axis system
                            %% should be modified >> how much loss to ambient !!
                            second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                            
                            
                            
                            fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), Power_Actual(ii)*second_hit_energy,ii);
                            
                            
                        elseif n==2
                            %with Roller
                            %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                            fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                            
                            
                            %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                            
                            %input are in   local axis system
                            
                            
                            
                            [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                            %% should be modified >> how much loss to ambient !!
                            second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_roller(beta2);
                            
                            
                            fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3), Power_Actual(ii)*second_hit_energy,ii);
                            
                        end
                        
                    elseif xyz_int_L
                        %          plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                        
                        % intersection with tape, the reflected ray comes from cylinder
                        %      fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                        
                        [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                        %% should be modified >> how much loss to ambient !!
                        second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                        
                        
                        
                        fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f %d \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                    elseif xyz_int_G_Roller
                        
                        %with Roller
                        %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                        
                        % only reflections
                        fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                        
                        
                        %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                        
                        [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                        %% should be modified >> how much loss to ambient !!
                        second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_roller(beta2);
                        
                        
                        fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f %d \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),Power_Actual(ii)*second_hit_energy,ii);
                        
                        
                    end
                    
                    
                end
                
                
                
                
                
            end
            
            
            
        end
        
        
        
        
        
    end
    
    
    
    
    
   fclose(fileID1);
fclose(fileID2);
fclose(fileID3); 
    
       fclose(fileID4);
fclose(fileID5);
fclose(fileID6); 
           fclose(fileID7);
fclose(fileID8);
fclose(fileID9); 


%        fclose(fileID10);
fclose(fileID11);
 
end












% grid  on;
% % axis([-2*R_cyl 2*R_cyl -2*R_cyl 2*R_cyl -z_cyl_end/2 1.5*z_cyl_end]);
% view([180 0]);




%first intersection are black stars
% legend('yellow for second intersections')
% reflection from cylinder are red
% reflection from the tape are blue



% toc



  
           
           
%     f = figure;
% c = uicontrol(h,'Style','popupmenu');
% c.Position = [0 100 90 40];
% c.String = {1:length(nip_point_M_all)};
% c.Callback = @show_step_outputs;
    
    
%         btn8 = uicontrol('Style', 'popupmenu', 'String', 'Clear Normal_lines',...
%         'Position', [900 20 90 40],...
%         'Callback', 'delete(Normal_lines)');
    



% The variables are specifically belongs to this analysis with these
% parameters,
% if a user change parameters, the computed analysis is not changed > that
% is why we should save these variables
  


%      save(strcat(dir,'\Optical_UOT.mat'),'Rot_Roller_axis_all','nip_point_M_all',...
%          'CV_mesh','LaserH_mesh','delta_new_all','pos_new_all');

% jobname
% ss

 if ~BRDF_mode
dir2save=UOT_pathfile;
     save(strcat(dir2save,'\Optical_UOT.mat'),'R_cyl','z_cyl_end','nip_point_M_all','Rot_Roller_axis_all',...
'N_tape', 'W_tape','R_tape','L_flat','tv','th_y','thick_T','deg_tape','W_R','theta_ind','H_indentation','tv_def_all',...
 'CV_mesh','jobname','Laser_head');

 else
    dir2save=UOT_pathfile;
     save(strcat(dir2save,'\Optical_BRDF_UOT.mat'),'R_cyl','z_cyl_end','nip_point_M_all','Rot_Roller_axis_all',...
'N_tape', 'W_tape','R_tape','L_flat','tv','th_y','thick_T','deg_tape','W_R','theta_ind','H_indentation','tv_def_all',...
 'CV_mesh','jobname','Laser_head');
     
 end





% 
% 
%    function show_step_outputs(src,jobname,nip_point_M_all,...
%                N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def,...
%                      rb1,rb2,rb3,rb4,rb6,rb7,rb8)
%         val = src.Value;
%         str = src.String;
% %         str{val};
% %         disp(['Selection: ' str{val}]);
%         
% Fig_handles = guidata(gcbo);
%         
%       ss=val;
% 
% % jobname='Example0';
%     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'r');   % file includes the xyz + ID
%     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'r');     % file includes the xyz + ID
%     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'r');
%     
%     
%     %%  instead of plotting in each loop, plot all afterwards
%     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'r');
%     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'r');
%     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'r');
%     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'r');
%     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'r');
%     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'r');
%     
%     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'r');  
%         
%        
% h=figure(1);
% 
% 
% 
% % delete([h_R])
% %  delete([Fig_handles.h_R',Fig_handles.h_T,Fig_handles.Graphics_laser,Fig_handles.Graphics_reflection',Fig_handles.reflection_int,Fig_handles.Laser_int_points,Fig_handles.Normal_lines']);
% % plot Mandrel
% % First is cylinder , second is bottom dome part, third is upper Dome part
% % Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);
% % assignin('base','Mandrel_plot',Mandrel_plot);
% 
% 
% % plot3(tv(1),tv(2),tv(3),'g*','Markersize',5)
% 
% 
%    Roller_Pos_TV=nip_point_M_all(ss,:,:);
%     
%     
%     tv(1)=Roller_Pos_TV(1)-(-tv_def(1));
%     tv(2)=Roller_Pos_TV(2)-tv_def(2);
%     tv(3)=Roller_Pos_TV(3)+tv_def(3);
%     
%      Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
% % output is real center of roller without deformation, end of arrow of
% % the local axis
%   [n1,n2,n3,x0,y0,z0,Fig_handles.h_R,Fig_handles.h_T,Center_roller]=Tape_plot_3D_Tshow(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);
% 
% assignin('base','h_R',Fig_handles.h_R);
% assignin('base','h_T',Fig_handles.h_T); 
%         
%         
%         
%            XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
%     XYZ_int=cell2mat(XYZ_int);
%     Fig_handles.Laser_int_points=plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
%     
%     XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     XYZ_Rays=cell2mat(XYZ_Rays);
%    Fig_handles.Graphics_laser= plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',2 );
%     
%     Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
%    Fig_handles.Normal_lines1=plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:');
%     
%     
%     Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
%     Fig_handles.Graphics_reflection1=plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--');
%     
%     
%     Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
%    Fig_handles.Normal_lines2= plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
%     
%     
%     Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
%     Fig_handles.Graphics_reflection2=plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--');
%     
%     Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Ref_int_xyz=cell2mat(Ref_int_xyz);
%     Fig_handles.reflection_int=plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');
%     
%     
%     %% Changing the graphical object, Visibility or deleting
%     
% %     Graphics_laser=findobj(h,'color','g','LineStyle',':');
%     assignin('base','Graphics_laser',Fig_handles.Graphics_laser);
%     % set(Graphics_laser,'Visible','off');
% %     Graphics_reflection1=findobj(h,'color','r');
% %     Graphics_reflection2=findobj(h,'color','b');
%     Fig_handles.Graphics_reflection=[Fig_handles.Graphics_reflection1;Fig_handles.Graphics_reflection2];
%     assignin('base','Graphics_reflection',Fig_handles.Graphics_reflection);
%     
% %     reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
%     assignin('base','reflection_int',Fig_handles.reflection_int);
%     
% %     Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
%     assignin('base','Laser_int_points',Fig_handles.Laser_int_points);
%     
% %     Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
% Fig_handles.Normal_lines=[Fig_handles.Normal_lines1;Fig_handles.Normal_lines2];
%     assignin('base','Normal_lines',Fig_handles.Normal_lines);
%     
%     % Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
%     % assignin('base','Thermal_points',Thermal_points);
%     
%     
%     
%     
%     
%     
%     % set(Graphics_reflection,'Visible','off');
%     
%     
%     
% 
%     set(rb1,'Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Fig_handles.Graphics_laser));
%     set(rb2,'Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Fig_handles.Graphics_reflection));
%     % mandrel is not updating
%      %set(rb3,'Callback', @(rb3,event) Visible_plotButtonPushed(rb1,Graphics_laser));
%     set(rb4,'Callback', @(rb4,event) Visible_plotButtonPushed(rb4,Fig_handles.h_R));
%     
%     
% %     rb3 = uicontrol('Style', 'radiobutton','Position',[280 50 90 40],'Value',true,...
% %         'String','Mandrel','Callback', @(rb3,event) Visible_plotButtonPushed(rb3,Mandrel_plot));
%     
% %     btn3 = uicontrol('Style', 'pushbutton', 'String', 'Clear Mandrel',...
% %         'Position', [280 20 90 40],...
% %         'Callback', 'delete(Mandrel_plot)');
%     
%     
% %     rb4 = uicontrol('Style', 'radiobutton','Position',[400 50 90 40],'Value',true,...
% %         'String','Roller','Callback', @(rb4,event) Visible_plotButtonPushed(rb4,h_R));
%     
% %     btn4 = uicontrol('Style', 'pushbutton', 'String', 'Clear Roller',...
% %         'Position', [400 20 90 40],...
% %         'Callback', 'delete(h_R)');
%     
% 
%     set(rb6,'Callback', @(rb6,event) Visible_plotButtonPushed(rb6,Fig_handles.reflection_int));
%     set(rb7,'Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Fig_handles.Laser_int_points));
%     set(rb8,'Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Fig_handles.Normal_lines));
%     
% %     rb6 = uicontrol('Style', 'radiobutton','Position',[640 50 90 40],'Value',true,...
% %         'String','1st reflections','Callback', @(rb6,event) Visible_plotButtonPushed(rb6,reflection_int));
% % 
% %     
% %     
% %     rb7 = uicontrol('Style', 'radiobutton','Position',[760 50 90 40],'Value',true,...
% %         'String','Laser_int','Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Laser_int_points));
%     
% 
%     
% %     rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
% %         'String','Normal lines','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Normal_lines));
%     
%     
% %     string=[]
% %     btn8 = uicontrol('Style', 'pushbutton', 'String', 'Clear Normal_lines',...
% %         'Position', [900 20 90 40],...
% %         'Callback', 'delete(Normal_lines)'); 
%         
%         
%         guidata(gcbo,Fig_handles) ;
%         
%    end


