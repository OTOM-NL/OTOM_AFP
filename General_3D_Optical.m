
% This program is able to calculate and combine of two axis systems to
% calculate intersectioin points and 1st reflection
%Input: Laser ray point and direction 3D, geometrical paratmers

%Transformation of rays and intersection points should be into axis-1 which
%is considered as a general axis system
% For caluclation in axis-2 all the parameters should be transformed into
% axis-2


% Calculate for 1st reflection


function [counter_ray,Rot_Roller_axis,tv]=General_3D_Optical(th_y,Tape_Sp,...
    Rxyz,R_cyl,z_cyl_end,Roller_Pos_TV,ID,Laser_head,L_xyz0,absorbtion_waste,W_R,H_indentation,...
    Graphic_chekbox,Laser_head_Rot,...
    BRDF_mode,Divergence_factor)

%   tic

% tv3 AND absorbtion_waste are useless here


if Graphic_chekbox(7)
    h1=figure(1);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h1,'Position', [0 0 100 100 ]);
    cameratoolbar(h1,'show');
end


if Graphic_chekbox(8)
    h2=figure(2);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h2,'Position', [0 0  100 100 ]);
end

if Graphic_chekbox(9)
    h3=figure(3);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h3,'Position', [0 0  100 100 ]);
    
    h100=figure(100);
    set(h100,'Position', [0 0 100 100 ]);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    
end


if  Graphic_chekbox(3) &&  Graphic_chekbox(5)
    h21=figure(21);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h21,'Position', [0 0  100 100 ]);
end








% set(h100,'Position', [680 558 100 100 ]);


N_tape=Tape_Sp(1); % EVEN number of points, should not be changed !!
W_tape=Tape_Sp(2) ; % width of the tape
R_tape=Tape_Sp(3);
L_flat=Tape_Sp(4);
thick_T=Tape_Sp(5); % thickness of the tape
deg_tape=Tape_Sp(6);


if Graphic_chekbox(7)
    h=figure(1);
    set(h,'Visible','off');
end


fileID1 = fopen('.\Cylinder_ints.txt','w');   % file includes the xyz + ID
fileID2 = fopen('.\Tape_ints.txt','w');     % file includes the xyz + ID
fileID3 = fopen('.\Roller_ints.txt','w');


%%  instead of plotting in each loop, plot all afterwards


fileID18 = fopen('.\beta_Tape.txt','w');
fileID19 = fopen('.\beta_Sub.txt','w');

fileID28 = fopen('.\beta2_Tape.txt','w');
fileID29 = fopen('.\beta2_Sub.txt','w');
%%

fileID4 = fopen('.\int_xyz.txt','w');
fileID5 = fopen('.\Laser_Rays.txt','w');
fileID6 = fopen('.\Normal_vectors_Mandrel.txt','w');
fileID7 = fopen('.\Reflection_vector_Mandrel.txt','w');
fileID8 = fopen('.\Normal_vectors_Tape.txt','w');
fileID9 = fopen('.\Reflection_vector_Tape.txt','w');

fileID11 = fopen('.\Ref_int_xyz.txt','w');


fprintf(fileID1,' X         ,Y         ,Z       ,ID  \r\n');  % 0 means which reflection
fprintf(fileID2,'X         ,Y         ,Z       ,ID \r\n');  % 0 means which reflection
fprintf(fileID3,'X         ,Y         ,Z       ,ID \r\n');

fprintf(fileID18,'Ray      \r\n');
fprintf(fileID19,'Ray      \r\n');
fprintf(fileID28,'Ref     \r\n');
fprintf(fileID29,'Ref    \r\n');
%ID = 0 intersection from laser, 1 from reflected rays

% Define the input ray
% Direction of laser head toward target
Rx= Rxyz(1); %-0.8239; %input('direction of ray in x-direction =');
Ry= Rxyz(2); %-0.1972;% input('direction of ray in y-direction =');
Rz=Rxyz(3); %-0.3455; %input('direction of ray in z-direction =');
normR=norm([Rx Ry Rz]); % normalize the Rx,Ry, Rz
Rx=Rx/normR;
Ry=Ry/normR;
Rz=Rz/normR;



%     laser_direction_G=[Rx;Ry;Rz];
% geometrical of cylinder

% plot Mandrel
% First is cylinder , second is bottom dome part, third is upper Dome part
if Graphic_chekbox(7)
Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);
assignin('base','Mandrel_plot',Mandrel_plot);
end
%% plot roller


% plot tape .

%tv  > translation vector
% tv=[0;R_cyl;z_cyl_end/2];  % for the tape
% tv_R=[0;R_cyl+R_tape+thick_T;z_cyl_end/2];  % for center of Roller


theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree

Nip_Mov =R_tape*sind(theta_ind);



% location of the Roller should be defined  and the gradient of the
% surface mandrel
%Thus:
%  1 translation + 2 rotation should be performed
% example


% tv is the translation vector
% tv=[R_cyl(1)/2;0 ; sqrt(3/4)*R_cyl(1) + z_cyl_end];  % for the Tape


if R_cyl (1) ~=0
    [Rot_Roller_axis,n_xyz]=Rot_Matrix_Finder_RollerTape(Roller_Pos_TV,th_y,R_cyl,z_cyl_end);
    
    
else
    
    %     R_cyl (1) ==0
    % Rotation around the x-axis
    Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];
    
    Rot_Roller_axis= [1 0 0; 0 cosd(90) -sind(90) ; 0 sind(90) cosd(90)]*Rot_y;
    
    
    % Rotation around the x-axis
    %     Rot_Roller_axis= [1 0 0; 0 cosd(90) -sind(90) ; 0 sind(90) cosd(90)]*Rot_y;
    
    % it was 0 0 -1  >> 17 juli 2019
    n_xyz=[0 ;0 ;1];
end



tv_def=(n_xyz*H_indentation);
%   -H_indentation+thick_T

% What about thickness of the Tape ..???!!

tv(1)=Roller_Pos_TV(1)+(-tv_def(1));
tv(2)=Roller_Pos_TV(2)-tv_def(2);
tv(3)=Roller_Pos_TV(3)-tv_def(3);
%%

% plot3(tv(1),tv(2),tv(3),'g*','Markersize',5)

% output is real center of roller without deformation, end of arrow of
% the local axis

[n1,n2,n3,x0,y0,z0,h_R,h_T,Center_roller]=Tape_plot_3D(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);
assignin('base','h_R',h_R);
assignin('base','h_T',h_T);

% tv_R=Center_roller;

tv_Tape_thickness=(n_xyz*(thick_T)); % SOME STUPID MANUPULATION IN THE NORMAL CALCULATIONS CAUSED TO HAVE A UNUSUAL NORMALS


% modified on 7/3/2019
% if R_cyl (1) ==0
% tv_Tape_thickness=(n_xyz*(thick_T)); % SOME STUPID MANUPULATION IN THE NORMAL CALCULATIONS CAUSED TO HAVE A UNUSUAL NORMALS
% 
%     
% end

% Before it was -tv_Tape_thickness(1) >>> modified on 16 April 2019

% it was  changed on 7/3/2019

tv_R(1)=Center_roller(1)-0*tv_def(1)+tv_Tape_thickness(1);
tv_R(2)=Center_roller(2)-0*tv_def(2)+tv_Tape_thickness(2);  % some speculation about + or - for tv_Tape_thickness
tv_R(3)=Center_roller(3)-0*tv_def(3)+tv_Tape_thickness(3);   % some speculation about + or - for tv_Tape_thickness

% plot3(tv_R(1),tv_R(2),tv_R(3),'c*','Markersize',8)


%n1,n2,n3,x0,y0,z0 are the point and normal of intersection of curve and
%flat part

axis equal;

% an ID which define the laser distribution

Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);

Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);



counter_ray=Laser_head_nx*Laser_head_ny;

% counter_ray=0;

% Gauss_Par_X=100;
% Gauss_Par_Y=100;

% ID=[2,Gauss_Par_X,Gauss_Par_Y]

laser_point_Actual_all=Laser_ray_generator (Rx,Ry,Rz,Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,L_xyz0,ID,Laser_head_Rot,Divergence_factor);


if Graphic_chekbox(7)
    h=figure(1);
    set(h,'Visible','off');
end

%
% absorbtion_waste=0.8;
% First_hit_energy=absorbtion_waste;  % should use better frensel equation or sth else
% second_hit_energy=1-First_hit_energy;


%% >> Modified on 11 Aug 2018  absorption fault
% p1 =   -6.55e-08 ;
% p2 =    1.56e-05  ;
% p3 =   -0.001348 ;
% p4 =     0.05033  ;
% p5 =       0.198 ;
% 
% Abs=@(x) p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;

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
    
%  Abs = @(x)  (p1*x^3 + p2*x^2 + p3*x + p4) /  (x^2 + q1*x + q2);
    
    
    %% Should be modified
%     if absorbtion_waste(1) ==-1
%         Abs_mandrel=@(x) (p1*x^3 + p2*x^2 + p3*x + p4) /  (x^2 + q1*x + q2); % it is assumed only 80% is absorbed due to second time melting ... which the quality is not as first time
%         wasted_energy_portion=absorbtion_waste(2);
%     else
%         Abs_mandrel=@(x) absorbtion_waste(1);
%         wasted_energy_portion=absorbtion_waste(2);  % it means less goes to Tape
%     end
%     

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


%% Specular MODEL
% there should be a loop, in each time, Lx Ly LZ will be gotten
%         from the laser_ray_generator

    No_graphic=1; % not to write reflections of A reflection
      
% BRDF=BRDF_mode;


% poolobj = gcp('nocreate');
% 
% if isempty(poolobj)
%     poolsize = 0;
%     parpool('local',7);
% else
%     poolsize = poolobj.NumWorkers
% end

% parpool('local',6);


if ~BRDF_mode
    for ii=1:counter_ray
        Lx=laser_point_Actual_all(ii,1);
        Ly=laser_point_Actual_all(ii,2);% input('y position of Laser =');
        Lz=laser_point_Actual_all(ii,3);%  input('z position of Laser =');
        
        laser_source_P_G=[Lx;Ly;Lz];
        
        
        
        % Define a direction for each point of laser head
        %  laser_direction_G=[Rx;Ry;Rz];
        laser_direction_G=laser_point_Actual_all(ii,5:7)';
        
        
        xyz_int_G=[];
        part_index=0;
        
        %     xyz_int_L=zeros(3,1);
        %     xyz_int_G_C=zeros(3,1);
        %     xyz_int_G_T=zeros(3,1);
        % xyz_int_G=[] ;%zeros(3,1);
        %     xyz_int_G_Roller=zeros(3,1);
        %     xyz_int_L_Roller=zeros(3,1);
        
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
            
            %                  fprintf(fileID1,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));  % 0 means which reflection
            fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);  % 0 means which reflection
            
            
            First_hit_energy=First_hit_energy+wasted_energy_portion;  % for reflection calculation recieve to the other parts
            
            % for paper 3
             fprintf(fileID19,' %12.8f %12.8f \r\n',beta,First_hit_energy);  % 0 means which reflection
            
            
            
            
        elseif n==2
            %with Tape
            xyz_int_G=xyz_int_G_T;
            [Ref_dir_L,beta]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,laser_direction_L,tv,Rot_Roller_axis,fileID8,fileID9); %input are in
            %     local axis system
            Ref_dir_G=Transformation_L2G (zeros(3,1),Rot_Roller_axis,Ref_dir_L);
            part_index=2;
            
            First_hit_energy= Abs(beta);
            
            %                 fprintf(fileID2,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
            fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);
            
            % for paper 3
             fprintf(fileID18,' %12.8f %12.8f \r\n',beta,First_hit_energy);  % 0 means which reflection
            
            
        elseif n==3
            %with Roller
            % THE AMOUNT OF ABSORBTION FOR THE ROLLER SHOULD BE
            % DIFFERENT
            xyz_int_G=xyz_int_G_Roller;
            [Ref_dir_G,beta]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,laser_direction_L,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R);
            
            part_index=3;
            
            First_hit_energy= Abs_roller(beta);
            
            %                 fprintf(fileID3,' %f %f %f , 0 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
            fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,laser_point_Actual_all(ii,4)*First_hit_energy);
            
            
            
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
            if ~isempty(xyz_int_G2)
                % plot3( xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3),'y.'  );
                
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                
                % text([xyz_int_G2(1)],[xyz_int_G2(2)],[xyz_int_G2(3)], ['P2']);
                % intersection with cylinder, the reflected ray comes from tape
                
                [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G2,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                
                second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                
                     fprintf(fileID29,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
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
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                     fprintf(fileID28,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                elseif n==2
                    %with Mandrel
                    %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                    
                    
                    
                    [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                    
                     fprintf(fileID29,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                    second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                    
                    fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                end
                
                
            elseif ~isempty(xyz_int_G_T) 
                %       plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                
                
                
                %                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                
                [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                %input are in   local axis system
                second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                 fprintf(fileID28,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                
                fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
            elseif ~isempty(xyz_int_G_M) 
                %with Mandrel
                %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                
                %                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                
                
                %  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
                
                [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G,z_cyl_end,fileID6,fileID7,No_graphic);
                
                second_hit_energy=(1-First_hit_energy)*Abs_mandrel(beta2);
                
                 fprintf(fileID29,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                %
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
                
            end
            
            
            
            
            
        elseif part_index==1;
            
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
                
                
                %%
                if n==1
                    %with Tape
                    %                  plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                    
                    
                    
                    %                fprintf(fileID2,' %f %f %f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                    
                    
                    [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                    %input are in   local axis system
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                    
                     fprintf(fileID28,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                    
                elseif n==2
                    %with Roller
                    %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    
                    %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    %input are in   local axis system
                    
                    
                    
                    [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=(1-First_hit_energy)*Abs_roller(beta2);
                    
                    
                    fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                end
                
            elseif ~isempty(xyz_int_L) 
                %          plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                
                % intersection with tape, the reflected ray comes from cylinder
                %      fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                
                [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                %% should be modified >> how much loss to ambient !!
                second_hit_energy=(1-First_hit_energy)*Abs(beta2);
                
                 fprintf(fileID28,' %12.8f %12.8f \r\n',beta2,second_hit_energy);  % 0 means which reflection
                     
                
                fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
            elseif ~isempty(xyz_int_G_Roller) 
                
                %with Roller
                %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                
                % only reflections
                fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                
                
                %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                
                [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                %% should be modified >> how much loss to ambient !!
                second_hit_energy=(1-First_hit_energy)*Abs_roller(beta2);
                
                
                fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
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
        Lx=laser_point_Actual_all(ii,1);
        Ly=laser_point_Actual_all(ii,2);% input('y position of Laser =');
        Lz=laser_point_Actual_all(ii,3);%  input('z position of Laser =');
        
        laser_source_P_G=[Lx;Ly;Lz];
        
        % Define a direction for each point of laser head
        %  laser_direction_G=[Rx;Ry;Rz];
        laser_direction_G=laser_point_Actual_all(ii,5:7)';
        
        
        xyz_int_G=[];
        part_index=0;
        
        %     xyz_int_L=zeros(3,1);
        %     xyz_int_G_C=zeros(3,1);
        %     xyz_int_G_T=zeros(3,1);
        % xyz_int_G=[] ;%zeros(3,1);
        %     xyz_int_G_Roller=zeros(3,1);
        %     xyz_int_L_Roller=zeros(3,1);
        
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
            fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);  % 0 means which reflection
            
            
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
            fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),laser_point_Actual_all(ii,4)*First_hit_energy);
            
            
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
            fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3) ,laser_point_Actual_all(ii,4)*First_hit_energy);
            
            
            
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
            if ~isempty(xyz_int_G2) 
                % plot3( xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3),'y.'  );
                
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                
                % text([xyz_int_G2(1)],[xyz_int_G2(2)],[xyz_int_G2(3)], ['P2']);
                % intersection with cylinder, the reflected ray comes from tape
                
                [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G2,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                
                
                second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                
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
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                    
                elseif n==2
                    %with Mandrel
                    %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                    
                    
                    
                    [~,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                    
                    
                    second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                    
                    fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                end
                
                
            elseif ~isempty(xyz_int_G_T) 
                %       plot3(xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                
                
                
                %                fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3));
                
                [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L2,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                %input are in   local axis system
                second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                
                
                fprintf(fileID2,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_T(1),xyz_int_G_T(2),xyz_int_G_T(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
            elseif ~isempty(xyz_int_G_M) 
                %with Mandrel
                %                     plot3(xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,'%12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                
                %                fprintf(fileID1,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3));
                
                
                %  [Temp_Ref_dir_G,beta2]=refl_cylin_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G);
                
                [Temp_Ref_dir_G,beta2]=refl_Mandrel_plot3D(R_cyl,xyz_int_G_M,Ref_dir_G(1:3,CC),z_cyl_end,fileID6,fileID7,No_graphic);
                
                second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_mandrel(beta2);
                
                %
                fprintf(fileID1,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_M(1),xyz_int_G_M(2),xyz_int_G_M(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
                
            end
            
            
            
            
            
        elseif part_index==1;
            
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
                
                
                %%
                if n==1
                    %with Tape
                    %                  plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                    
                    
                    
                    %                fprintf(fileID2,' %f %f %f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                    
                    
                    [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                    %input are in   local axis system
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                    
                    
                    
                    fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                    
                elseif n==2
                    %with Roller
                    %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                    fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    
                    %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                    
                    %input are in   local axis system
                    
                    
                    
                    [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                    %% should be modified >> how much loss to ambient !!
                    second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_roller(beta2);
                    
                    
                    fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                    
                end
                
            elseif ~isempty(xyz_int_L) 
                %          plot3(xyz_int_G(1),xyz_int_G(2),xyz_int_G(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3));
                
                % intersection with tape, the reflected ray comes from cylinder
                %      fprintf(fileID2,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G(1),xyz_int_G(2),xyz_int_G(3));
                
                [Temp_Ref_dir_L,beta2]=tape_ray_reflection_3D(R_tape,S1,S2,xyz_int_L,x0,y0,Ref_dir_L2,tv,Rot_Roller_axis,fileID8,fileID9,No_graphic);
                %% should be modified >> how much loss to ambient !!
                second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs(beta2);
                
                
                
                fprintf(fileID2,' %12.8f %12.8f %12.8f  %12.8f \r\n',xyz_int_G2(1),xyz_int_G2(2),xyz_int_G2(3), laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
            elseif ~isempty(xyz_int_G_Roller) 
                
                %with Roller
                %                     plot3(xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),'y.');    % it should be done in general m-file, after transfomation into new axis system
                
                % only reflections
                fprintf(fileID11,' %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                
                
                %                fprintf(fileID3,' %12.8f %12.8f %12.8f , 1 \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3));
                
                [Temp_Ref_dir_L,beta2]=refl_Roller_plot3D(R_tape,xyz_int_L_Roller,Ref_dir_L2,tv_R,Rot_Roller_axis,fileID8,fileID9,W_R,No_graphic);
                %% should be modified >> how much loss to ambient !!
                second_hit_energy=Ref_dir_G(4,CC)*(1-First_hit_energy)*Abs_roller(beta2);
                
                
                fprintf(fileID3,' %12.8f %12.8f %12.8f %12.8f \r\n',xyz_int_G_Roller(1),xyz_int_G_Roller(2),xyz_int_G_Roller(3),laser_point_Actual_all(ii,4)*second_hit_energy);
                
                
            end
            
            
        end
        
        
        
        
        
    end
    
    
    
end





end






















% grid  on;
% % axis([-2*R_cyl 2*R_cyl -2*R_cyl 2*R_cyl -z_cyl_end/2 1.5*z_cyl_end]);
% view([180 0]);



fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
%%


fclose(fileID18);
fclose(fileID19);
fclose(fileID28);
fclose(fileID29);


%first intersection are black stars
% legend('yellow for second intersections')
% reflection from cylinder are red
% reflection from the tape are blue



% toc


%% PLOT section of ALL: Rays, intersections
% for reading the Data

fileID4 = fopen('.\int_xyz.txt','r');
fileID5 = fopen('.\Laser_Rays.txt','r');

fileID6 = fopen('.\Normal_vectors_Mandrel.txt','r');
fileID7 = fopen('.\Reflection_vector_Mandrel.txt','r');
fileID8 = fopen('.\Normal_vectors_Tape.txt','r');
fileID9 = fopen('.\Reflection_vector_Tape.txt','r');

fileID11 = fopen('.\Ref_int_xyz.txt','r');



if Graphic_chekbox(7)
    
    XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
    XYZ_int=cell2mat(XYZ_int);
    plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
    
    
    
    XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    XYZ_Rays=cell2mat(XYZ_Rays);
    plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',2 );
    
    
    
    Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
    plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:');
    
    
    
    
    Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
    plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--');
    
    
    
    
    
    Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
    plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
    
    
    
    
    Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
    plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--');
    
    
    
    Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
    Ref_int_xyz=cell2mat(Ref_int_xyz);
    plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');
    
    
    %% Changing the graphical object, Visibility or deleting
    
    Graphics_laser=findobj(h,'color','g','LineStyle',':');
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
    
    
    
    %% trasnparency 
    
    sld_transparency= uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',1,...
       'Units','normalized'...
    ,'Position', [0.09 0.98 0.08 0.02],'Callback',  @(sld_transparency,event) change_transp(sld_transparency)); 



    
    
end

function Visible_plotButtonPushed(rb1,Graphics_laser)
%         x = linspace(0,2*pi,100);
%         y = sin(x);
%         plot(ax,x,y)
if rb1.Value ==false
    set(Graphics_laser,'Visible','off');
else
    set(Graphics_laser,'Visible','on');
end



    function change_transp (sld_transparency)
        data=get(gca,'ch');
        
        for i=4:10
            data(i).Color(4) = get(sld_transparency,'value');
            
        end
        
    

