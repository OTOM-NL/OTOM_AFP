
%%
% Last modification: 26 March 2018
%  fix long substrate length, in-out surface
% Temperature of the mandrel contact with substrate


%%





function thermal_domain_generator_Tr(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,No_dev_L ,w,thick_sub,No_dev,R_cyl,z_cyl_end,...
    tv3,W_R,materials_Tape,materials_sub,Velocity,Total_energy,...
    Rxyz,ID,Laser_head,L_xyz0,absorbtion_waste,...
    node_space,Angular_space,L_flat_space,...
    Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
           Graphic_chekbox,Input_reader_mode,Transient_ID,manufacturing_type,Laser_head_Rot,...
           BRDF_mode,text_status,Divergence_factor,nonlinearMaterials_Tape_Sub)

       
       Live_output=Transient_ID(2);
       init_Temp=Transient_ID(3);
       Inc_transient=Transient_ID(4);
       Time_transient=Transient_ID(5);
       Video_transient=Transient_ID(6);
       
       
       % modified on 9th April 2019
       if manufacturing_type{5}
           
       else
           manufacturing_type{5}='2100;1000';
           
       end

   
       
       
       if Input_reader_mode
           Total_energy=1;
              Graphic_chekbox(1)=1;
%            Graphic_chekbox(3:end)=0;
           
            fileID_input_reader = fopen('.\Inputs_mode.txt');
        [Inputs_Vel_power] = textscan(fileID_input_reader,'%f %f','delimiter',',');
        Vel_Tank=cell2mat(Inputs_Vel_power(1,1)); %/1000;  % should be changed latter for general case >>>>   /1000
        Power_Tank=cell2mat(Inputs_Vel_power(1,2));
        
        
        fclose(fileID_input_reader);
           
       else
    Power_Tank=1;
    Vel_Tank=Velocity;
           
       end
       
       
       
       % Check_MB=get(Output_checks.checkbox1,'Value');
% Check_Tecplot=get(Output_checks.checkbox2,'Value');
% Check_Tape_profile=get(Output_checks.checkbox3,'Value');
% Check_Tape_contour=get(Output_checks.checkbox4,'Value');
% 
% Check_Sub_Profile=get(Output_checks.checkbox5,'Value');
% Check_Sub_contour=get(Output_checks.checkbox6,'Value');
% Check_Config=get(Output_checks.checkbox7,'Value');
% Check_combined_intensity=get(Output_checks.checkbox8,'Value');
% 
% Check_Combined_Temp=get(Output_checks.checkbox9,'Value');
       
     
       
       

if Graphic_chekbox(1)
    fid13 = fopen('.\Supp_files\Postprocess_Parameter.txt');
    if fid13 ~=-1
        out = textscan(fid13,'%s','delimiter',',');
        
        
    else
        msgdialog ('NO parameter !!!');
    end
    fclose(fid13);
    
    Measure_Box_Tape=str2num(out{1}{2}); % c, Blx, Bly
    Measure_Box_Sub=str2num(out{1}{4});
else
      Measure_Box_Tape=[]; % c, Blx, Bly
    Measure_Box_Sub=[];
    
end
%%







fileID21 = fopen('.\Ponits_in_domain.txt','w');
fprintf(fileID21,' X         ,Y         ,Z       \r\n');

%  if Node_z_3D_thermal >1
%          figure(36);
%          figure(46);
%          figure(56);
%      end




% set(h36,'Position', [680 558 100 100 ]);
% set(h46,'Position', [680 558 100 100 ]);
% set(h56,'Position', [680 558 100 100 ]);



assignin('base','Total_energy',Total_energy);
assignin('base','Tape_speed',Velocity);
assignin('base','Tape_Tension',10);

node_num_Z=(2*No_dev)-1;



%  z1 > starting position of the substarte

%  No_dev_L   >>step_size_angle
% node_num_Z >> No_dev  >>
%%
% All the input values should be turned off
%%

% fig_name=get(figure(1),'name');
%
% if get(figure(1),'name') ~='Laser HO finder'

close(figure(1));
% end


close(figure(2));
close(figure(3));
close(figure(10));


%%objectives
% std_Nip_point_Temp_T
% std_Nip_point_Temp_sub
% Nip_point_Temp_sub, Nip_point_Temp_T are the outputs



N_tape = 360; % EVEN number of points, should not be changed !!

% deg_tape=360/8; % degrees for the cylinder section

Tape_Sp=[N_tape;W_tape;R_tape;L_flat;thick_T;deg_tape];  %Tape_Specification



% ID=[2,Gauss_Par_X,Gauss_Par_Y]


%          New optical model
[counter_ray,Rot_Roller_axis,tv]=General_3D_Optical(th_y,Tape_Sp,...
    Rxyz,R_cyl,z_cyl_end,tv3,ID,Laser_head,L_xyz0,absorbtion_waste,W_R,H_indentation,...
      Graphic_chekbox,Laser_head_Rot,BRDF_mode,Divergence_factor);

%%






% >> Pause to let the computer write data into txt file
pause(0.01); %% was 0.05

% for the substrate
fileID1 = fopen('.\Cylinder_ints.txt');
% C = textscan(fileID1,'%d %d %d %d','HeaderLines',1,'Delimiter',',') ;
node = textscan(fileID1,' %f %f %f %f ','Delimiter',',','HeaderLines',1) ;
Nodes_optical=cell2mat(node);

%% for the tape
fileID2 = fopen('.\Tape_ints.txt');
% C = textscan(fileID1,'%d %d %d %d','HeaderLines',1,'Delimiter',',') ;
node_T = textscan(fileID2,' %f %f %f %f ','Delimiter',',','HeaderLines',1) ;
Nodes_optical_T=cell2mat(node_T);


%%

if Graphic_chekbox(1)
fileID_Temp_Red_Box_Sub = fopen('.\Sub3D_Temp_RedBox.txt','w');
fileID_Temp_Red_Box_Tape = fopen('.\Tape_Temp_RedBox.txt','w');
else
   fileID_Temp_Red_Box_Sub = [];
fileID_Temp_Red_Box_Tape = [];
end


Nodes_optical(:,4)=Nodes_optical(:,4)*Total_energy/counter_ray;   %% laser distribution should be implied here
Nodes_optical_T(:,4)=Nodes_optical_T(:,4)*Total_energy/counter_ray;


% points=zeros(No_dev_L*node_num_Z,4);  % store 100 points

h1=figure(1);
set(h1,'Visible','off');
hold on;



theta_ind=acosd((R_tape-H_indentation)/R_tape);  %theta_ind in degree

Nip_Mov =R_tape*sind(theta_ind);


%% Generating thermal points and decide which optical point is in the domain
%we already know that the domain for the tape, do not need to calculate

%% Now the program should decide which part is exposed for laying substrate which depends on Nip-point

% if tv(3) is between 0 and z_cylinder_end
% For cylinder part

%  [points, Boarders]=helical_3D_points(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);

% more than z_cylinder_end or less than zero
% For Dome part
%   [points, Boarders]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);




%%
if R_cyl(1) ~= 0
    
    if tv (3) >=0 && tv (3) <=z_cyl_end
        
        % [points, Boarders,starting_index]=helical_3D_points_4free_onSub(R_cyl(1),z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
        [points, ~ ]=helical_3D_points(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
        
    elseif tv (3) < 0
        % bottom dome part - 1st Dome
        
        [points, ~,~]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
        
    elseif tv (3) > z_cyl_end
        
        % Upper Dome part - 2nd Dome
        [points, ~,~]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
        
        
    end
    
else
    
    [points, ~,~]=sub_on_flat(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
    
    
end


%%




% UpperBoarder= Boarders(:,UpperBoarder_index);
%  BottomBoarder= Boarders(:,BottomBoarder_index);




points=points';

[Len_Nodes_optical,jerk]=size(Nodes_optical);

Points_in_domain=zeros(Len_Nodes_optical,3);
counter=0;

P_temp=zeros(4,3);
E_points=zeros(length(points),1);
E_points_4_adv_func=zeros(length(points),1);
dist_opt_sub=length(points);



%%

No_dev_L_Cut=length(E_points)/node_num_Z;
xnode=node_num_Z;
ynode=No_dev_L_Cut;


X=zeros(xnode,ynode);
Y=zeros(xnode,ynode);
Z=zeros(xnode,ynode);
% Intensity=zeros(xnode,ynode);

for kk=1:xnode
    index=((1:ynode)*xnode)-kk+1;
    X(kk,1:end)=points(index,1)';
    Y(kk,1:end)=points(index,2)';
    Z(kk,1:end)=points(index,3)';
%     Intensity(kk,1:end)=points(index,4)'; % RHS in the thermal model
end



S(:,:,1) = X;
S(:,:,2) = Y;
S(:,:,3) = Z;

T = Nodes_optical(:,1:3); %[Points.XData',Points.YData',Points.ZData'];

d = 1.1*norm(squeeze(S(1,1,:))'-squeeze(S(2,2,:))');

% tic
inSurf =zeros(1,length(T(:,1)));
[s1,s2,~]=size(S);    
Gi = 1;
Gj = 1;
%Main code
for i=1:(s1-Gi)
    for j = 1:(s2-Gj)
        p1 =squeeze(S(i,j,:))';
        p2 =squeeze(S(i,j+Gj,:))';
        p3 =squeeze(S(i+Gi,j,:))';
        p4 =squeeze(S(i+Gi,j+Gj,:))';
        for z=1:size(T,1);
            if inSurf(z) == 0 && norm(T(z,:)-(0.5*(p1+p3)))< d;
                if SideCheck(p1,p2,p3,T(z,:))==1 || SideCheck(p4,p3,p2,T(z,:))==1; 
                    inSurf(z)= z;
                end;
            end
        end;
    end;
end;

inSurf(inSurf==0)=[];











  
for ii=1:length(inSurf)
    
    ind=inSurf(ii); 
     
    
  
 
    
%     if (THETA_all >=359 & THETA_all<=365)
        
        
        counter=counter+1;
        Points_in_domain(counter,1:3)=[Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3)];
        %             plot3(Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3),'cp');
        
        
        fprintf(fileID21,' %12.8f %12.8f %12.8f \r\n',Nodes_optical(ind,1),Nodes_optical(ind,2),Nodes_optical(ind,3));
        
        
        
        % find four neighbouring points
        for dd=1:length(points)
            dist_opt_sub(dd)=norm(abs(Points_in_domain(counter,1:3)- points(dd,1:3)  ));
            
        end
        
        [m,n]=sort(dist_opt_sub);
        % find nearest Z in specified theta angle
        
        P_temp(1,1:3)=points(n(1),1:3)  ;
        P_temp(2,1:3)= points(n(2),1:3) ;
%         P_temp(3,1:3)= points(n(3),1:3) ;
%         P_temp(4,1:3)=points(n(4),1:3) ;
        
        
        
        w1_temp=1/norm([Nodes_optical(ind,1)-P_temp(1,1),...
            Nodes_optical(ind,2)-P_temp(1,2),...
            Nodes_optical(ind,3)-P_temp(1,3)]);
        
        w2_temp=1/norm([Nodes_optical(ind,1)-P_temp(2,1),...
            Nodes_optical(ind,2)-P_temp(2,2),...
            Nodes_optical(ind,3)-P_temp(2,3)]);
        
%         w3_temp=1/norm([Nodes_optical(ind,1)-P_temp(3,1),...
%             Nodes_optical(ind,2)-P_temp(3,2),...
%             Nodes_optical(ind,3)-P_temp(3,3)]);
%         
%         w4_temp=1/norm([Nodes_optical(ind,1)-P_temp(4,1),...
%             Nodes_optical(ind,2)-P_temp(4,2),...
%             Nodes_optical(ind,3)-P_temp(4,3)]);
        
        
        cte=w1_temp+w2_temp; %+w3_temp+w4_temp;
        
        % normalized weighting values for each point based on the distance from the
        % optical point
        w1=w1_temp/cte;
        w2=w2_temp/cte;
%         w3=w3_temp/cte;
%         w4=w4_temp/cte;
        
        % The percentage of energy should be saved in the txt file for each position
        % because of the superposition and some node can get more than 1 time
        % energy, the data should be stored first in a vector and added, then
        % saved in the text file
        
        % then next step will assigning this energy to the corresponding unfolded
        % nodes
        
        
        E=Nodes_optical(ind,4);
        % E vector is filled based on the generation of points which was used
        %  in this file
        E_points(n(1) )=E*w1 +  E_points(n(1));
        E_points(n(2))=E*w2 +        E_points(n(2));
%         E_points( n(3)  )=E*w3 + E_points( n(3));
%         E_points( n(4) )=E*w4 + E_points( n(4));
        
        
%     else
%         
%         %         disp(THETA_all);
%     end
    
    
end

% axis equal;
view([180 20]);
title('Optical-Thermal Model');
hold on;
%transfer energy from optical to the thermal matrices



points(:,4)=E_points;
counter=0;

No_dev_L_Cut=length(E_points)/node_num_Z;

% starting_index=floor(No_dev_L*(Nip_Mov/L_prim)); % after this position , the nip-point will be defined
% No=(No_dev_L-starting_index+1)*No_row_nodes;
% N_aLL=(No_dev_L-starting_index+1)*node_num_Z;
N_aLL=No_dev_L_Cut*node_num_Z;


for jj=1:No_dev_L_Cut
    for ii=node_num_Z:-1:1
        
        INDEX=N_aLL- (jj*node_num_Z)+ii;
        %      INDEX_3D=N_aLL- counter;
        
        
        counter=counter+1;
        E_points_4_adv_func(counter)=E_points(INDEX);
        %         E_points_4_adv_func_3D(counter)=E_points(INDEX_3D);
    end
end

% thus we have the XYZ, and energy for the thermal model now. We jus unfold
% the tape and nodes to get rid of Z coordinate, but we need Z for
% illustartion of final countour on the geometry

% for ii=1:counter
% plot3(Points_in_domain(ii,1),Points_in_domain(ii,2),Points_in_domain(ii,3),'cp')
%
% end
%now using the points inside the boarders to further assignment
%%Function to assign the optical points to the thermal poitns
% change the view angle to the problem to find out this direction of
% looking is the right way for the thermal model without further difficulty


%% unfolded thermal model with heat generation
% X in the unfolded thermal domain is between 0 to w, because of the change
% of angle view,90 degrees rotation, node_num_Z
% Y in the unfolded thermal domain is between 0 to L' ,  length  (th)




Lx=2*w;
Ly=L_prim-Nip_Mov;

% Y_unfolded= linspace(0,w,node_num_Z);
% 
% X_unfolded=linspace(0,L_prim-Nip_Mov,No_dev_L_Cut);

delta_x=Lx/(xnode-1); % it was Lx/(xnode-1);
delta_y=Ly/(ynode-1);
% 
% delta_x=X_unfolded(2)-X_unfolded(1);
% delta_y=Y_unfolded(2)-Y_unfolded(1);
% xnode=node_num_Z;
% ynode=No_dev_L_Cut;
% 
% 
% X=zeros(xnode,ynode);
% Y=zeros(xnode,ynode);
% Z=zeros(xnode,ynode);
Intensity=zeros(xnode,ynode);

for kk=1:xnode
    index=((1:ynode)*xnode)-kk+1;
%     X(kk,1:end)=points(index,1)';
%     Y(kk,1:end)=points(index,2)';
%     Z(kk,1:end)=points(index,3)';
    Intensity(kk,1:end)=points(index,4)'; % RHS in the thermal model
end


   sub_ele_size=delta_x*delta_y;


Intensity=Intensity/sub_ele_size;


% assignin('base','Intensity_Sub',Intensity);


% N=xnode*ynode;

if Graphic_chekbox(8)
h2=figure(2);
set(h2,'Visible','off');

% subplot(1,2,1)
surf(X,Y,Z,Intensity);
% title('Intensity of the optical part-substrate')
hold on;
view([180 0]);


end
% %% CAlculation of the Substrate temperature
% 
% 
% Lx=w;
% Ly=L_prim-Nip_Mov;
% 
% Rel_v_theta=sui+th_y;
% 
% % T_amb_mandrel=55;  % for each layer should be different !!!
% % Temp_Right_sub = T_incoming
% 
% 
% % if Node_z_3D_thermal >1
%     % 3D
%     znode=Node_z_3D_thermal;
%     Lz=thick_sub;
%     
%     
%         Live_output=Transient_ID(2);
%        init_Temp=Transient_ID(3);
%        Inc_transient=Transient_ID(4);
%        Time_transient=Transient_ID(5);
% 
%     
%     E_points_COPY=E_points;
%     
%     for ss=1:length(Power_Tank)
%         
%         E_points=Power_Tank(ss)*E_points_COPY;
%         Velocity=  Vel_Tank(ss);
% 
%             T_matrix_3D=FDM_3D_transient_combined(Ly,Lx,Lz,ynode,xnode,znode,E_points,Velocity,materials_sub,Temp_Right_sub(1),h_conv_Sub,...
%           T_amb_mandrel,Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,Graphic_chekbox(2),Graphic_chekbox(5),Graphic_chekbox(6),...
%               Live_output,init_Temp,Inc_transient,Time_transient);
%         
%  
%         % Linear Material model
% %         T_matrix_3D=thermal_advection_combined_3D(Ly,Lx,Lz,ynode,xnode,znode,E_points,Velocity,materials_sub,Temp_Right_sub(1),h_conv_Sub,...
% %           T_amb_mandrel,Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,Graphic_chekbox(2),Graphic_chekbox(5),Graphic_chekbox(6));
%         
%   
%         % Nonlinear material model
%         % T_matrix_3D=thermal_advection_combined_3D_NONLinear(Ly,Lx,Lz,ynode,xnode,znode,E_points,Velocity,materials_sub,Temp_Right_sub,h_conv_Sub,T_amb_mandrel,Rel_v_theta);
%         
%     end
%     
%     T_matrix=T_matrix_3D;
%     
%     
%     
%      for ii=1:xnode
%         for jj=1:ynode
% %             counter=counter+1;
% %             index=(jj)+((ii-1)*ynode);
%             T_matrix(ii,jj)=T_matrix_3D(ii,ynode-jj+1);
%             
%         end
%      end
%     
%     
% 
%      
% %     
% % else
%     % 2D model
% %     E_points_4_adv_func_COPY=E_points_4_adv_func;
% %     for ss=1:length(Power_Tank)
% %         
% %         E_points_4_adv_func=Power_Tank(ss)*E_points_4_adv_func_COPY;
% %         Velocity=  Vel_Tank(ss);
% %         
% %         T_matrix=thermal_advection_comb_FEM(Ly,Lx,ynode,xnode,E_points_4_adv_func,Velocity,materials_sub,Temp_Right_sub(1),h_conv_Sub,thick_sub,Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,...
% %             Graphic_chekbox(2),Graphic_chekbox(5),Graphic_chekbox(6));
% %         
% %     end
% %     
% %     
% % end
% 
% % assignin('base','Temperature_sub',T_matrix);
% 
% %% Objective for optimization
% % Nip_point_Temp_sub=T_matrix(:,end);
% 
% % std_Nip_point_Temp_sub=std(Nip_point_Temp_sub);
% %%
% 
% 
% if Graphic_chekbox(9)
% 
% 
% % T_matrix(:,end)
% h3=figure(3);
% set(h3,'Visible','off');
% 
% surf(X,Y,Z,T_matrix);
% title('Temperature distribution from thermal part-2D');
% axis equal;
% hold on;
% colorbar;
% view([180 0]);
% end
% 
%  if Graphic_chekbox(7)
% 
% h1=figure(1);
% set(h1,'Visible','off');
% 
% Surf_Sub=surf(X,Y,Z,T_matrix);
%     end

%% Now calculation for Tape
%define the thermal points for the tape
% all the optical points of the tape are in the doimain

% tv=[0;R_cyl(1)-H_indentation;tv3];


R_tape_th=R_tape+thick_T;


h1=figure(1);
set(h1,'Visible','off');

% node_space=0.55;
% Angular_space=3; % in degree
% L_flat_space=node_space;

% th_v indicate which degree is the tangent for the Tape part

[points_T_G,z,dl,th_v]=Tape_points (deg_tape,R_tape_th,W_tape,L_flat,tv,Rot_Roller_axis,W_R,...
    node_space,Angular_space,L_flat_space,theta_ind);

E_points_T=zeros(length(points_T_G),1);
E_points_4_adv_func_T=zeros(length(points_T_G),1);


% limit the search scope
% X_max=max(points_T_G(:,1));
% X_min=min(points_T_G(:,1));
% delta_X=(X_max-X_min)/20;
%
% Y_max=max(points_T_G(:,2));
% Y_min=min(points_T_G(:,2));
% delta_Y=(Y_max-Y_min)/20;
%
% Z_max=max(points_T_G(:,3));
% Z_min=min(points_T_G(:,3));
% delta_Z=(Z_max-Z_min)/20;

points_T_G=points_T_G';
Temp_norm=zeros(length(points_T_G),1);

[row_T,col]=size(Nodes_optical_T);

for ii=1:row_T
    
    % limit the scope search, search in these limits
    
    %     points_T_G(points_T_G(:,1)>)
    %  X_limit= X_op-delta_X  <X_op<X_op+delta_X;
    %  Y_limit=   Y_op-delta_Y  <Y_op<Y_op+delta_Y;
    % Z_limit=     Z_op-delta_Z  <Z_op<Z_op+delta_Z;
    
    %X_limit,Y_limit, and Z_limit should have a same size
    
    % XYZ_limits=[X_limit,Y_limit,Z_limit]
    
    for jj=1:length(points_T_G)
        Temp_norm(jj)=norm(points_T_G(jj,1:3)-Nodes_optical_T(ii,1:3));
    end
    
    [Temp,n]=sort(Temp_norm);
    % store n(1),n(2),n(3),n(4)
    
    P_temp(1,1:3)=points_T_G(n(1),1:3);
    P_temp(2,1:3)= points_T_G(n(2),1:3);
%     P_temp(3,1:3)= points_T_G(n(3),1:3);
%     P_temp(4,1:3)=points_T_G(n(4),1:3);
    
    
    %%
    
    w1_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(1,1),...
        Nodes_optical_T(ii,2)-P_temp(1,2),...
        Nodes_optical_T(ii,3)-P_temp(1,3)]);
    
    w2_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(2,1),...
        Nodes_optical_T(ii,2)-P_temp(2,2),...
        Nodes_optical_T(ii,3)-P_temp(2,3)]);
    
%     w3_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(3,1),...
%         Nodes_optical_T(ii,2)-P_temp(3,2),...
%         Nodes_optical_T(ii,3)-P_temp(3,3)]);
%     
%     w4_temp=1/norm([Nodes_optical_T(ii,1)-P_temp(4,1),...
%         Nodes_optical_T(ii,2)-P_temp(4,2),...
%         Nodes_optical_T(ii,3)-P_temp(4,3)]);
    
    
    cte=w1_temp+w2_temp; %+w3_temp+w4_temp;
    
    % normalized weighting values for each point based on the distance from the
    % optical point
    w1=w1_temp/cte;
    w2=w2_temp/cte;
%     w3=w3_temp/cte;
%     w4=w4_temp/cte;
    
    
    
    E=Nodes_optical_T(ii,4);
    % E vector is filled based on the generation of points which was used
    %  in this file
    E_points_T( n(1))=E*w1 +  E_points_T(n(1)  );
    E_points_T( n(2))=E*w2 +        E_points_T( n(2));
%     E_points_T( n(3))=E*w3 +  E_points_T( n(3));
%     E_points_T( n(4))=E*w4 +  E_points_T( n(4));
    

    
end



node_num_Z_T=length(z);


counter=0;

% for jj=1:length  (th_v)+length(dl)
%     for ii=node_num_Z_T:-1:1
for ii=length(points_T_G):-1:1
    
    counter=counter+1;
    %         index=(jj-1)*node_num_Z_T +ii;
    E_points_4_adv_func_T(counter)=E_points_T(ii);
    
    %     end
end







%% unfolded thermal model with heat generation
% X in the unfolded thermal domain is between 0 to w, because of the change
% of angle view,90 degrees rotation, node_num_Z
% Y in the unfolded thermal domain is between 0 to L' ,  length  (th)



deg=(deg_tape);

if deg_tape< theta_ind
    deg=(theta_ind);
    
end





L_prim_T=R_tape*(deg-theta_ind)*(pi/180) +L_flat;

X_unfolded_T= linspace(0,W_tape,node_num_Z_T);
Y_unfolded_T=linspace(0,L_prim_T,length  (th_v)+length(dl));
Lx_T=W_tape;
Ly_T=L_prim_T;

delta_x_T=X_unfolded_T(2)-X_unfolded_T(1);
delta_y_T=Y_unfolded_T(2)-Y_unfolded_T(1);

% make number of xnode and ynode for Tape Thermal analysis
xnode_T=node_num_Z_T;
ynode_T=length  (th_v)+length(dl);



X_T=zeros(xnode_T,ynode_T);
Y_T=zeros(xnode_T,ynode_T);
Z_T=zeros(xnode_T,ynode_T);
Intensity_T=zeros(xnode_T,ynode_T);

for kk=1:xnode_T
    index=((1:ynode_T)*xnode_T)-kk+1;
    X_T(kk,1:end)=points_T_G(index,1)';
    Y_T(kk,1:end)=points_T_G(index,2)';
    Z_T(kk,1:end)=points_T_G(index,3)';
    Intensity_T(kk,1:end)=E_points_T(index)';
end

Tape_ele_size=delta_x_T*delta_y_T;

Intensity_T=Intensity_T/Tape_ele_size;

% assignin('base','Intensity_Tape',Intensity_T);

% N_T=xnode_T*ynode_T;

if Graphic_chekbox(8)

h2=figure(2);
set(h2,'Visible','off');

% subplot(1,2,2)
surf(X_T,Y_T,Z_T,Intensity_T);
title('Power Intensity (W)');
axis equal;
hold on;
colorbar;
% view([180 0])
view([-180 -60]);


end


%% CAlculation of the temperature



% indicate how many rows of the Tape is tangent to the Roller to have a
% different convection

Row_number_Tape_Roller_tangent=length  (th_v);

% h_conv_T=100;

%%
% Linear Material model


E_points_4_adv_func_T_COPY=E_points_4_adv_func_T;

% for ss=1:length(Power_Tank)
    
    % In this power, initial input power should be unity as 1.00
%     E_points_4_adv_func_T=Power_Tank(ss)*E_points_4_adv_func_T_COPY;
%     Velocity=  Vel_Tank(ss);
%     
% %     if transient_code
%     
%         T_matrix_T=FDM_2D_transient_combined_Tape(Ly_T,Lx_T,ynode_T,xnode_T,E_points_4_adv_func_T,Velocity,materials_Tape,Temp_Right_T,h_conv_T,thick_T,Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box_Tape,...
%               Graphic_chekbox (3),Graphic_chekbox (4),Graphic_chekbox(2));
          
%     else
          
    
%     T_matrix_T=thermal_advection_comb_Tape_FEM(Ly_T,Lx_T,ynode_T,xnode_T,E_points_4_adv_func_T,Velocity,materials_Tape,Temp_Right_T,h_conv_T,thick_T,Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box_Tape,...
%               Graphic_chekbox (3),Graphic_chekbox (4),Graphic_chekbox(2));
          
%     end
          
          
          
    
    
    % Nonlinear Material model
    % T_matrix_T=thermal_advection_comb_Tape_FEM_NONLinear(Ly_T,Lx_T,ynode_T,xnode_T,E_points_4_adv_func_T,Velocity,materials_Tape,Temp_Right_T,h_conv_T,thick_T,Row_number_Tape_Roller_tangent);
    
% end

%%




%% CAlculation of the Substrate temperature


% Lx=2*w;
% Ly=L_prim-Nip_Mov;

Rel_v_theta=-sui+th_y;  % - or + of sui is changed

% T_amb_mandrel=55;  % for each layer should be different !!!
% Temp_Right_sub = T_incoming


% if Node_z_3D_thermal >1
    % 3D
    znode=Node_z_3D_thermal;
    Lz=thick_sub;
    
    
%         Live_output=Transient_ID(2);
%        init_Temp=Transient_ID(3);
%        Inc_transient=Transient_ID(4);
%        Time_transient=Transient_ID(5);

    
%     E_points_COPY=E_points;
    
%     for ss=1:length(Power_Tank)
        
%         E_points=Power_Tank(ss)*E_points_COPY;
        Velocity=  Vel_Tank;
Velocity_T=Velocity;

%    E_points_4_adv_func_T=Power_Tank(ss)*E_points_4_adv_func_T_COPY;

    
%     if transient_code
% linear material model

 


        
       T_matrix_Sub_Tape= FDM_2D_3D_transient_TapeANDSub_combined(Ly,Lx,Lz,ynode,xnode,znode,E_points,Velocity,materials_sub,Temp_Right_sub(1),h_conv_Sub,...
  T_amb_mandrel,Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,Graphic_chekbox(2),Graphic_chekbox(5),Graphic_chekbox(6),...
             Live_output,init_Temp,Inc_transient,Time_transient,...
             Ly_T,Lx_T,ynode_T,xnode_T,E_points_4_adv_func_T,Velocity_T,materials_Tape,Temp_Right_T_Roller,h_conv_T,thick_T,Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box_Tape,...
     Graphic_chekbox (3),Graphic_chekbox (4),Graphic_chekbox(2),...
           Power_Tank, Video_transient,...
             manufacturing_type,W_R,R_tape,...
             text_status,...
            nonlinearMaterials_Tape_Sub );
  
      
        
        
%     end
    
    T_matrix_3D=T_matrix_Sub_Tape{1};
    T_matrix_T=T_matrix_Sub_Tape{2};
    
    
     for ii=1:xnode
        for jj=1:ynode

            T_matrix(ii,jj)=T_matrix_3D(ii,ynode-jj+1);
            
        end
     end
    
    

   
%%


if Graphic_chekbox(9)


% T_matrix(:,end)
h3=figure(3);
set(h3,'Visible','off');

surf(X,Y,Z,T_matrix);
title('Temperature distribution from thermal part-2D');
axis equal;
hold on;
colorbar;
view([180 0]);
end

 if Graphic_chekbox(7)

h1=figure(1);
set(h1,'Visible','off');

Surf_Sub=surf(X,Y,Z,T_matrix);
    end





















if Graphic_chekbox (1)  
fclose (fileID_Temp_Red_Box_Sub);
fclose(fileID_Temp_Red_Box_Tape);
end







% assignin('base','Temperature_Tape',T_matrix_T);



% disp(T_matrix_T);


% Nip_point_Temp_T=T_matrix_T(:,1);   % temperature along nip_point ~ end of the boarder
% std_Nip_point_Temp_T=std(Nip_point_Temp_T);



if Graphic_chekbox(9)
h3=figure(3);
set(h3,'Visible','off');
% subplot(1,2,2)
surf(X_T,Y_T,Z_T,T_matrix_T);
title('Temperature distribution-Tape');
hold on;

view([-180 -60]);


end

if Graphic_chekbox(7)
h1=figure(1);
set(h1,'Visible','off');

Surf_T=surf(X_T,Y_T,Z_T,T_matrix_T);

colorbar;


%%


fileID21 = fopen('.\Ponits_in_domain.txt','r');

XYZ_int = textscan(fileID21,'%f %f %f','Delimiter',',','HeaderLines',1) ;
XYZ_int=cell2mat(XYZ_int);
plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'cp');

% assignin('base','Surf_T',Surf_T);
% assignin('base','Surf_Sub',Surf_Sub);

end



fclose(fileID1);
fclose(fileID2);





%% These two function for the case when sensor parameters are changing
% Define sensor parameter
%            E_points_4_adv_func*(new_Total_energy/Total_energy)
%            Velocity
%          thick_T
%          Temp_Right_sub
%          Temp_Right_T


% T_matrix=thermal_advection_comb_FEM(Ly,Lx,ynode,xnode,E_points_4_adv_func,Velocity,materials_sub,Temp_Right_sub,h_conv_Sub,thick_T);
% T_matrix_T=thermal_advection_comb_Tape_FEM(Ly_T,Lx_T,ynode_T,xnode_T,E_points_4_adv_func_T,Velocity,materials_Tape,Temp_Right_T,h_conv_T,thick_T);

%%

%% for radio bottun
% Normal_lines=findobj(h,'color','y','Marker',':'); % intersection from the reflection





if Graphic_chekbox(7)

% % Surf_Sub
% rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
%     'String','Normal lines','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Normal_lines));
rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,...
    'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Surf_T,Surf_Sub]));

btn5 = uicontrol('Style', 'pushbutton', 'String', 'Clear Tape-Sub',...
    'Position', [520 20 90 40],...
    'Callback', 'delete([Surf_T,Surf_Sub])');


Graphics_point_in_domain=findobj(h1,'color','c','Marker','p');
assignin('base','Graphics_point_in_domain',Graphics_point_in_domain);

rb9 = uicontrol('Style', 'radiobutton','Position',[1040 50 90 40],'Value',true,...
    'String','P_in_domain','Callback', @(rb9,event) Visible_plotButtonPushed(rb9,Graphics_point_in_domain));

btn9 = uicontrol('Style', 'pushbutton', 'String', 'Clear P_in_domain',...
    'Position', [1040 20 90 40],...
    'Callback', 'delete(Graphics_point_in_domain)');



h1=figure(1);
set(h1,'Position', [310 28 860 720 ],'color','w');
else
        close(figure(1));

end

%      cameratoolbar(h1,'Show');
if Graphic_chekbox(8)
h2=figure(2);
view([180 20]);
set(h2,'Position', [620 328 560 420],'color','w');
else
 
close(figure(2));

end

if Graphic_chekbox(4)
% h3=figure(3);
% view([180 20]);
% set(h3,'Position', [630 328 560 420],'color','w');

h100=figure(100);
set(h100,'Position', [650 328 560 420 ],'color','w');

else
        close(figure(100));

% close(figure(3));
end




if Graphic_chekbox(9)
    h3=figure(3);
view([180 20]);
set(h3,'Position', [630 328 560 420],'color','w');
else
    close(figure(3));
end



if  Graphic_chekbox(3) ||  Graphic_chekbox(5)
    h21=figure(21);
    set(h21,'Position', [440 128 1160 420 ],'color','w');
else


close(figure(21));

end



if Graphic_chekbox(6)
   h100=figure(100);
set(h100,'Position', [650 328 560 420 ],'color','w'); 
end


% if 
%       h21=figure(21);
%      set(h21,'Position', [440 128 1160 420 ],'color','w');
% end






if Node_z_3D_thermal >1
    h36= figure(36);
  javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    h46=figure(46);
    javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

    h56=figure(56);
    javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h36,'Position', [660 328 260 220],'color','w');
    set(h46,'Position', [670 328 560 420 ],'color','w');
    set(h56,'Position', [680 328 560 420 ],'color','w');
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




% pcshow   for representing cloud point
% ptCloud = pointCloud(points(:,1:3))
% pcshow(ptCloud)
