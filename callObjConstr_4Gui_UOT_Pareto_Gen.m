



function callObjConstr_4Gui_UOT_Pareto_Gen(handles,...
    th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
            R_cyl,...
            W_R,...
            H_indentation,...
                nip_point_M_all,Rot_Roller_axis_all,...
            CV_mesh,Laser_head,...
            Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data)

step_SF(1)=str2num(handles.step_UOT_optimization{1});
step_SF(2)=str2num(handles.step_UOT_optimization{2});
%% Inputs
%% return struc data into numeric and Run


%% Geometrical Parameters

sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
% L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
% th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
% w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
thick_sub=str2double(handles.Geometrical_parameters{10}); % thiness of substrate




%% Process Parameters
materials_Tape=str2num(handles.Process_parameters{1});
materials_sub=str2num(handles.Process_parameters{2});
Velocity=str2double(handles.Process_parameters{3});
Total_energy=str2double(handles.Process_parameters{4});
ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
absorbtion_waste=str2num(handles.Process_parameters{6});
Temp_Right_sub=str2num(handles.Process_parameters{7}); % Includes also mandrel Temperature
% Temp_Right_T=str2num(handles.Process_parameters{8});
 Temp_Right_T_Roller=str2num(handles.Process_parameters{8});
h_conv_Sub=str2double(handles.Process_parameters{9});
h_conv_T=str2num(handles.Process_parameters{10});
Roller_Force=str2double(handles.Process_parameters{11});


% assignin('base','Consolidation_Force',Roller_Force);
% 
% if Roller_Force ~=0
%     
%     if isempty(handles.fit_Func)
%         
%         uiwait(warndlg('No Force-displacement data, Please load the data.'));
%         H_indentation=0;
%         %         Roller_def_Callback(hObject, eventdata, handles);
%     else
%         
%         fitresult= handles.fit_Func;
%         H_indentation=fitresult(Roller_Force);
%         figure(50);
%         plot(Roller_Force,H_indentation,'ko');
%         text (Roller_Force,H_indentation,' Maximum normal deformation');
%         javaFrame    = get(gcf,'JavaFrame');
%         iconFilePath = 'OTOM-icon.png';
%         javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%         
%     end
%     
% else
%     H_indentation=0;
%     
% end
% 
% if H_indentation > R_tape
%     h=errordlg('Normal displacement is more than Roller radius!');
%     H_indentation=R_tape;
% end



%% Computational parameters

step_size_angle=str2double(handles.Computational_parameters{1});  % in radian
No_dev=str2double(handles.Computational_parameters{2});
node_space=str2double(handles.Computational_parameters{3});  % change to the number of devision
Angular_space=str2double(handles.Computational_parameters{4});
L_flat_space=str2double(handles.Computational_parameters{5});  % change to the number of devision
% change to the number of devision
%%

%% Temperature of the mandrel
if length (Temp_Right_sub) ==2
    T_amb_mandrel=Temp_Right_sub(2);
else
    T_amb_mandrel=20;
end
%%  Check graphical outputs
if isempty(handles.checkbox1)
    Graphic_chekbox=ones(1,9);
else
    Graphic_chekbox=[handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,...
        handles.checkbox5,handles.checkbox6,handles.checkbox7,handles.checkbox8,...
        handles.checkbox9];
end


if handles.checkbox_transient_code
    Transient_ID=[handles.checkbox_transient_code, handles.checkbox_Live_output,...
        str2double({handles.edit_init_Temp,...
        handles.edit_inc, handles.edit_Time,handles.checkbox_Video})];
else
      Transient_ID=[0, 0,20,        1, 1,0];
end



Node_z_3D_thermal=handles.Node_z_3D_thermal;
dos('animator_19_frame_camtes.exe -i  &');
%%
Temp_desired_Sub=str2double(get(handles.OP_edit1,'String'));   % desired temperature of the Substrate
Temp_desired_T=str2double(get(handles.OP_edit2,'String'));   % desired temperature of the Tape


Power_desired_Sub=str2double(get(handles.edit_Absorbed_power,'String'));   % desired temperature of the Substrate
Power_desired_Tape=str2double(get(handles.edit_Absorbed_power_T,'String')); 



            
            

% ID_parameter=get(handles.OP_popupmenu1,'Value') ;  % show which parameter will be optimized

Exact_Temp_checkBox=get(handles.OP_checkbox1,'Value') ;
Uniform_Temp_checkBox=get(handles.OP_checkbox2,'Value') ;

%%

% use or not use divergence of the laser
Laser_div_OnOff=get(handles.checkbox_laser_divergence,'Value');

if Laser_div_OnOff
    
    fid23 = fopen('.\Supp_files\Laser_characteristic.txt','r');
    out = textscan(fid23,'%s ','delimiter',',');
    fclose(fid23);
    
    Divergence_factor=str2double(out{1}{2});
    
else
    Divergence_factor=0;
    
end

%%
%%
% Pre analysis

manufacturing_type=handles.manufacturing_type;

oldData=[];


% if Graphic_chekbox(1)
%     fid13 = fopen('.\Supp_files\Postprocess_Parameter.txt');
%     if fid13 ~=-1
%         out = textscan(fid13,'%s','delimiter',',');
%         
%         Measure_Box_Tape=str2num(out{1}{2}); % c, Blx, Bly
%         Measure_Box_Sub=str2num(out{1}{4});
%         
%     else
%         msgdialog ('NO parameter !!!');
%         Measure_Box_Tape=[]; % c, Blx, Bly
%         Measure_Box_Sub=[];
%     end
%     fclose(fid13);
%     
%     
% else
%     Measure_Box_Tape=[]; % c, Blx, Bly
%     Measure_Box_Sub=[];
%     
% end


% if Tecplot_check
fileID_Tecplot_Tape = fopen('.\UOT_Temp_2D_Tape_Tr_Optimized.plt','w');
% end
% if Tecplot_check_T
fileID_Tecplot_Sub3D = fopen('.\UOT_Temp_3D_Sub_Tr_Optimized.plt','w');
fileID_Tecplot_PowerIntensity_UOT_sub = fopen('.\UOT_PowerIntensity_Sub_Optimized.plt','w');
fileID_Tecplot_PowerIntensity_UOT_Tape = fopen('.\UOT_PowerIntensity_Tape_Optimized.plt','w');




%     function  Pre_analysis_UOT_OPT;
%     
%     end
%%
% Precalculation for the Thermal temperature calculation
znode=Node_z_3D_thermal;

if znode==0
    znode=4;
end

% Time=Time_transient; % from user
% inc=Inc_transient;% from user


% node_num_Z=(2*No_dev)-1;

% Rec_moment=round(inc/No_record_location); % when the data should be recorded
%%

mat_size=size(CV_mesh);
ynode=mat_size(4);
xnode=mat_size(3);

node_num_Z=ynode;

% modified 7  Jan 2020
%  L_prim=xnode*norm([CV_mesh(1,1,2,1)-CV_mesh(1,1,1,1) ,CV_mesh(1,2,2,1)-CV_mesh(1,2,1,1),CV_mesh(1,3,2,1)-CV_mesh(1,3,1,1) ] );
%   w=ynode*norm([CV_mesh(1,1,1,2)-CV_mesh(1,1,1,1) ,CV_mesh(1,2,1,2)-CV_mesh(1,2,1,1),CV_mesh(1,3,1,2)-CV_mesh(1,3,1,1) ] );
%  
% modified 31  march 2020
 L_prim=(xnode-1)*norm([CV_mesh(1,1,2,1)-CV_mesh(1,1,1,1) ,CV_mesh(1,2,2,1)-CV_mesh(1,2,1,1),CV_mesh(1,3,2,1)-CV_mesh(1,3,1,1) ] );
  w=(ynode-1)*norm([CV_mesh(1,1,1,2)-CV_mesh(1,1,1,1) ,CV_mesh(1,2,1,2)-CV_mesh(1,2,1,1),CV_mesh(1,3,1,2)-CV_mesh(1,3,1,1) ] );
 



% No_dev_L_Cut=13;
% ynode=No_dev_L_Cut;

% xnode=node_num_Z;
N=ynode*xnode*znode;
N_plane=ynode*xnode;






Temp_Incoming=Temp_Right_sub(1);


% T_left=[ T_nip_Tape; T_nip_Sub];
T_right=0;
T_left=Temp_Incoming;
T_top=0;
T_bottom=0;

BC_Left=[];
BC_top=[];
BC_bottom=[];
BC_right=[];

% for ii=0:znode-1
%     BC_Left=[BC_Left (1:ynode) + (N_plane *ii)];
%     BC_top=[BC_top (ynode:ynode:N_plane)+(N_plane *ii)];
%     BC_bottom=[BC_bottom (1:ynode:N_plane-1)+(N_plane *ii)];
%     BC_right=[BC_right (N_plane:-1:N_plane-ynode+1)+(N_plane *ii) ];
%
% end
%
% BC_bot_Z=1:N_plane;
% BC_top_Z=N_plane*(znode-1)+1:N;



for ii=1:znode-2
    BC_Left=[BC_Left (2:ynode-1) + (N_plane *ii)];
    BC_top=[BC_top (ynode+ynode:ynode:N_plane-ynode)+(N_plane *ii)];
    BC_bottom=[BC_bottom (1+ynode:ynode:N_plane-1-ynode)+(N_plane *ii)];
    BC_right=[BC_right (N_plane-1:-1:N_plane-ynode+2)+(N_plane *ii) ];
    
    
end

BC_bot_Z=[];
BC_top_Z=[];

for jj=1:xnode-2
    BC_bot_Z=[BC_bot_Z (2:ynode-1) +  (ynode *jj)];
    BC_top_Z=[BC_top_Z (2:ynode-1) + (N_plane *(znode-1)) +  (ynode *jj)];
end


% boundary conditions of lines
BC_L_Xleft_Zbot=2:ynode-1;
BC_L_Xleft_Ztop=(2:ynode-1) + (N_plane *(znode-1));
BC_L_Xright_Zbot=(N_plane-1:-1:N_plane-ynode+2);
BC_L_Xright_Ztop=(N_plane-1:-1:N_plane-ynode+2) +(N_plane *(znode-1));

BC_L_Ybot_Zbot=(1+ynode:ynode:N_plane-1-ynode);
BC_L_Ytop_Zbot=(ynode+ynode:ynode:N_plane-ynode);
BC_L_Ybot_Ztop=(1+ynode:ynode:N_plane-1-ynode) +(N_plane *(znode-1));
BC_L_Ytop_Ztop=(ynode+ynode:ynode:N_plane-ynode) +(N_plane *(znode-1));

BC_L_Xleft_Ybot=1 + (N_plane *(1:znode-2));
BC_L_Xleft_Ytop=ynode + (N_plane *(1:znode-2));
BC_L_Xright_Ybot=(N_plane-ynode+1)+(N_plane *(1:znode-2));
BC_L_Xright_Ytop=(N_plane)+(N_plane *(1:znode-2));

% points
BC_P_Xleft_Ybot_Z_bot=1;
BC_P_Xleft_Ytop_Z_bot=ynode;
BC_P_Xleft_Ybot_Z_top= (1) + (N_plane *(znode-1));
BC_P_Xleft_Ytop_Z_top=(ynode) + (N_plane *(znode-1));

BC_P_Xright_Ybot_Z_bot=(N_plane-ynode+1);
BC_P_Xright_Ytop_Z_bot=(N_plane);
BC_P_Xright_Ybot_Z_top= (N_plane-ynode+1) +(N_plane *(znode-1));
BC_P_Xright_Ytop_Z_top=N ; %(N_plane) +(N_plane *(znode-1))






counter=0;
for kk=xnode-1:-1:1
    counter=counter+1;
    index_middle_long(counter)=floor(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=1; %  row number
index_nip_point=(ynode*(row_nip-1)+1):((row_nip)*(ynode));


index_middle_long_Bott=index_middle_long +(N-N_plane);
index_nip_point_Bott=index_nip_point+(N-N_plane);

% roller properties
if manufacturing_type{5}
    
else
manufacturing_type{5}='2100;1000';
    
end





CM_ind_C_i_jp1_k=[];
CM_ind_C_i_jm1_k=[];
CM_ind_C_i_j_km1=[];
CM_ind_C_i_j_kp1=[];
CM_ind_0=[];
CM_ind_C_im1_j_k=[];
CM_ind_C_ip1_j_k=[];

%            first coloumn i index, second j index,

for kk=1:N
    
    
    if mod(kk,ynode)~=0
        if kk <N
            %                 CM(kk,kk+1)=C_i_jp1_k;
            %                 CM(kk+1,kk)=C_i_jm1_k;
            
            CM_ind_C_i_jp1_k=[CM_ind_C_i_jp1_k; kk,kk+1];
            CM_ind_C_i_jm1_k=[CM_ind_C_i_jm1_k;kk+1,kk ];
            
        end
        
    end
    
    if kk <N-ynode+1
        if mod(kk,N_plane)~=0
            
            %                 CM(kk,kk+ynode)=C_ip1_j_k;
            %                 CM(kk+ynode,kk)=C_im1_j_k;
            
            CM_ind_C_ip1_j_k=[CM_ind_C_ip1_j_k; kk,kk+ynode];
            CM_ind_C_im1_j_k=[CM_ind_C_im1_j_k;kk+ynode,kk ];
            
        else
            %                 CM(kk-ynode+1:kk,(kk-ynode+1:kk)+ynode)=0;
            %                 CM((kk-ynode+1:kk)+ynode,kk-ynode+1:kk)=0;
            
            
            
            CM_ind_0=[CM_ind_0;  [ [(kk-ynode+1:kk) , (kk-ynode+1:kk)+ynode] ;
                [(kk-ynode+1:kk)+ynode,  (kk-ynode+1:kk)] ]' ];
            %                  CM_ind_0=[CM_ind_0;(kk-ynode+1:kk)+ynode,kk-ynode+1:kk ];
            
            
        end
    end
    
    if kk <N-(ynode*xnode)+1
        %             CM(kk,kk+(ynode*xnode))=C_i_j_kp1;
        %             CM(kk+(ynode*xnode),kk)=C_i_j_km1;
        
        
        CM_ind_C_i_j_kp1=[CM_ind_C_i_j_kp1;kk,kk+(ynode*xnode)];
        CM_ind_C_i_j_km1=[CM_ind_C_i_j_km1;kk+(ynode*xnode),kk];
        
        
    end
    
end

CM_indexes=[ [1:N;1:N]'  ...
    ;CM_ind_C_i_j_km1; CM_ind_C_i_j_kp1;CM_ind_0;CM_ind_C_im1_j_k;CM_ind_C_ip1_j_k;CM_ind_C_i_jm1_k;CM_ind_C_i_jp1_k];

Len_CM_index=zeros(1,7);
Len_CM_index(1)=length(CM_ind_C_i_j_km1);
Len_CM_index(2)=length(CM_ind_C_i_j_kp1);
Len_CM_index(3)=length(CM_ind_0);
Len_CM_index(4)=length(CM_ind_C_im1_j_k);
Len_CM_index(5)=length(CM_ind_C_ip1_j_k);
Len_CM_index(6)=length(CM_ind_C_i_jm1_k);
Len_CM_index(7)=length(CM_ind_C_i_jp1_k);







% ^Should be based on Roller force! >Nip_Mov
Nip_Mov=0;

Ly=1*w;
Lx=L_prim-Nip_Mov;


% it is assumed we go along the substrate/ certain path, no deviation !
Rel_v_theta= 0; %-sui+th_y;  % - or + of sui is changed



%     znode=Node_z_3D_thermal;
Lz=thick_sub;

delta_x=Lx/(xnode-0); % it was Lx/(xnode-1);
delta_y=Ly/(ynode-0);
delta_z=Lz/(znode-0);


% for 3D representation
if Graphic_chekbox(6)
 advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode);
 
end

fileID1 = fopen('.\nodes.txt');

node = textscan(fileID1,' %*d %f %f %f ','Delimiter',',') ;
Nodes=cell2mat(node);





theta_ind=0;

deg=(deg_tape);

if deg_tape< theta_ind
    deg=(theta_ind);
    
end





th_v=Tape_points_Data{1,4};
%     points_T_G=Tape_points_Data{1,1};
z= Tape_points_Data{1,2};
dl= Tape_points_Data{1,3};



node_num_Z_T=length(z);



L_prim_T=R_tape*(deg-theta_ind)*(pi/180) +L_flat;

%     X_unfolded_T= linspace(0,W_tape,node_num_Z_T);
%     Y_unfolded_T=linspace(0,L_prim_T,length  (th_v)+length(dl));
Ly_T=W_tape;
Lx_T=L_prim_T;

%     delta_x_T=X_unfolded_T(2)-X_unfolded_T(1);
%     delta_y_T=Y_unfolded_T(2)-Y_unfolded_T(1);

% make number of xnode and ynode for Tape Thermal analysis
ynode_T=node_num_Z_T;
xnode_T=length  (th_v)+length(dl);

delta_x_T=Lx_T/xnode_T;
delta_y_T=Ly_T/ynode_T;

N_T=ynode_T*xnode_T;




Temp_Incoming_T=Temp_Right_T_Roller(1);


%   T_amb=25;  % Temperature for convection from surface, not edge!!
T_right_T=0;
T_left_T=Temp_Incoming_T;  % incoming velocity
T_top_T=00;
T_bottom_T=0;

% ^Should be modified >>>>>>>>>>
T_amb_T=20; % out of plane temperature






BC_Left_T=1:ynode_T;
BC_top_T=ynode_T:ynode_T:N_T;
BC_bottom_T=1:ynode_T:N_T-1;
BC_right_T=N_T:-1:N_T-ynode_T+1 ;

indexxy_T=[N_T,N_T-ynode_T+1, 1 ynode_T];

index_middle_long_T=zeros(1,xnode_T-1);
for kk=1:xnode_T-1
    index_middle_long_T(kk)=floor(ynode_T/2)+(ynode_T*kk);
end
% Define a nip-point line
row_nip_T=0; % number before the last row
index_nip_point_T=(ynode_T*(xnode_T-row_nip_T-1)+1):((xnode_T-row_nip_T)*(ynode_T));




% if manufacturing_type{1}
%         Continious_pipe_winding && Layer_number
%     
%     fileID_T_Final_Tape = fopen(manufacturing_type{3});
%     T_Final_Tape = textscan(fileID_T_Final_Tape,' %f ','Delimiter',',','HeaderLines',1);
%     T_Final_Tape=cell2mat(T_Final_Tape);
%     
%     
%     fclose(fileID_T_Final_Tape);
%     
%     if T_Final_Sub_3D
%     if length(T_Final_Tape) ==N_T
%         Temp_T=T_Final_Tape;
%     else
%         Temp_T=init_Temp*ones(N_T,1);
%     end
%     
%     
% else
%     Temp_T=init_Temp*ones(N_T,1);
%     
%     
% end



CM_T_C_i_jp1_T=[];
CM_T_C_i_jm1_T=[];
CM_T_C_ip1_j_T=[];
CM_T_C_im1_j_T=[];


for kk=1:N_T
    
    if mod(kk,ynode_T)~=0
        if kk <N_T
            %                 CM_T(kk,kk+1)=C_i_jp1_T;
            %                 CM_T(kk+1,kk)=C_i_jm1_T;
            
            CM_T_C_i_jp1_T=[CM_T_C_i_jp1_T; kk,kk+1];
            CM_T_C_i_jm1_T=[CM_T_C_i_jm1_T;kk+1,kk ];
        end
        
        %                RHS(kk,1)=-g_dot/delta_x/delta_y;
    end
    
    if kk <N_T-ynode_T+1
        %             CM_T(kk,kk+ynode_T)=C_ip1_j_T;
        %             CM_T(kk+ynode_T,kk)=C_im1_j_T;
        
        CM_T_C_ip1_j_T=[CM_T_C_ip1_j_T; kk,kk+ynode_T];
        CM_T_C_im1_j_T=[CM_T_C_im1_j_T;kk+ynode_T,kk ];
        
        
    end
    
end

CM_indexes_T=[ [1:N_T;1:N_T]'  ...
    ;CM_T_C_im1_j_T; CM_T_C_ip1_j_T;CM_T_C_i_jm1_T;CM_T_C_i_jp1_T];





Len_CM_T_index=zeros(1,4);
Len_CM_T_index(1)=length(CM_T_C_im1_j_T);
Len_CM_T_index(2)=length(CM_T_C_ip1_j_T);
Len_CM_T_index(3)=length(CM_T_C_i_jm1_T);
Len_CM_T_index(4)=length(CM_T_C_i_jp1_T);

%%
    
    if Graphic_chekbox(1)
        fileID_Temp_Red_Box_Sub = fopen('.\Sub3D_Temp_RedBox.txt','w');
        fileID_Temp_Red_Box_Tape = fopen('.\Tape_Temp_RedBox.txt','w');
    else
        fileID_Temp_Red_Box_Sub = [];
        fileID_Temp_Red_Box_Tape = [];
    end
    
    
    
    if Graphic_chekbox(1)
    fid13 = fopen('.\Supp_files\Postprocess_Parameter.txt');
    if fid13 ~=-1
        out = textscan(fid13,'%s','delimiter',',');
        
        Measure_Box_Tape=str2num(out{1}{2}); % c, Blx, Bly
        Measure_Box_Sub=str2num(out{1}{4});
        
    else
        msgdialog ('NO parameter !!!');
        Measure_Box_Tape=[]; % c, Blx, Bly
        Measure_Box_Sub=[];
    end
    fclose(fid13);
    
    
else
    Measure_Box_Tape=[]; % c, Blx, Bly
    Measure_Box_Sub=[];
    
end


    fileID4 = fopen(strcat(handles.UOT_pathfile, sprintf('T_Sub_3D_%d.txt',step_SF(1))),'r');
    fileID5 = fopen(strcat(handles.UOT_pathfile, sprintf('T_Tape_%d.txt',step_SF(1))),'r');
    
    
    T= textscan(fileID4,' %f ','Delimiter',',','HeaderLines',0) ;
    T=cell2mat(T);
    
    Temp_T= textscan(fileID5,' %f ','Delimiter',',','HeaderLines',0) ;
    Temp_T=cell2mat(Temp_T);
    
    fclose(fileID5);
    fclose(fileID4);

steps=step_SF(2)-step_SF(1);
mkdir(handles.UOT_pathfile);


echo off;
% set(handles.text_Opt_Res,'Visible','On')
Objective_num=str2num(handles.Objective_num{:});


fileID_Selected_LaserID = fopen(strcat(handles.UOT_pathfile, sprintf('Selected_Laser_ID.txt')),'w');


% figure(59);
% hold on;
Video = VideoWriter('.\Temp_evoloution_optimized','MPEG-4');
% video.FrameRate = 10;
Video.FrameRate=25; Video.Quality = 100;

open(Video);


   fileID32 = fopen(strcat(handles.UOT_pathfile, sprintf('Optimization_UOT_ObjS.txt')),'w');

for ss=step_SF(1):step_SF(2)
    
    
    % for generating new Temperature data
%          fileID15=fopen(strcat(handles.UOT_pathfile, sprintf('T_Sub_3D_%d.txt',ss)),'w');
% fileID16=fopen(strcat(handles.UOT_pathfile, sprintf('T_Tape_%d.txt',ss)),'w');

    
    % which objective has pririority
    
       fileID31 = fopen(strcat(handles.UOT_pathfile, sprintf('Optimization_UOT_Obj_%d.txt',ss)),'r');
        Objectives= textscan(fileID31,' %f %f %f %f %f %f ','Delimiter',',','HeaderLines',0) ;
    Objectives=cell2mat(Objectives);
    
    % find whcih design is selected 
%     if length(Objective_num)>1
% %         use mean of two columns
%         Mean_Objectives=mean(Objectives (:,Objective_num)')';
%          [m,I]= min(Mean_Objectives (:,:));
%     else
%       [m,I]= min(Objectives (:,Objective_num));
%     end
     if length(Objective_num)>1
        %         use mean of two columns
        % make a scale
        Objectives(:,1)=Objectives(:,1)/mean(Objectives(:,1));

        
         if max(Objectives(:,2))>0                      
        Objectives(:,2)=Objectives(:,2)/mean(Objectives(:,2));
         end
        
        
        if max(Objectives(:,5))>0                      
        Objectives(:,5)=Objectives(:,5)/mean(Objectives(:,5));
        end
          if max(Objectives(:,6))>0                      
        Objectives(:,6)=Objectives(:,6)/mean(Objectives(:,6));
        end

        
        Mean_Objectives=mean(Objectives (:,Objective_num)')';
        [m,I]= min(Mean_Objectives (:,:));
    else
        [m,I]= min(Objectives (:,Objective_num));
    end
    
    
    
    
    
    
       
      fileID30 = fopen(strcat(handles.UOT_pathfile, sprintf('Optimization_UOT_Var_%d.txt',ss)),'r');
   Var= textscan(fileID30,' %f %f %f %f %f %f %f ','Delimiter',',','HeaderLines',0) ;
    Var=cell2mat(Var);
    
    % for paper
%     Std_nip=
%     disp(Objectives(I(1),2));
    
    ID=[3 Var(I(1),:)];
    
    % random the first choice
%         ID=Var(1,:);
    
    
    fclose(fileID30);
      fclose(fileID31);
      
%       Write selected ID file
      
      fprintf(fileID_Selected_LaserID,' %f  %f %f %f %f  %f %f %f  \r\n',ID);
    
    set(handles.text_status,'String',[num2str(((ss-step_SF(1))/steps)*100) '%']);
    drawnow;
    
    % to write optimization summary 
%     fileID22=fopen(strcat(handles.UOT_pathfile, sprintf('Optimization_UOT_%d.txt',ss)),'w');
%     
%  
    
    
    % read start step temperatiure of the tape and substrate
    
    
    
    %%
%     Lb=str2num(handles.Range{1}); % for multiple input variable, it is a combination of all variable
%     Ub=str2num(handles.Range{2});
%     
%     
%     % [x, fval]=ga(@myObjective,1,[] ,[], [], [],Lb,Ub,@ellipseparabola);
%     
%     
%     A = []; b = []; % No linear inequality constraints
%     Aeq = []; beq = []; % No linear equality constraints
% %     'PlotFcns',@gaplotpareto,
% options =optimoptions('ga','FunctionTolerance',1e-3,'ConstraintTolerance',1e-4,'MaxStallGenerations',3,...
% 'MaxStallTime',2);

% 'MaxIterations',1500


%     options = gaoptimset('FunctionTolerance',1e-1,'ConstraintTolerance',1e-1);
%     nvar=length(Lb);
%     
%     [x,Fval,exitFlag,Output]=gamultiobj(@myObjective,nvar,A, ...
%         b,Aeq,beq,Lb,Ub,options);
%     fprintf('The number of points on the Pareto front was: %d\n', size(x,1));
%     fprintf('The number of generations was : %d\n', Output.generations);
% %     figure;
% %     javaFrame    = get(gcf,'JavaFrame');
% %     iconFilePath = 'OTOM-icon.png';
% %     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
%     
% %     plot3(x,Fval(:,1),Fval(:,2),'r*');
%     s=sprintf('Number of design= %d, Pareto design= %d',Output.funccount,size(x,1));
% %     title(['Multi-obj-optimization, ',s]);
%     set(handles.text_Opt_Res,'String',s)
%     xlabel('Var');
%     ylabel('Obj1 (mean difference)');
%     zlabel('Obj2 (STD)');
%     grid on;
    
    %%
    
    % Decide about the best design !!!




%     Tecplot_check=Graphic_chekbox(2);
    %
 
    
    
%     if Tecplot_check_T
%         
%         fclose(fileID_Tecplot_Tape);
%     end
    
    
    
 
    
    
  set(handles.OP_uitable1, 'Data', ID);
%         oldData = get(handles.OP_uitable1,'Data');

        % back to original 
        if isempty(handles.checkbox1)
    Graphic_chekbox=zeros(1,9);
else
    Graphic_chekbox=[handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,...
        handles.checkbox5,handles.checkbox6,handles.checkbox7,handles.checkbox8,...
        handles.checkbox9];
end
        
    
    
%     co=1;
%     ID= x(co,:);
%     while ~ (ID >= Lb & ID <=Ub)
%         co=co+1;
%        ID= x(co,:);
%         
%     end
    
    % Revision 21 Dec 2019
    Total_energy=ID(end);
   
    %%
%     generate data for the next step  based on the optimized value
     delta_t=delta_x/Velocity;
     
     T_old=T;
T_old_T=Temp_T;
    
 [ Temp_tape, Temp_sub,T,Temp_T,~,~] = thermal_domain_generator_UOT_OPT(steps,ss,th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
            sui,L_prim,w,thick_sub,No_dev,R_cyl,...
            W_R,materials_Tape,materials_sub,Total_energy,Velocity,...
            node_space,Angular_space,L_flat_space,...
            Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
            Graphic_chekbox,Transient_ID,handles.manufacturing_type,...
            nip_point_M_all,Rot_Roller_axis_all,...
            handles.UOT_pathfile,CV_mesh,Laser_head,ID(1:end-1),...
            Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data,...
            T,Temp_T,...
               fileID_Tecplot_Tape,fileID_Tecplot_Sub3D,...
        delta_x,delta_y,delta_z,N,N_plane,...
        Nodes,T_right,T_left,T_top,T_bottom,BC_Left,BC_top,BC_bottom,BC_right,...
        BC_bot_Z,BC_top_Z,index_middle_long_Bott,index_nip_point_Bott,CM_indexes,delta_x_T,delta_y_T,...
        T_left_T,T_top_T,T_bottom_T,T_amb_T,BC_Left_T,BC_top_T,BC_bottom_T,BC_right_T,index_middle_long_T,...
        index_nip_point_T,CM_indexes_T,...
        Len_CM_index,Len_CM_T_index,CM_ind_0,...
        index_middle_long,index_nip_point,...
        BC_L_Xleft_Zbot, BC_L_Xleft_Ztop,BC_L_Xright_Zbot,BC_L_Xright_Ztop,...
        BC_L_Ybot_Zbot,BC_L_Ytop_Zbot,BC_L_Ybot_Ztop,BC_L_Ytop_Ztop,...
        BC_L_Xleft_Ybot,BC_L_Xleft_Ytop,BC_L_Xright_Ybot,BC_L_Xright_Ytop,...
        BC_P_Xleft_Ybot_Z_bot, BC_P_Xleft_Ytop_Z_bot,BC_P_Xleft_Ybot_Z_top, BC_P_Xleft_Ytop_Z_top,...
        BC_P_Xright_Ybot_Z_bot, BC_P_Xright_Ytop_Z_bot,BC_P_Xright_Ybot_Z_top,BC_P_Xright_Ytop_Z_top,...
        indexxy_T,...
        Measure_Box_Tape,fileID_Temp_Red_Box_Tape,Lx_T,Ly_T,xnode_T,ynode_T,Lx,Ly,Lz,xnode,ynode,znode,...
        Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,...
        fileID_Tecplot_PowerIntensity_UOT_sub,fileID_Tecplot_PowerIntensity_UOT_Tape,delta_t);
        
    
    
    % Saving new Temperature data
%     fprintf(fileID15,'%f \r\n',T );
% 
% fprintf(fileID16,'%f \r\n',Temp_T );
% %
%    fclose(fileID15);
%  fclose(fileID16);  
 
 
    %%
    figure(59);
cla;
plot(linspace(-Lx,0,length(T(index_middle_long))),T(index_middle_long),'b');
hold on;
plot(linspace(-Lx_T,0,length(Temp_T(index_middle_long_T))),Temp_T(index_middle_long_T),'r');
legend ('Substrate','Tape','Location','northwest');
axis([-Lx_T 0 0 600]);
xlabel('Length [m]');
ylabel('Temperature [^{\circ}C]');
nip=(T(index_middle_long(end))+Temp_T(index_middle_long_T(end)))/2;
str=sprintf(' Temp Nip= %f, time step= %d',nip,ss);
title(str);
% pause(0.1);
Mov(ss-step_SF(1)+1)=getframe(gcf);

% T_old=T;
% T_old_T=Temp_T;

%%
  Temp_sub_nip=Temp_sub(:,end);
        Temp_tape_nip= Temp_tape(:,1);
 f=zeros(6,1);

[T_sub_max,index_max]=max( T_old(index_middle_long));
% should be generalized
[T_tape_max,index_max_T]=max( T_old_T(index_middle_long_T));
%%

S2=(gradient( T_old(index_middle_long))./(delta_x));
S1=(gradient( T(index_middle_long))./(delta_x));

S2_T=(gradient( T_old_T(index_middle_long_T))./(delta_x_T));
S1_T=(gradient( Temp_T(index_middle_long_T))./(delta_x_T));

nabla_rate_sub=  mean(S2(S2>0))-mean(S1(S1>0));
nabla_rate_Tape=  mean(S2_T(S2_T>0))-mean(S1_T(S1_T>0));




%         [T_sub_max,index_max]=max( T_old(index_middle_long));
% should be generalized
%         T_tape_max=max( T_old_T(index_middle_long_T));

% depends on the cooling section and convection with the  toolings
up_slop_Sub=45/(5*delta_x);
low_slop_Sub=-2/(5*delta_x);
sub_slope=(T_sub_max-mean(Temp_sub_nip))/((xnode-index_max)*delta_x);

      

up_slop_T=35/(5*delta_x_T);
low_slop_T=-5/(5*delta_x_T);
Tape_slope=(T_tape_max-mean(Temp_tape_nip))/((xnode_T-index_max_T)*delta_x_T);
        
        
        
if nabla_rate_sub>0 & T_sub_max >1*Temp_desired_Sub
    if sub_slope >up_slop_Sub 
        f(5)=abs(T_sub_max-Temp_desired_Sub)/((xnode-index_max)*delta_x);
%         | 1.0*Temp_desired_Sub-T_sub_max>0
%     c(2)=0.95*Temp_desired_Sub-T_sub_max;
    end
else
     if sub_slope <low_slop_Sub 
%          1.0*Temp_desired_Sub-T_sub_max>0 
        f(5)=-abs((T_sub_max-Temp_desired_Sub)/((xnode-index_max)*delta_x));
%         | 1.0*Temp_desired_Sub-T_sub_max>0
%     c(2)=0.95*Temp_desired_Sub-T_sub_max;
    end
end

if nabla_rate_Tape>0 & T_tape_max >1*Temp_desired_T
    if Tape_slope >up_slop_T
%         T_tape_max-1.2*Temp_desired_T>0 
f(6)=abs(T_tape_max-Temp_desired_T)/((xnode_T-index_max_T)*delta_x_T);
    end
else
    if Tape_slope <low_slop_T
%         1.0*Temp_desired_T-T_tape_max>0
    f(6)=-abs((T_tape_max-Temp_desired_T)/((xnode_T-index_max_T)*delta_x_T));
    end
end


  fprintf(fileID32,' %f  %f %f %f %f  %f   \r\n',f);
      
end


fclose(fileID32);

%%
%    Tecplot_check_T=Graphic_chekbox(2);
   
   
    if Graphic_chekbox(2)
        
        fclose(fileID_Tecplot_Tape);
         fclose(fileID_Tecplot_Sub3D);
      
%         fclose(fileID_Tecplot_Roller_Q_Temp);
        fclose(fileID_Tecplot_PowerIntensity_UOT_sub);
        fclose(fileID_Tecplot_PowerIntensity_UOT_Tape);
   
    end
    
       if Graphic_chekbox (1)
        fclose (fileID_Temp_Red_Box_Sub);
        fclose(fileID_Temp_Red_Box_Tape);
    end

   fclose(fileID_Selected_LaserID);
  
%Kill loading show success!
dos('taskkill /im animator_19_frame_camtes.exe');
writeVideo(Video,Mov);

 
   fclose all;



% fprintf(fileID22,'  %d %f %f  %f %f %f %f    %f %f %f %f \r\n',newRow );


    
end