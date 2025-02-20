
%%
% Last modification: 20 Dec 2019

%% Transient
% (th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
%     sui,L_prim ,w,thick_sub,No_dev,R_cyl,...
%     W_R,materials_Tape,materials_sub,Total_energy,Velocity,...
%     node_space,Angular_space,L_flat_space,...
%     Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
%     Graphic_chekbox,Transient_ID,manufacturing_type,...
%     nip_point_M_all,Rot_Roller_axis_all,...
%     UOT_pathfile,CV_mesh,Laser_head,ID,...
%     text_status,Points_in_domain_all_sub,Points_in_domain_all_Tape,Tape_points_Data,...
%      T,Temp_T)


function [Temp_tape, Temp_sub,T_old, T_old_T,Power_Abs_Sub,Power_Abs_Tape]=thermal_domain_generator_UOT_OPT(steps,ss,th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,w,thick_sub,No_dev,R_cyl,...
    W_R,materials_Tape,materials_sub,Total_energy,Velocity,...
    node_space,Angular_space,L_flat_space,...
    Temp_Right_sub,Temp_Right_T_Roller,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,....
    Graphic_chekbox,Transient_ID,manufacturing_type,...
    nip_point_M_all,Rot_Roller_axis_all,...
    UOT_pathfile,CV_mesh,Laser_head,ID,...
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
    fileID_Tecplot_PowerIntensity_UOT_sub,fileID_Tecplot_PowerIntensity_UOT_Tape,delta_t)


% 

% tic



Live_output=Transient_ID(2);
init_Temp=Transient_ID(3);
Inc_transient=Transient_ID(4);
Time_transient=Transient_ID(5);
Video_transient=Transient_ID(6);


 Tecplot_check=Graphic_chekbox(2);
%
% Tecplot_check_T=Graphic_chekbox(2);



% UOT_pathfile





Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);
Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);


counter_ray=Laser_head_nx*Laser_head_ny;

%  ID=[3,50,50, 50,50, 0.01,0.01];

% it was modified  on 21 Dec 2019 > problem not normalized 
Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

% steps=length(nip_point_M_all);
%
%
%
% % Precalculation for the Thermal temperature calculation
% znode=Node_z_3D_thermal;
%
% if znode==0
%     znode=4;
% end

% Start_Time= 0; % should be obtained from the last time of previous winding layer

% No_record_location=str2double(manufacturing_type{4});


% % if Tecplot_check
% fileID_Tecplot_Tape = fopen('.\UOT_Temp_2D_Tape_Tr.plt','w');
% % end
% % if Tecplot_check_T
% fileID_Tecplot_Sub3D = fopen('.\UOT_Temp_3D_Sub_Tr.plt','w');
% fileID_Tecplot_PowerIntensity_UOT_sub = fopen('.\UOT_PowerIntensity_Sub.plt','w');
% fileID_Tecplot_PowerIntensity_UOT_Tape = fopen('.\UOT_PowerIntensity_Tape.plt','w');
% end

%
% % Time=Time_transient; % from user
% % inc=Inc_transient;% from user
%

% node_num_Z=(2*No_dev)-1;

% Rec_moment=round(inc/No_record_location); % when the data should be recorded
%%

mat_size=size(CV_mesh);
ynode=mat_size(4);
xnode=mat_size(3);

node_num_Z=ynode;
%
%
%
%
% % No_dev_L_Cut=13;
% % ynode=No_dev_L_Cut;
%
% % xnode=node_num_Z;
% N=ynode*xnode*znode;
% N_plane=ynode*xnode;
%
%
%
%
%
%
% Temp_Incoming=Temp_Right_sub(1);
%
%
% % T_left=[ T_nip_Tape; T_nip_Sub];
% T_right=0;
% T_left=Temp_Incoming;
% T_top=0;
% T_bottom=0;
%
% BC_Left=[];
% BC_top=[];
% BC_bottom=[];
% BC_right=[];
%
% % for ii=0:znode-1
% %     BC_Left=[BC_Left (1:ynode) + (N_plane *ii)];
% %     BC_top=[BC_top (ynode:ynode:N_plane)+(N_plane *ii)];
% %     BC_bottom=[BC_bottom (1:ynode:N_plane-1)+(N_plane *ii)];
% %     BC_right=[BC_right (N_plane:-1:N_plane-ynode+1)+(N_plane *ii) ];
% %
% % end
% %
% % BC_bot_Z=1:N_plane;
% % BC_top_Z=N_plane*(znode-1)+1:N;
%
%
%
% for ii=1:znode-2
%     BC_Left=[BC_Left (2:ynode-1) + (N_plane *ii)];
%     BC_top=[BC_top (ynode+ynode:ynode:N_plane-ynode)+(N_plane *ii)];
%     BC_bottom=[BC_bottom (1+ynode:ynode:N_plane-1-ynode)+(N_plane *ii)];
%     BC_right=[BC_right (N_plane-1:-1:N_plane-ynode+2)+(N_plane *ii) ];
%
%
% end
%
% BC_bot_Z=[];
% BC_top_Z=[];
%
% for jj=1:xnode-2
%     BC_bot_Z=[BC_bot_Z (2:ynode-1) +  (ynode *jj)];
%     BC_top_Z=[BC_top_Z (2:ynode-1) + (N_plane *(znode-1)) +  (ynode *jj)];
% end
%
%
% % boundary conditions of lines
% BC_L_Xleft_Zbot=2:ynode-1;
% BC_L_Xleft_Ztop=(2:ynode-1) + (N_plane *(znode-1));
% BC_L_Xright_Zbot=(N_plane-1:-1:N_plane-ynode+2);
% BC_L_Xright_Ztop=(N_plane-1:-1:N_plane-ynode+2) +(N_plane *(znode-1));
%
% BC_L_Ybot_Zbot=(1+ynode:ynode:N_plane-1-ynode);
% BC_L_Ytop_Zbot=(ynode+ynode:ynode:N_plane-ynode);
% BC_L_Ybot_Ztop=(1+ynode:ynode:N_plane-1-ynode) +(N_plane *(znode-1));
% BC_L_Ytop_Ztop=(ynode+ynode:ynode:N_plane-ynode) +(N_plane *(znode-1));
%
% BC_L_Xleft_Ybot=1 + (N_plane *(1:znode-2));
% BC_L_Xleft_Ytop=ynode + (N_plane *(1:znode-2));
% BC_L_Xright_Ybot=(N_plane-ynode+1)+(N_plane *(1:znode-2));
% BC_L_Xright_Ytop=(N_plane)+(N_plane *(1:znode-2));
%
% % points
% BC_P_Xleft_Ybot_Z_bot=1;
% BC_P_Xleft_Ytop_Z_bot=ynode;
% BC_P_Xleft_Ybot_Z_top= (1) + (N_plane *(znode-1));
% BC_P_Xleft_Ytop_Z_top=(ynode) + (N_plane *(znode-1));
%
% BC_P_Xright_Ybot_Z_bot=(N_plane-ynode+1);
% BC_P_Xright_Ytop_Z_bot=(N_plane);
% BC_P_Xright_Ybot_Z_top= (N_plane-ynode+1) +(N_plane *(znode-1));
% BC_P_Xright_Ytop_Z_top=N ; %(N_plane) +(N_plane *(znode-1))
%




%%



%
% counter=0;
% for kk=xnode-1:-1:1
%     counter=counter+1;
%     index_middle_long(counter)=floor(ynode/2)+(ynode*kk);
% end
% % Define a nip-point line
% row_nip=1; %  row number
% index_nip_point=(ynode*(row_nip-1)+1):((row_nip)*(ynode));
%
%
% index_middle_long_Bott=index_middle_long +(N-N_plane);
% index_nip_point_Bott=index_nip_point+(N-N_plane);

% T=init_Temp*ones(N,1);

% if manufacturing_type{1}
%
%
%
%     % READ Temperature from last layer
%     fileID_Temp_final_sub = fopen( manufacturing_type{2});
%     T_Final_Sub_3D = textscan(fileID_Temp_final_sub,' %f ','Delimiter',',','HeaderLines',1);
%     T_Final_Sub_3D=cell2mat(T_Final_Sub_3D);
%
%     fclose(fileID_Temp_final_sub);
%
%
%
%     % if T_Final_Sub_3D
%     if length (T_Final_Sub_3D) ==N
%         T=T_Final_Sub_3D;
%     else
%         T=init_Temp*ones(N,1);
%     end
%     % end
%
% else
%     T=init_Temp*ones(N,1);
%
%
% end


% % roller properties
% if manufacturing_type{5}
%
% else
% manufacturing_type{5}='2100;1000';
%
% end
%
%
%
%
%
% CM_ind_C_i_jp1_k=[];
% CM_ind_C_i_jm1_k=[];
% CM_ind_C_i_j_km1=[];
% CM_ind_C_i_j_kp1=[];
% CM_ind_0=[];
% CM_ind_C_im1_j_k=[];
% CM_ind_C_ip1_j_k=[];
%
% %            first coloumn i index, second j index,
%
% for kk=1:N
%
%
%     if mod(kk,ynode)~=0
%         if kk <N
%             %                 CM(kk,kk+1)=C_i_jp1_k;
%             %                 CM(kk+1,kk)=C_i_jm1_k;
%
%             CM_ind_C_i_jp1_k=[CM_ind_C_i_jp1_k; kk,kk+1];
%             CM_ind_C_i_jm1_k=[CM_ind_C_i_jm1_k;kk+1,kk ];
%
%         end
%
%     end
%
%     if kk <N-ynode+1
%         if mod(kk,N_plane)~=0
%
%             %                 CM(kk,kk+ynode)=C_ip1_j_k;
%             %                 CM(kk+ynode,kk)=C_im1_j_k;
%
%             CM_ind_C_ip1_j_k=[CM_ind_C_ip1_j_k; kk,kk+ynode];
%             CM_ind_C_im1_j_k=[CM_ind_C_im1_j_k;kk+ynode,kk ];
%
%         else
%             %                 CM(kk-ynode+1:kk,(kk-ynode+1:kk)+ynode)=0;
%             %                 CM((kk-ynode+1:kk)+ynode,kk-ynode+1:kk)=0;
%
%
%
%             CM_ind_0=[CM_ind_0;  [ [(kk-ynode+1:kk) , (kk-ynode+1:kk)+ynode] ;
%                 [(kk-ynode+1:kk)+ynode,  (kk-ynode+1:kk)] ]' ];
%             %                  CM_ind_0=[CM_ind_0;(kk-ynode+1:kk)+ynode,kk-ynode+1:kk ];
%
%
%         end
%     end
%
%     if kk <N-(ynode*xnode)+1
%         %             CM(kk,kk+(ynode*xnode))=C_i_j_kp1;
%         %             CM(kk+(ynode*xnode),kk)=C_i_j_km1;
%
%
%         CM_ind_C_i_j_kp1=[CM_ind_C_i_j_kp1;kk,kk+(ynode*xnode)];
%         CM_ind_C_i_j_km1=[CM_ind_C_i_j_km1;kk+(ynode*xnode),kk];
%
%
%     end
%
% end
%
% CM_indexes=[ [1:N;1:N]'  ...
%     ;CM_ind_C_i_j_km1; CM_ind_C_i_j_kp1;CM_ind_0;CM_ind_C_im1_j_k;CM_ind_C_ip1_j_k;CM_ind_C_i_jm1_k;CM_ind_C_i_jp1_k];
%
% Len_CM_index=zeros(1,7);
% Len_CM_index(1)=length(CM_ind_C_i_j_km1);
% Len_CM_index(2)=length(CM_ind_C_i_j_kp1);
% Len_CM_index(3)=length(CM_ind_0);
% Len_CM_index(4)=length(CM_ind_C_im1_j_k);
% Len_CM_index(5)=length(CM_ind_C_ip1_j_k);
% Len_CM_index(6)=length(CM_ind_C_i_jm1_k);
% Len_CM_index(7)=length(CM_ind_C_i_jp1_k);
%




% %^Should be based on Roller force! >Nip_Mov
% Nip_Mov=0;
%
% Ly=2*w;
% Lx=L_prim-Nip_Mov;
%
%
% % it is assumed we go along the substrate/ certain path, no deviation !
% Rel_v_theta= 0; %-sui+th_y;  % - or + of sui is changed
%
%
%
% %     znode=Node_z_3D_thermal;
% Lz=thick_sub;
%
% delta_x=Lx/(xnode-0); % it was Lx/(xnode-1);
% delta_y=Ly/(ynode-0);
% delta_z=Lz/(znode-0);
%
%
% % for 3D representation
% if Graphic_chekbox(6)
%  advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode);
%
% end
%
% fileID1 = fopen('.\nodes.txt');
%
% node = textscan(fileID1,' %*d %f %f %f ','Delimiter',',') ;
% Nodes=cell2mat(node);



% convection_term=-h*A/(xnode*ynode);%*((delta_x*delta_y*delta_z)^2);
% convection_term=convection_term/(delta_x*delta_y*delta_z);% because delta_z is constant


%  figure(61);  % transient data output representation
% title('Substrate and Tape');
% figure(61);  % transient data output representation
%     javaFrame    = get(gcf,'JavaFrame');
%     iconFilePath = 'OTOM-icon.png';
%     javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));



%
% theta_ind=0;
%
% deg=(deg_tape);
%
% if deg_tape< theta_ind
%     deg=(theta_ind);
%
% end
%
%
%
%
%
% th_v=Tape_points_Data{1,4};
% %     points_T_G=Tape_points_Data{1,1};
% z= Tape_points_Data{1,2};
% dl= Tape_points_Data{1,3};


%
%
%
%
%
% node_num_Z_T=length(z);
%
%
%
% L_prim_T=R_tape*(deg-theta_ind)*(pi/180) +L_flat;
%
% %     X_unfolded_T= linspace(0,W_tape,node_num_Z_T);
% %     Y_unfolded_T=linspace(0,L_prim_T,length  (th_v)+length(dl));
% Ly_T=W_tape;
% Lx_T=L_prim_T;
%
% %     delta_x_T=X_unfolded_T(2)-X_unfolded_T(1);
% %     delta_y_T=Y_unfolded_T(2)-Y_unfolded_T(1);
%
% % make number of xnode and ynode for Tape Thermal analysis
% ynode_T=node_num_Z_T;
% xnode_T=length  (th_v)+length(dl);
%
% delta_x_T=Lx_T/xnode_T;
% delta_y_T=Ly_T/ynode_T;
%
% N_T=ynode_T*xnode_T;
%
%
%
%
% Temp_Incoming_T=Temp_Right_T_Roller(1);
%
%
% %   T_amb=25;  % Temperature for convection from surface, not edge!!
% T_right_T=0;
% T_left_T=Temp_Incoming_T;  % incoming velocity
% T_top_T=00;
% T_bottom_T=0;
%
% % ^Should be modified >>>>>>>>>>
% T_amb_T=20; % out of plane temperature
%
%
%
%
%
%
% BC_Left_T=1:ynode_T;
% BC_top_T=ynode_T:ynode_T:N_T;
% BC_bottom_T=1:ynode_T:N_T-1;
% BC_right_T=N_T:-1:N_T-ynode_T+1 ;
%
% indexxy_T=[N_T,N_T-ynode_T+1, 1 ynode_T];
%
% index_middle_long_T=zeros(1,xnode_T-1);
% for kk=1:xnode_T-1
%     index_middle_long_T(kk)=floor(ynode_T/2)+(ynode_T*kk);
% end
% % Define a nip-point line
% row_nip_T=0; % number before the last row
% index_nip_point_T=(ynode_T*(xnode_T-row_nip_T-1)+1):((xnode_T-row_nip_T)*(ynode_T));
%
%

% load T_Final_Tape;

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




%
%
% CM_T_C_i_jp1_T=[];
% CM_T_C_i_jm1_T=[];
% CM_T_C_ip1_j_T=[];
% CM_T_C_im1_j_T=[];
%
%
% for kk=1:N_T
%
%     if mod(kk,ynode_T)~=0
%         if kk <N_T
%             %                 CM_T(kk,kk+1)=C_i_jp1_T;
%             %                 CM_T(kk+1,kk)=C_i_jm1_T;
%
%             CM_T_C_i_jp1_T=[CM_T_C_i_jp1_T; kk,kk+1];
%             CM_T_C_i_jm1_T=[CM_T_C_i_jm1_T;kk+1,kk ];
%         end
%
%         %                RHS(kk,1)=-g_dot/delta_x/delta_y;
%     end
%
%     if kk <N_T-ynode_T+1
%         %             CM_T(kk,kk+ynode_T)=C_ip1_j_T;
%         %             CM_T(kk+ynode_T,kk)=C_im1_j_T;
%
%         CM_T_C_ip1_j_T=[CM_T_C_ip1_j_T; kk,kk+ynode_T];
%         CM_T_C_im1_j_T=[CM_T_C_im1_j_T;kk+ynode_T,kk ];
%
%
%     end
%
% end
%
% CM_indexes_T=[ [1:N_T;1:N_T]'  ...
%     ;CM_T_C_im1_j_T; CM_T_C_ip1_j_T;CM_T_C_i_jm1_T;CM_T_C_i_jp1_T];
%
%
%
%
%
% Len_CM_T_index=zeros(1,4);
% Len_CM_T_index(1)=length(CM_T_C_im1_j_T);
% Len_CM_T_index(2)=length(CM_T_C_ip1_j_T);
% Len_CM_T_index(3)=length(CM_T_C_i_jm1_T);
% Len_CM_T_index(4)=length(CM_T_C_i_jp1_T);






% set initial Temp of tape, substrate
T_old=T;
T_old_T=Temp_T;





%
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

% close(figure(1));
% close(figure(2));
% close(figure(3));
% close(figure(10));


%%objectives
% std_Nip_point_Temp_T
% std_Nip_point_Temp_sub
% Nip_point_Temp_sub, Nip_point_Temp_T are the outputs



%     fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'r');
%     fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'r');





%%   Estimate the temperature of the roller during the process
%

% Temp_roller=Temp_Right_T_Roller(2);
if manufacturing_type{5}
    Roller_rho_cp=str2num(manufacturing_type{5});
    ro_roller=Roller_rho_cp(1); %2180;
    cp_roller=Roller_rho_cp(2);
    
else
    ro_roller=2100; %2180;
    cp_roller=1e3;
end

% delta_t=delta_x/v_original;


Q_roller=zeros(ss+steps,1);
Temp_roller_current=Temp_Right_T_Roller(2); %+ (4*40);


% for ss=1:steps
% fileID_Q_roller= fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'r');
% %      fileID3 = fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'r');
% Q_roller_temp = textscan(fileID_Q_roller,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1);
% Q_roller_temp=cell2mat(Q_roller_temp);
% fclose(fileID_Q_roller);
%
% Ray_index=Q_roller_temp(:,2);
%
% % incoming heat flux at every step
% % sum( Q_roller_temp(:,1)'*  Power_Actual(Ray_index) );
% % Automatically sum >>
% Q_roller(ss)=Q_roller_temp(:,1)'*  Power_Actual(Ray_index) ;
% end


%    if Tecplot_check
%        fileID_Tecplot_Roller_Q_Temp = fopen('.\UOT_Roller_Q_Temp.txt','w');
%        fprintf(fileID_Tecplot_Roller_Q_Temp,'Q(heat flux), Temp \r\n');
%
%     end










% This loop is a ttransient loop, for every step, movement, every step in
% time !!

% for ss=1:steps

%     set(text_status,'String',[num2str((ss/steps)*100) '%']);
%     drawnow;


% Temp Roller emtimation increase
% time=delta_x/Velocity;
Temp_roller_increase=0; %Roller_temp_estimation(ro_roller,R_tape,W_R,cp_roller,time,Q_roller(ss));
Temp_roller_current=Temp_roller_current+Temp_roller_increase ;  %+ linspace(0,Temp_roller_increase,inc)+;







th_v=Tape_points_Data{ss,4};
points_T_G=Tape_points_Data{ss,1};
z= Tape_points_Data{ss,2};
dl= Tape_points_Data{ss,3};






% >> Pause to let the computer write data into txt file
%     pause(0.01); %% was 0.05

% for the substrate
% fileID1 = fopen('.\Cylinder_ints.txt');
% C = textscan(fileID1,'%d %d %d %d','HeaderLines',1,'Delimiter',',') ;

% S_ind >> source index
%     node = textscan(fileID1,' %f %f %f %f %f','Delimiter',',','HeaderLines',1) ;
%     Nodes_optical=cell2mat(node);

%      S_ind_Nodes_optical=cell2mat(node(5));

%% for the tape
% fileID2 = fopen('.\Tape_ints.txt');
% C = textscan(fileID1,'%d %d %d %d','HeaderLines',1,'Delimiter',',') ;
%     node_T = textscan(fileID2,' %f %f %f %f %f','Delimiter',',','HeaderLines',1) ;
%     Nodes_optical_T=cell2mat(node_T);






Nodes_optical=Points_in_domain_all_sub{ss}  ;
Nodes_optical_T= Points_in_domain_all_Tape{ss} ;


%     if Graphic_chekbox(1)
%         fileID_Temp_Red_Box_Sub = fopen('.\Sub3D_Temp_RedBox.txt','w');
%         fileID_Temp_Red_Box_Tape = fopen('.\Tape_Temp_RedBox.txt','w');
%     else
%         fileID_Temp_Red_Box_Sub = [];
%         fileID_Temp_Red_Box_Tape = [];
%     end

%%
Nodes_optical(:,4)=Nodes_optical(:,4)*Total_energy/counter_ray;   %% laser distribution should be implied here
Nodes_optical_T(:,4)=Nodes_optical_T(:,4)*Total_energy/counter_ray;
%%


% points=zeros(No_dev_L*node_num_Z,4);  % store 100 points

%     h1=figure(1);
%     set(h1,'Visible','off');
%     hold on;



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
% Boarders contains only X and Z coordinates ?!


% if R_cyl(1) ~= 0
%
%     if tv (3) >=0 && tv (3) <=z_cyl_end
%
%         % [points, Boarders,starting_index]=helical_3D_points_4free_onSub(R_cyl(1),z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
%         [points, Boarders]=helical_3D_points(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
%
%     elseif tv (3) < 0
%         % bottom dome part - 1st Dome
%
%         [points, Boarders,~]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
%
%     elseif tv (3) > z_cyl_end
%
%         % Upper Dome part - 2nd Dome
%         [points, Boarders,~]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
%
%
%     end
%
% else
%
%     [points, Boarders,~]=sub_on_flat(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv);
%
%
% end


mat_size=size(CV_mesh);

X=reshape(CV_mesh(ss,1,:,:),mat_size(3:4));
Y= reshape(CV_mesh(ss,2,:,:),mat_size(3:4));
Z=reshape(CV_mesh(ss,3,:,:),mat_size(3:4));

% all Eulerian points
points=[X(:),Y(:),Z(:)];

% make boarder points > CCW
%         Ps1=[X(1:end,1),Y(1:end,1),Z(1:end,1) ];
%         Ps2=[X(end,2:end)' ,Y(end,2:end)' ,Z(end,2:end)' ];
%         Ps3=[X(end-1:-1:1,end),Y(end-1:-1:1,end),Z(end-1:-1:1,end)];
%         Ps4=[X(1,end-1:-1:2)',Y(1,end-1:-1:2)',Z(1,end-1:-1:2)' ];

%    Boarders=[Ps1;Ps2;Ps3;Ps4]';


%    plot3(Boarders(1,:),Boarders(2,:),Boarders(3,:),'c--')
%    axis equal;
%    hold on;
%%




%     if R_cyl(1) ~= 0
%
%         UpperBoarder_index=Boarders(2,:)>0;
%         BottomBoarder_index=Boarders(2,:)<0;
%
%         UpperBoarder= Boarders(:,UpperBoarder_index);
%         BottomBoarder= Boarders(:,BottomBoarder_index);
%
% %     elseif R_cyl(1) == 0
% %         %         UpperBoarder_index=Boarders;
% %         %  BottomBoarder_index=Boarders;
% %
% %         %
% %         UpperBoarder= Boarders;
% %         BottomBoarder= Boarders;
%
%     end



% UpperBoarder= Boarders(:,UpperBoarder_index);
%  BottomBoarder= Boarders(:,BottomBoarder_index);




%     points=points';

[Len_Nodes_optical,jerk]=size(Nodes_optical);
%
%     Points_in_domain=zeros(Len_Nodes_optical,3);
%     counter=0;
%
%     P_temp=zeros(4,3);
E_points=zeros(length(points),1);
E_points_4_adv_func=zeros(length(points),1);
%     dist_opt_sub=length(points);


for ii=1:Len_Nodes_optical
    
    % read neighboring indices and weight functions
%         out= Nodes_optical(ii,6:13);
        out= Nodes_optical(ii,6:9);
%         n=out(1:4);
%         Weight=out(5:8);
         n=out(1:2);
        Weight=out(3:4);
    
    
    
    % Apply different laser distribution
    % 26 March-2019
    % find  index of ray source
    Ray_index=Nodes_optical(ii,5);
    
    E=Nodes_optical(ii,4) * Power_Actual(Ray_index);
    
    
    
    % E vector is filled based on the generation of points which was used
    %  in this file
    E_points(n(1))=E*Weight(1) +  E_points(n(1));
    E_points(n(2))=E*Weight(2) +   E_points(n(2));
%     E_points(n(3))=E*Weight(3) +  E_points( n(3));
%     E_points(n(4))=E*Weight(4) +  E_points( n(4));
    
    
    
    
end

% axis equal;
% view([180 20]);
% title('Optical-Thermal Model');
% hold on;
%transfer energy from optical to the thermal matrices



%     points(:,4)=E_points;
counter=0;


%%   >>>>>>>>>>>>>>>>>>>>>>>
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


%% unfolded thermal model with heat generation
% X in the unfolded thermal domain is between 0 to w, because of the change
% of angle view,90 degrees rotation, node_num_Z
% Y in the unfolded thermal domain is between 0 to L' ,  length  (th)




%% Show part

%     Y_unfolded= linspace(0,w,node_num_Z);
%
%     X_unfolded=linspace(0,L_prim-Nip_Mov,No_dev_L_Cut);



%     delta_x=X_unfolded(2)-X_unfolded(1);
%     delta_y=Y_unfolded(2)-Y_unfolded(1);




%     X=zeros(xnode,ynode);
%     Y=zeros(xnode,ynode);
%     Z=zeros(xnode,ynode);
%     Intensity=zeros(xnode,ynode);

%     for kk=1:xnode
%         index=((1:ynode)*xnode)-kk+1;
%         X(kk,1:end)=points(index,1)';
%         Y(kk,1:end)=points(index,2)';
%         Z(kk,1:end)=points(index,3)';
%         Intensity(kk,1:end)=E_points (index)'; % RHS in the thermal model
%     end

% assignin('base','Intensity_Sub',Intensity);


% N=xnode*ynode;

% if Graphic_chekbox(8)
% h2=figure(2);
% set(h2,'Visible','off');
%
% % subplot(1,2,1)
% surf(X,Y,Z,Intensity);
% % title('Intensity of the optical part-substrate')
% hold on;
% view([180 0]);
%
%
% end
% %% CAlculation of the Substrate temperature



%% Now calculation for Tape
%define the thermal points for the tape
% all the optical points of the tape are in the doimain

% th_v indicate which degree is the tangent for the Tape part


%% Get from Optical_UOT.mat

%       Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
%       Roller_Pos_TV=nip_point_M_all(ss,:,:);


%     [points_T_G,z,dl,th_v]=Tape_points (deg_tape,R_tape_th,W_tape,L_flat,Roller_Pos_TV,Rot_Roller_axis,W_R,...
%         node_space,Angular_space,L_flat_space,theta_ind);

E_points_T=zeros(length(points_T_G),1);
E_points_4_adv_func_T=zeros(length(points_T_G),1);




points_T_G=points_T_G';
%     Temp_norm=zeros(length(points_T_G),1);

[row_T,col]=size(Nodes_optical_T);

for ii=1:row_T
    
    
 %         out= Nodes_optical_T(ii,6:13);
%         
%         %                 out= Nodes_optical(ii,6:13);
%         n=out(1:4);
%         Weight=out(5:8);
        out= Nodes_optical_T(ii,6:9);
        
        %                 out= Nodes_optical(ii,6:13);
        n=out(1:2);
        Weight=out(3:4);
    
    
    
    Ray_index=Nodes_optical_T(ii,5);
    
    E=Nodes_optical_T(ii,4) * Power_Actual(Ray_index);
    
    
    % E vector is filled based on the generation of points which was used
    %  in this file
    E_points_T( n(1))=E*Weight(1) +  E_points_T(n(1)  );
    E_points_T( n(2))=E*Weight(2) +  E_points_T( n(2));
%     E_points_T( n(3))=E*Weight(3) +  E_points_T( n(3));
%     E_points_T( n(4))=E*Weight(4) +  E_points_T( n(4));
    
    
    
end






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







%% Show part






%     X_T=zeros(xnode_T,ynode_T);
%     Y_T=zeros(xnode_T,ynode_T);
%     Z_T=zeros(xnode_T,ynode_T);
%     Intensity_T=zeros(xnode_T,ynode_T);
%
%     for kk=1:xnode_T
%         index=((1:ynode_T)*xnode_T)-kk+1;
%         X_T(kk,1:end)=points_T_G(index,1)';
%         Y_T(kk,1:end)=points_T_G(index,2)';
%         Z_T(kk,1:end)=points_T_G(index,3)';
%         Intensity_T(kk,1:end)=E_points_T(index)';
%     end


% assignin('base','Intensity_Tape',Intensity_T);

 N_T=xnode_T*ynode_T;

% if Graphic_chekbox(8)
%
% h2=figure(2);
% set(h2,'Visible','off');
%
% % subplot(1,2,2)
% surf(X_T,Y_T,Z_T,Intensity_T);
% title('Power Intensity (W)');
% axis equal;
% hold on;
% colorbar;
% % view([180 0])
% view([-180 -60]);
%
%
% end


%% CAlculation of the temperature



% indicate how many rows of the Tape is tangent to the Roller to have a
% different convection

Row_number_Tape_Roller_tangent=length  (th_v);

% h_conv_T=100;

%%
% Linear Material model


% E_points_4_adv_func_T_COPY=E_points_4_adv_func_T;





%     Velocity=  Vel_Tank;
Velocity_T=Velocity;

% end


%     Temp_roller=Temp_Right_T_Roller(2);

    if Tecplot_check
         sub_ele_size=delta_x*delta_y;
Tape_ele_size=delta_x_T*delta_y_T;

    Tecplot_Matlab_creator_2D (ss,N_plane,xnode,ynode,Velocity,Lx,Ly,E_points_4_adv_func/sub_ele_size,fileID_Tecplot_PowerIntensity_UOT_sub);
    Tecplot_Matlab_creator_2D (ss,N_T,xnode_T,ynode_T,Velocity,Lx_T,Ly_T,E_points_4_adv_func_T/Tape_ele_size,fileID_Tecplot_PowerIntensity_UOT_Tape);
%      fprintf(fileID_Tecplot_Roller_Q_Temp,' %12.8f , %12.8f  \r\n',Q_roller(ss),Temp_roller_current);


    end
% ss=1;

%       disp max intensity
%     sub_ele_size=delta_x*delta_y;
% Tape_ele_size=delta_x_T*delta_y_T;
% 
% Intensity_max=([max(E_points_4_adv_func/sub_ele_size),max(E_points_4_adv_func_T/Tape_ele_size)]);
%     str=sprintf('%f %f', Intensity_max)

% in order to see effect on the nip-point in advance !!

T_amb=20;

[T_matrix_Sub_Tape,T_old, T_old_T]= UOT_FDM_2D_3D_transient_TapeANDSub_combined(Lx,Ly,Lz,xnode,ynode,znode,E_points,Velocity,materials_sub,Temp_Right_sub(1),h_conv_Sub,...
    T_amb_mandrel,Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box_Sub,Graphic_chekbox(2),Graphic_chekbox(5),Graphic_chekbox(6),...
    Live_output,init_Temp,Inc_transient,Time_transient,...
    Lx_T,Ly_T,xnode_T,ynode_T,E_points_4_adv_func_T,Velocity_T,materials_Tape,Temp_roller_current,h_conv_T,thick_T,Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box_Tape,...
    Graphic_chekbox (3),Graphic_chekbox (4),Graphic_chekbox(2),...
    Total_energy, Video_transient,...
    manufacturing_type,ss,...
    T_old, T_old_T,fileID_Tecplot_Tape,fileID_Tecplot_Sub3D,...
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
    indexxy_T,delta_t,T_amb);


%Target
% Temp_sub= T_matrix_Sub_Tape{1}(:,end);
% Temp_tape= T_matrix_Sub_Tape{2}(:,1);
Temp_sub= T_matrix_Sub_Tape{1}(:,:);
Temp_tape= T_matrix_Sub_Tape{2}(:,:);


Power_Abs_Sub=sum(E_points);

Power_Abs_Tape=sum(E_points_T);


%     Temp_tape=T_old_T;
%    Temp_sub=  T_old;



%     To know each step temperature
% fileID15= fopen(sprintf('.\T_Sub_3D%d.txt',ss),'w');
% fprintf(fileID15,'xnode= %d, ynode=%d, znode=%d \r\n', xnode,ynode,znode );


% fileID15=fopen(strcat(UOT_pathfile, sprintf('T_Sub_3D_%d.txt',ss)),'w');



% fprintf(fileID15,'%f \r\n',T_old );
%
%
% fileID16=fopen(strcat(UOT_pathfile, sprintf('T_Tape_%d.txt',ss)),'w');
%
% fprintf(fileID16,'%f \r\n',T_old_T );
%
%    fclose(fileID15);
%  fclose(fileID16);




% end


% if Tecplot_check
%
%     fclose(fileID_Tecplot_Sub3D);
%     fclose(fileID_Tecplot_Roller_Q_Temp);
%     fclose(fileID_Tecplot_PowerIntensity_UOT_sub);
%     fclose(fileID_Tecplot_PowerIntensity_UOT_Tape);
%
% end
%
% if Tecplot_check_T
%
%     fclose(fileID_Tecplot_Tape);
% end


%
% if Graphic_chekbox (1)
%     fclose (fileID_Temp_Red_Box_Sub);
%     fclose(fileID_Temp_Red_Box_Tape);
% end


%% These two function for the case when sensor parameters are changing





