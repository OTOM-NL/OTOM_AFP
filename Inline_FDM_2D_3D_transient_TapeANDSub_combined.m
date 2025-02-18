

function T_matrix=Inline_FDM_2D_3D_transient_TapeANDSub_combined(Lx,Ly,Lz,xnode,ynode,znode,E_points,materials,Temp_Incoming,h,Temp_Conv_inside,...
    Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box,Tecplot_check,Graphic_profile,Grpahic_contour,...
    Live_output,init_Temp,Inc_transient,Time_transient,...
    Lx_T,Ly_T,xnode_T,ynode_T,E_points_T,materials_T,Temp_Right_T_Roller,h_T,thickness_T,...
    Row_number_Tape_Roller_tangent_T,Measure_Box_Tape_T,fileID_Temp_Red_Box_T,Graphic_Profile_T,Graphic_contour_T,Tecplot_check_T,...
    video_ID,manufacturing_type,ADS_Vars,TwinCat_file_name,...
    W_R,R_tape)


%8 April 2019
T_amb=20;

T_roller= Temp_Right_T_Roller(2);


% modified on 5 April 2019
Rel_v_theta=0;

% This function
% 1- recieve variable from TwinCat  >> input collection
% 2-Compute the temperature
% 3- Send back variables and new Temperature to TwinCat >> feedback section

asm = NET.addAssembly(strcat(pwd,'\TwinCAT.Ads.dll'));
import TwinCAT.Ads.*;
% Create instance of class TcAdsClient
adsClt=TwinCAT.Ads.TcAdsClient;

% tcClient.Connect(851)
adsClt.Connect(851);

% 1 variable to indicate the process is running/Stop

err=[];

[m,n]=size(ADS_Vars);
readValue=zeros(m,1);
hVar=zeros(m,1);
Data_type=cell(m,1);

% for saving handles and initialization

Active_Var=0;
for ii=1:m
    
    if ADS_Vars{ii}
        
        Active_Var=Active_Var+1;
        try
            
            hVar(ii) = adsClt.CreateVariableHandle( strcat(TwinCat_file_name,'.',ADS_Vars{ii,2})); %name of the variable
        catch err
            adsClt.Dispose();
            %             disp(err.message);
            msgbox ('Error');
        end
        
        
        Var=adsClt.ReadSymbolInfo(strcat(TwinCat_file_name,'.',ADS_Vars{ii,2}));
        Data_type{ii}=char(Var.Datatype);
        
        
    else
        break
    end
    
end
%%

% check whether manufacturing_type{1} is existed

% manufacturing_type{1}=objects(2).Value;    % Use cooling
% manufacturing_type{2}=objects(8).String;  % Sub Initial
% manufacturing_type{3}=objects(9).String;  % Tape inital
% manufacturing_type{4}=objects(13).String;  % No locations

% 9th April 2019
if manufacturing_type{5}
    
else
manufacturing_type{5}='2100;1000';
    
end




% decision for using Colling Temp or not
Layer_number=1;  %0 for winding on mandrel surface
Start_Time= 0; % should be obtained from the last time of previous winding layer

No_record_location=str2double(manufacturing_type{4});

Continious_pipe_winding=1;

Roller_rho_cp=str2num(manufacturing_type{5});

%   fprintf(fileID_85,'Layer number= %d , Continious_pipe_winding=%d ,Time_transient=%d, Inc_transient=%d,delta_t=%d \r\n',Layer_number+1,Continious_pipe_winding,Time_transient,Inc_transient,delta_t );


Continious_pipe_winding=1;


% Temp_Cooling
%   fileID_91=fopen('.\cooling\Temp_Cooling.txt','r');

% fclose(fileID_91);
%% for cooling simulation
% if ~ fileID_91
%
%     V_ave=mean(v_original);
% %     load T_nip_Tape; %from  txt file Tape_Temp_Global
% %     load T_nip_Sub;  % from txt file Sub3D_Temp_Global
% %     NO_Lap=40;   % no function at the moment
%     H=10;
%
% %     fileID_79=fopen(manufacturing_type{1},'r');  % Sub
% %      fileID_78=fopen(manufacturing_type{2},'r') ;  % Tape
%
%           T_nip_Sub = dlmread(manufacturing_type{1}) ;
%      T_nip_Tape = dlmread(manufacturing_type{2}) ;
%
%     Temp_Conv_inside=20;
%
%      fileID_90=fopen('Analysis_log.txt','r');
%
%      Time_last_layer = textscan(fileID_90,'  %f  ','Delimiter',',','HeaderLines',2) ;
% Time_last_layer=cell2mat(Time_last_layer);
%     fclose(fileID_90);
%
%     Temp_Cooling_pipe(Layer_number,V_ave,T_nip_Tape,T_nip_Sub,R_cyl,materials, Lx,Ly,H,Wind_dir,Temp_Conv_inside,Time_last_layer,delay_layers,Time_transient)
%
% else
%
% %       Time_last_layer = textscan(fileID_91,'  %f  ','Delimiter',',','HeaderLines',1)
% %       fileID_91=fopen('.\cooling\Temp_Cooling.txt','r');
% % fclose(fileID_91);
% end




fileID_71 = fopen('.\Tape_Temp_Global.txt','w');
fileID_72 = fopen('.\Sub3D_Temp_Global.txt','w');

fileID_84 = fopen('.\Analysis_log.txt','r');
Record_Time = textscan(fileID_84,' %f ','Delimiter',',','HeaderLines',2);
Record_Time=cell2mat(Record_Time);
fclose(fileID_84);


fileID_85 = fopen('.\Analysis_log.txt','w');





if Tecplot_check
    fileID_Tecplot_Tape = fopen('.\Temp_2D_Tape_Transient.plt','w');
end
if Tecplot_check_T
    fileID_Tecplot_Sub3D = fopen('.\Temp_3D_Sub_Transient.plt','w');
end



%% Should be modified
Time=Time_transient; % from user
inc=Inc_transient;% from user
delta_t=10e-3;  % the frequency that information is recieved

% if length(Power_Tank) >1
%     if inc ~=length(Power_Tank)
%         warndlg('Increment and number of data should be matched !');
%         inc=length(Power_Tank);
%         delta_t=Time/inc;
%
%     end
% else
%     v_original(1:inc)=v_original;
%     v_original_T(1:inc)=v_original_T;
%     Power_Tank(1:inc)=Power_Tank;
% end


Rec_moment=round(inc/No_record_location); % when the data should be recorded


fprintf(fileID_85,'Layer number= %d , Continious_pipe_winding=%d ,Time_transient=%d, Inc_transient=%d,delta_t=%d \r\n',Layer_number+1,Continious_pipe_winding,Time_transient,Inc_transient,delta_t );
fprintf(fileID_85,'Time of Recording Data (location should be replaced later) \r\n');


%% Substrate pre-analysis
% load initial value
% load T_initial;

% if video_ID % from user
% Video_Sub = VideoWriter('Sub-3D.avi','Uncompressed AVI');
% video.FrameRate = 30;

% open(Video_Sub );
%
% Video_Tape = VideoWriter('Tape-2D.avi','Uncompressed AVI');
% % video.FrameRate = 30;
%
% open(Video_Tape);
% end


A=Lx*Ly;

k=  materials(1);
%density
ro= materials(2);
%sprecific heat capacity
cp= materials(3);
%thermal diffusivity

if length(materials)>3
    ky=materials(4);
    kz=materials(5);
else
    ky=k;
    kz=k;
    
end



alpha=(k/ro*cp);






delta_x=Lx/(xnode-0); % it was Lx/(xnode-1);
delta_y=Ly/(ynode-0);
delta_z=Lz/(znode-0);

N=ynode*xnode*znode;
N_plane=ynode*xnode;

advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode)

%% How the heat flux is implemented should be investigated !?





fileID1 = fopen('.\nodes.txt');

node = textscan(fileID1,' %*d %f %f %f ','Delimiter',',') ;
Nodes=cell2mat(node);

%% after layer 1 , it should be applicable if Continious_pipe_winding

%   fileID_91=fopen('.\cooling\Temp_Cooling.txt','r');
% Data_shit = textscan(fileID_91,' %d ') ;

%%  ********************** Should be modified  Temp_incoming
% if manufacturing_type{1}   % if the file exist and there is a data in it
%     %     Temp_Incoming=Temp_initial(Time,inc,Layer_number,Record_Time); % read cooling file generated in .\cooling\
% else
%     
%     %     Temp_Incoming(1:inc)= Temp_Incoming;
% end



%%

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







%%
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




g_dot=1/((xnode-0)*(ynode-0));
g_dot= g_dot/(delta_x*delta_y*delta_z);

g_dot=ones(N_plane,1)*g_dot;
g_dot([BC_L_Xleft_Zbot,BC_L_Xright_Zbot,BC_L_Ybot_Zbot,BC_L_Ytop_Zbot])=g_dot([BC_L_Xleft_Zbot,BC_L_Xright_Zbot,BC_L_Ybot_Zbot,BC_L_Ytop_Zbot])*2;
g_dot([ BC_P_Xleft_Ybot_Z_bot,    BC_P_Xleft_Ytop_Z_bot,BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ytop_Z_bot])=...
    g_dot([ BC_P_Xleft_Ybot_Z_bot,    BC_P_Xleft_Ytop_Z_bot,BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ytop_Z_bot])*4;

g_dot_plane=g_dot; %(1:N_plane);
%%


transient_term=ro*cp/delta_t;


if manufacturing_type{1}
    
    
    
    %     Reading This file
    %      manufacturing_type{3} =
    %  D:\Profiles\zaamia\Desktop\OTOM-Dev-29June2018\T_Final_Sub_3D.txt
    % fileID_Temp_final_sub = fopen('T_Final_Sub_3D.txt');
    % T_Final_Sub_3D = textscan(fileID_Temp_final_sub,' %f ','Delimiter',',');
    
    %  manufacturing_type{4} =
    %  D:\Profiles\zaamia\Desktop\OTOM-Dev-29June2018\T_Final_Tape.txt
    
    % Temp=zeros(1,N);
    % load('T_Final_Sub_3D.mat');
    % READ Temperature from last layer
    fileID_Temp_final_sub = fopen( manufacturing_type{2});
    T_Final_Sub_3D = textscan(fileID_Temp_final_sub,' %f ','Delimiter',',','HeaderLines',1);
    T_Final_Sub_3D=cell2mat(T_Final_Sub_3D);
    
    fclose(fileID_Temp_final_sub);
    
    
    
    % if T_Final_Sub_3D
    if length (T_Final_Sub_3D) ==N
        T=T_Final_Sub_3D;
    else
        T=init_Temp*ones(N,1);
    end
    % end
    
else
    T=init_Temp*ones(N,1);
    
    
end

% CM_orig=sparse(1:N,1:N,ones(1,N),N,N);





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
            
            CM_ind_C_i_jp1_k=[CM_ind_C_i_jp1_k; kk,kk+1];
            CM_ind_C_i_jm1_k=[CM_ind_C_i_jm1_k;kk+1,kk ];
            
        end
        
    end
    
    if kk <N-ynode+1
        if mod(kk,N_plane)~=0
            
            
            CM_ind_C_ip1_j_k=[CM_ind_C_ip1_j_k; kk,kk+ynode];
            CM_ind_C_im1_j_k=[CM_ind_C_im1_j_k;kk+ynode,kk ];
            
        else
            
            
            CM_ind_0=[CM_ind_0;  [ [(kk-ynode+1:kk) , (kk-ynode+1:kk)+ynode] ;
                [(kk-ynode+1:kk)+ynode,  (kk-ynode+1:kk)] ]' ];
            %                  CM_ind_0=[CM_ind_0;(kk-ynode+1:kk)+ynode,kk-ynode+1:kk ];
            
            
        end
    end
    
    if kk <N-(ynode*xnode)+1
        
        
        CM_ind_C_i_j_kp1=[CM_ind_C_i_j_kp1;kk,kk+(ynode*xnode)];
        CM_ind_C_i_j_km1=[CM_ind_C_i_j_km1;kk+(ynode*xnode),kk];
        
        
    end
    
end

CM_indexes=[ [1:N;1:N]'  ...
    ;CM_ind_C_i_j_km1; CM_ind_C_i_j_kp1;CM_ind_0;CM_ind_C_im1_j_k;CM_ind_C_ip1_j_k;CM_ind_C_i_jm1_k;CM_ind_C_i_jp1_k];


%%

convection_term=-h*A/(xnode*ynode);%*((delta_x*delta_y*delta_z)^2);

convection_term=convection_term/(delta_x*delta_y*delta_z);% because delta_z is constant



fig_inline_long_width=figure('Name','Live Tape & Substrate','NumberTitle','off','Units','normalized', 'position',[0.05 0.55 0.4 0.4 ]);  % transient data output representation
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

fig_inline_long_width.ToolBar = 'none';
fig_inline_long_width.MenuBar = 'none';
set(fig_inline_long_width,'color','0.8 0.8 0.8');

%   hold on;
p1=[];
p2=[];
p3=[];
p4=[];


subplot(2,3,1);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Sub-Up-length','Color','k');
xlabel('Length');
ylabel('Temperature');

hold(ax);


subplot(2,3,2);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Sub-Bot-Le','Color','k');
xlabel('Length');
ylabel('Temperature');
hold(ax);

subplot(2,3,3);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Tape-Le','Color','k');
xlabel('Length');
ylabel('Temperature');
hold(ax);
%

subplot(2,3,4);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Sub-Up-Nip','Color','k');
xlabel('Width');
ylabel('Temperature');
hold(ax);


subplot(2,3,5);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Sub-Bot-Nip','Color','k');
xlabel('Width');
ylabel('Temperature');
hold(ax);


subplot(2,3,6);
ax=gca;
% set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'k';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Tape-Nip','Color','k');
xlabel('Width');
ylabel('Temperature');
hold(ax);


% Creat a GIF
axes('Position' , [0.1 0.92 0.05 0.05]);
gifplayer;

title('Substrate & Tape');


%% Tape pre-analysis
% Video = VideoWriter('animation-2D.avi','Uncompressed AVI');
% video.FrameRate = 60;
% %
% open(Video);


k_T=  materials_T(1);
%density
ro_T= materials_T(2);
%sprecific heat capacity
cp_T= materials_T(3);



if length(materials_T)>3
    ky_T=materials_T(4);
    kz_T=materials_T(5);
else
    ky_T=k_T;
    kz_T=k_T;
    
end



%thermal diffusivity
alpha_T=(k_T/ro_T*cp_T);
A_T=Lx_T*Ly_T;

delta_x_T=Lx_T/xnode_T;
delta_y_T=Ly_T/ynode_T;

N_T=ynode_T*xnode_T;



%   Temp_Conv_inside=25;  % Temperature for convection from surface, not edge!!
T_right_T=0;
T_left_T=Temp_Right_T_Roller(1);  % incoming velocity
T_top_T=0;
T_bottom_T=0;

T_amb_T=20; % out of plane temperature


BC_Left_T=1:ynode_T;
BC_top_T=ynode_T:ynode_T:N_T;
BC_bottom_T=1:ynode_T:N_T-1;
BC_right_T=N_T:-1:N_T-ynode_T+1 ;

indexxy_T=[N_T,N_T-ynode_T+1, 1 ynode_T];



for kk=1:xnode_T-1
    index_middle_long_T(kk)=floor(ynode_T/2)+(ynode_T*kk);
end
% Define a nip-point line
row_nip_T=0; % number before the last row
index_nip_point_T=(ynode_T*(xnode_T-row_nip_T-1)+1):((xnode_T-row_nip_T)*(ynode_T));

% Heat flux is constant here
g_dot_T=1/thickness_T;
g_dot_T=g_dot_T/(delta_x_T*delta_y_T);  % for the lement size
g_dot_T=g_dot_T/(xnode_T*ynode_T);  % total energy devided per number of nodes for each node

g_dot_T=ones(N_T,1)*g_dot_T;
g_dot_T([BC_bottom_T BC_top_T])=g_dot_T([BC_bottom_T BC_top_T])*2;
g_dot_T([BC_right_T BC_Left_T])=g_dot_T([BC_right_T BC_Left_T])*2;




transient_term_T=ro_T*cp_T/delta_t;


% load T_Final_Tape;

if manufacturing_type{1}
    %     Continious_pipe_winding && Layer_number
    
    fileID_T_Final_Tape = fopen(manufacturing_type{3});
    T_Final_Tape = textscan(fileID_T_Final_Tape,' %f ','Delimiter',',','HeaderLines',1);
    T_Final_Tape=cell2mat(T_Final_Tape);
    
    
    fclose(fileID_T_Final_Tape);
    
    % if T_Final_Sub_3D
    if length(T_Final_Tape ==N_T)
        Temp_T=T_Final_Tape;
    else
        Temp_T=init_Temp*ones(N_T,1);
    end
    
    
else
    Temp_T=init_Temp*ones(N_T,1);
    
    
end


% %      figure(1);
%      pos=get(gcf,'Position');
% figure(62);  % transient data output representation

%   hold on;
p1_T=[];
p2_T=[];
p3_T=[];
p4_T=[];



CM_T_C_i_jp1_T=[];
CM_T_C_i_jm1_T=[];
CM_T_C_ip1_j_T=[];
CM_T_C_im1_j_T=[];


for kk=1:N_T
    
    if mod(kk,ynode_T)~=0
        if kk <N_T
            
            
            CM_T_C_i_jp1_T=[CM_T_C_i_jp1_T; kk,kk+1];
            CM_T_C_i_jm1_T=[CM_T_C_i_jm1_T;kk+1,kk ];
        end
        
        %                RHS(kk,1)=-g_dot/delta_x/delta_y;
    end
    
    if kk <N_T-ynode_T+1
        
        CM_T_C_ip1_j_T=[CM_T_C_ip1_j_T; kk,kk+ynode_T];
        CM_T_C_im1_j_T=[CM_T_C_im1_j_T;kk+ynode_T,kk ];
        
        
    end
    
end

CM_indexes_T=[ [1:N_T;1:N_T]'  ...
    ;CM_T_C_im1_j_T; CM_T_C_ip1_j_T;CM_T_C_i_jm1_T;CM_T_C_i_jp1_T];



%%
fig_Inline_Monitoring=figure('Name','Live monitoring Vs Time','NumberTitle','off','Units','normalized', 'position',[0.05 0.05 0.4 0.4 ]);
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

% fig = gcf; % current figure handle
fig_Inline_Monitoring.ToolBar = 'none';
fig_Inline_Monitoring.MenuBar = 'none';
set(fig_Inline_Monitoring,'color','k');
% fig.Title='Live monitoring';



subplot(2,2,1);
ax=gca;
set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'w';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Input variable 1','Color','w');
hold(ax);
% ax.OuterPosition=[0 0 0.1 0.1];
% ax.ActivePositionProperty ='outerposition';

%

%

subplot(2,2,2);
ax=gca;
set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'w';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Output variable 1','Color','w');
hold(ax);

%

%

subplot(2,2,3);
ax=gca;
set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'w';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Input variable 2','Color','w');
hold(ax);
%

%


subplot(2,2,4);
ax=gca;
set(ax,'color',[0.2 0.2 0.1]);
ax.GridColor = 'w';
ax.XColor = 'y';
ax.YColor = 'y';
grid(ax,'on');
title('Output variable 2','Color','w');
hold(ax);

%

%


% Creat a GIF
axes('Position' , [0.475 0.92 0.05 0.05]);
gifplayer;

%
%   figure(79);

Div=40;

y=ones(1,Div)*20;
y2=y;

t=1:Div;
% figure;
subplot(2,2,1);
xlabel('Time');
ylabel('Temperature');
title('Tape');
PL_Tape=plot(t,y);

subplot(2,2,2);
xlabel('Time');
ylabel('Temperature');
title('Substrate');
PL_Sub=plot(t,y2);




Y_in1=zeros(1,Div);
Y_in2=zeros(1,Div);

subplot(2,2,3);
xlabel('Time');
ylabel('Velocity');
title('Inputs');
PL_Vel=plot(t,Y_in1);

subplot(2,2,4);
xlabel('Time');
ylabel('Power');
title('Inputs');
Pl_Pow=plot(t,Y_in2);

%% Pre-computation
% % >> Estimation of heat flux inside the Roller
rho_roller=Roller_rho_cp(1); %2180;
r_roller=R_tape;
w_roller=W_R;
cp_roller=Roller_rho_cp(2);
% time=delta_t;
fileID_Q_roller= fopen('.\Roller_ints.txt');
Q_roller = textscan(fileID_Q_roller,' %*f %*f %*f %f','Delimiter',',','HeaderLines',1);
Q_roller=cell2mat(Q_roller);
fclose(fileID_Q_roller);
Q=sum(Q_roller);
Temp_roller_current=T_roller;
Temp_roller_increase=Roller_temp_estimation(rho_roller,r_roller,w_roller,cp_roller,delta_t,Q);


%% Computations




counter=0;


p=1; % just to enter the loop

while p >0
    
    tic;
    
    for ii=1:Active_Var
        
        % for getting data Value from TwinCat
        
        switch  Data_type{ii}(6)
            case 'I'
                %  >>> These sentences to read value every time
                dataStreamRead = AdsStream(2); %number of y
                binRead = AdsBinaryReader(dataStreamRead);
                adsClt.Read(hVar(ii),dataStreamRead);
                readValue(ii) = binRead.ReadInt16();
                
                %                 msgbox(sprintf('Value= %f',readValue(ii)));
                
            case 'R'
                
                % for Real Datatype
                
                %  >>> These sentences to read value every time
                dataStreamRead = AdsStream(8); %number of y
                binRead = AdsBinaryReader(dataStreamRead);
                adsClt.Read(hVar(ii),dataStreamRead);
                readValue(ii) = binRead.ReadDouble();
                
                %                 msgbox(sprintf('Value= %f',readValue(ii)));
                
        end
        
    end
    
    p=readValue(1);  % still in the process or not
    v_original=readValue(2);
    Power_Tank=readValue(3);
    
    v_original_T=v_original;
    
    % for tt=1:inc
    
    
    
    T_old=T;
    RHS=zeros(N,1);
    %         RHS=sparse(N,1);
    
    
    v0=ro*cp*v_original;%/(delta_z^1);25*25*
    v=v0*abs(cosd(Rel_v_theta));
    v_2=v0*abs(sind(Rel_v_theta));
    pec=(v_original*delta_x)/alpha;
    
    RHS(1:N_plane)=-Power_Tank*E_points.*g_dot_plane*N_plane;
    
    diag_term=-2*k* (1/delta_x^2) - 2*ky*(1/delta_y^2)- 2*kz* (1/delta_z^2) -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    C_ip1_j_k=(k*(1/delta_x)^2);
    C_im1_j_k=(k*(1/delta_x)^2)+(v*(1/delta_x));
    C_i_jp1_k=(ky*(1/delta_y)^2);
    C_i_jm1_k=(ky*(1/delta_y)^2)+(v_2*(1/delta_y));
    C_i_j_kp1=(kz*(1/delta_z)^2);
    C_i_j_km1=(kz*(1/delta_z)^2);
    
    
    
    CM_Values=[ ones(N,1) *diag_term  ;
        ones(length(CM_ind_C_i_j_km1),1)*C_i_j_km1;...
        ones(length(CM_ind_C_i_j_kp1),1)*C_i_j_kp1;...
        zeros(length(CM_ind_0),1);...
        ones(length(CM_ind_C_im1_j_k),1)*C_im1_j_k;...
        ones(length(CM_ind_C_ip1_j_k),1)*C_ip1_j_k;...
        ones(length(CM_ind_C_i_jm1_k),1)*C_i_jm1_k;...
        ones(length(CM_ind_C_i_jp1_k),1)*C_i_jp1_k;...
        ];
    
    
    CM=sparse(CM_indexes(:,1),CM_indexes(:,2),CM_Values);
    
    
    %% for continous boundary x
    % for example in right side of the cube >> T_right =0
    % and T_i+1,j =0, and T_i,j=-1
    % based on this formulation, T_i+1,j is automatically will be zero if
    % the corresponding ghost node is zero
    %       diag_term_CB_=-1*alpha*((1/delta_x^2)+(2/delta_y^2)+(2/delta_z^2))-(v/delta_x);
    %     Cx_R >> continious in x-dir in right side
    
    %         T_left=0;
    %         %     diag_term_Cx_R=-1*k* ( (delta_y*delta_z)^2 + 2*(delta_x*delta_z)^2 + 2*(delta_x*delta_y)^2 )-(v*delta_x*(delta_y*delta_z)^2);
    %         diag_term_Cx_R=-1*k* ( (1/delta_x)^2 +2*(1/delta_y)^2+ 2*(1/delta_z)^2 )-(v*(1/delta_x));
    %
    %         for ii=1:length(BC_Left)
    %             CM(BC_Left(ii),BC_Left(ii))=diag_term_Cx_R;  % only for the right side
    %
    %         end
    
    
    %%  %% for continous boundary z
    %     T_top_Z=0;
    %     T_bot_Z=0;
    %
    %     diag_term_Cz_R=-k*2*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    %
    %
    %     for ii=1:length(BC_bot_Z)
    %         CM(BC_bot_Z(ii),BC_bot_Z(ii))=diag_term_Cz_R;  % only for the right side
    %
    %         CM(BC_top_Z(ii),BC_top_Z(ii))=diag_term_Cz_R;  % only for the right side
    %
    %     end
    
    
    T_top_Z=0;
    T_bot_Z=0;
    
    diag_term_Cz_R=-k*2*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    
    for ii=1:length(BC_bot_Z)
        CM(BC_bot_Z(ii),BC_bot_Z(ii))=diag_term_Cz_R;  % only for the right side
        
        CM(BC_top_Z(ii),BC_top_Z(ii))=diag_term_Cz_R;  % only for the right side
        
    end
    
    
    
    diag_term_Cx_R=-k*1*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    diag_term_Cy_R=-k*2*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    for ii=1:length(BC_Left)
        CM(BC_Left(ii),BC_Left(ii))=diag_term_Cx_R;  % only for the right side
        CM(BC_right(ii),BC_right(ii))=diag_term_Cx_R;  % only for the right side
        
    end
    % ALONG y-axis
    for ii=1:length(BC_top)
        CM(BC_top(ii),BC_top(ii))=diag_term_Cy_R;  % only for the right side
        CM(BC_bottom(ii),BC_bottom(ii))=diag_term_Cy_R;  % only for the right side
        
    end
    
    
    %         along lines
    
    diag_term_Cxz=-k*1*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    for ii=1:length(BC_L_Xleft_Zbot)
        CM(BC_L_Xleft_Zbot(ii),BC_L_Xleft_Zbot(ii))=diag_term_Cxz;
        CM(BC_L_Xleft_Ztop(ii),BC_L_Xleft_Ztop(ii))=diag_term_Cxz;
        CM(BC_L_Xright_Zbot(ii),BC_L_Xright_Zbot(ii))=diag_term_Cxz;
        CM(BC_L_Xright_Ztop(ii),BC_L_Xright_Ztop(ii))=diag_term_Cxz;
        
    end
    
    
    diag_term_Cyz=-k*2*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    for ii=1:length(BC_L_Ybot_Zbot)
        CM(BC_L_Ybot_Zbot(ii),BC_L_Ybot_Zbot(ii))=diag_term_Cyz;
        CM(BC_L_Ytop_Zbot(ii),BC_L_Ytop_Zbot(ii))=diag_term_Cyz;
        CM(BC_L_Ybot_Ztop(ii),BC_L_Ybot_Ztop(ii))=diag_term_Cyz;
        CM(BC_L_Ytop_Ztop(ii),BC_L_Ytop_Ztop(ii))=diag_term_Cyz;
        
    end
    
    diag_term_Cxy=-k*1*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    for ii=1:length(BC_L_Xleft_Ybot)
        CM(BC_L_Xleft_Ybot(ii),BC_L_Xleft_Ybot(ii))=diag_term_Cxy;
        CM(BC_L_Xleft_Ytop(ii),BC_L_Xleft_Ytop(ii))=diag_term_Cxy;
        CM(BC_L_Xright_Ybot(ii),BC_L_Xright_Ybot(ii))=diag_term_Cxy;
        CM(BC_L_Xright_Ytop(ii),BC_L_Xright_Ytop(ii))=diag_term_Cxy;
        
    end
    
    
    % points BC
    
    
    diag_term_Cxyz=-k*1*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
    
    CM(BC_P_Xleft_Ybot_Z_bot,BC_P_Xleft_Ybot_Z_bot)=diag_term_Cxyz;
    CM(BC_P_Xleft_Ytop_Z_bot,BC_P_Xleft_Ytop_Z_bot)=diag_term_Cxyz;
    CM(BC_P_Xleft_Ybot_Z_top,BC_P_Xleft_Ybot_Z_top)=diag_term_Cxyz;
    CM(BC_P_Xleft_Ytop_Z_top,BC_P_Xleft_Ytop_Z_top)=diag_term_Cxyz;
    
    
    CM(BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ybot_Z_bot)=diag_term_Cxyz;
    CM(BC_P_Xright_Ytop_Z_bot,BC_P_Xright_Ytop_Z_bot)=diag_term_Cxyz;
    CM(BC_P_Xright_Ybot_Z_top,BC_P_Xright_Ybot_Z_top)=diag_term_Cxyz;
    CM(BC_P_Xright_Ytop_Z_top,BC_P_Xright_Ytop_Z_top)=diag_term_Cxyz;
    
    
    
    
    % convection air upprt and bottom surface
    
    Conv_inside=convection_term*Temp_Conv_inside;
    
     Conv_amb=convection_term*T_amb;
    RHS(1:N_plane)=RHS(1:N_plane)+ Conv_amb;
    % Create separate index and value arrays.
    % Call sparse to assemble the index and value arrays.
    
    for hh=(znode-1)*N_plane+1:N    % for convection of Top_z surface
        CM(hh,hh)=CM(hh,hh)+convection_term;
        RHS(hh)=RHS(hh)+ Conv_inside;
    end
    
    
    for ii=1:length(CM_ind_0)
        CM(CM_ind_0(ii,1),CM_ind_0(ii,2))=0;
    end
    
    
    % Note: The effect of continoius boundary is not that much in case of
    % high velocity
    
    % T_left should be matrix of the incloming conditions Temperature T(BC_Left)
    
    
    
    %     RHS(BC_right)=RHS(BC_right)-((k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    %     RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
    %
    %     RHS(BC_bottom)=RHS(BC_bottom)-((ky*(1/delta_y)^2)+(v_2*(1/delta_y)))*T_bottom;
    %     RHS(BC_top)=RHS(BC_top)-(ky*(1/delta_y)^2)*T_top;
    %
    %     RHS(BC_top_Z)=RHS(BC_top_Z)-(kz*(1/delta_z)^2)*T_top_Z;
    %     RHS(BC_bot_Z)=RHS(BC_bot_Z)-(kz*(1/delta_z)^2)*T_bot_Z;
    %
    RHS(BC_right)=RHS(BC_right)-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    
    % 2*2 lines
    RHS([BC_L_Xright_Zbot BC_L_Xright_Ztop])=RHS([BC_L_Xright_Zbot BC_L_Xright_Ztop])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    RHS([BC_L_Xright_Ybot, BC_L_Xright_Ytop])=RHS([BC_L_Xright_Ybot, BC_L_Xright_Ytop])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    % 4 corner points
    RHS([BC_P_Xright_Ybot_Z_bot,  BC_P_Xright_Ytop_Z_bot, BC_P_Xright_Ybot_Z_top, BC_P_Xright_Ytop_Z_top])=RHS([BC_P_Xright_Ybot_Z_bot,  BC_P_Xright_Ytop_Z_bot, BC_P_Xright_Ybot_Z_top, BC_P_Xright_Ytop_Z_top])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    
    RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
    
    RHS(BC_bottom)=RHS(BC_bottom)-((ky*(1/delta_y)^2)+(v_2*(1/delta_y)))*T_bottom;
    RHS(BC_top)=RHS(BC_top)-(ky*(1/delta_y)^2)*T_top;
    
    RHS(BC_top_Z)=RHS(BC_top_Z)-(kz*(1/delta_z)^2)*T_top_Z;
    RHS(BC_bot_Z)=RHS(BC_bot_Z)-(kz*(1/delta_z)^2)*T_bot_Z;
    
    
    
    RHS=RHS - (transient_term*T_old);
    
    Temp= (RHS'/CM) ;
    
    %    [L,U,P] = lu(CM);
    % T= U\(L\(P*RHS));
    
    % [L,U,p] = lu(CM,'vector');
    % T = U\(L\(RHS(p,:)));
    
    T=full(Temp)' ;
    
    
    %% for extracting the measured box
    if fileID_Temp_Red_Box
        Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box_Sub,delta_x,delta_y,index_middle_long);
        fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
    end
    
    
    %% Matlab Tecplot creator animation
    %     if Tecplot_check
    %         Tecplot_Matlab_creator (tt,N,xnode,ynode,znode,v_original,Lx,Ly,Lz,T,fileID_Tecplot_Sub3D);
    %     end
    
    %%
    if Live_output
        figure(fig_inline_long_width);  % transient data output representation
        
        delete([p1 p2 p3 p4]);
        % along long and width
        %    title(sprintf('Time = %f',(tt/inc)*Time));
        subplot(2,3,1);
        %     hold on;
        
        % title('Tmepereture along middle of the tape');
        %     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
        %     ylabel('Temp');
        p1=plot(linspace(0,Lx,length(index_middle_long)),T(index_middle_long),'b');
        %     legend ('Sub-Upper-length');
        subplot(2,3,2)
        p2=plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
        %     legend ('Sub-Bottom-length');
        %     drawnow;
        
        
        %     figure(31);
        subplot(2,3,4);
        %     hold on;
        
        %     ylabel('Temp');
        % title('Tmepereture along width nip-point of the tape');
        %     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
        p3=plot(linspace(0,Ly,length(index_nip_point)),T(index_nip_point),'b');
        %     legend ('Sub-Upper-Nip');
        subplot(2,3,5);
        p4=plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
        %     legend ('Sub-Bottom-Nip');
        %     drawnow;
        
    end
    
    
    %% Tape computations
    
    T_old_T=Temp_T;
    CM_T=sparse(N_T);
    %    CM_T=zeros(N);
    
    
    
    % RHS=zeros(N,1);
    RHS_T=sparse(N_T,1);
    RHS_T(1:end)=-Power_Tank*E_points_T.*g_dot_T *N_T; % /N because of intensity input power come as a matrix
    
    v_T=ro_T*cp_T*v_original_T;
    pec_T=(v_original_T*delta_x_T)/alpha_T;
    
    
    
    h1_T=h_T(1);
    h2_T=h_T(2);
    
    convection_term_T=-h1_T*A_T/(xnode_T*ynode_T);  % A/(xnode*ynode); is the area of each element !
    convection_term_T=convection_term_T/(delta_x_T*delta_y_T*thickness_T);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3
    
    convection_term_Roller_T=-h2_T*A_T/(xnode_T*ynode_T);  % A/(xnode*ynode); is the area of each element !
    convection_term_Roller_T=convection_term_Roller_T/(delta_x_T*delta_y_T*thickness_T);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3
    
    diag_term_T=-2*k_T*(1/delta_x_T^2)-2*ky_T*(1/delta_y_T^2)-(v_T/delta_x_T)+convection_term_T-transient_term_T;
    diag_term_Roller_contact_T=-2*k_T*(1/delta_x_T^2)-2*ky_T*(1/delta_y_T^2)-(v_T/delta_x_T)+convection_term_Roller_T-transient_term_T;
    index_Roller_con_T=1+(Row_number_Tape_Roller_tangent_T*ynode_T);
    
    C_ip1_j_T=(k_T/delta_x_T^2);
    C_im1_j_T=(k_T/delta_x_T^2)+(v_T/delta_x_T);
    C_i_jp1_T=(ky_T/delta_y_T^2);
    C_i_jm1_T=(ky_T/delta_y_T^2);
    
    CM_Values_T=[ ones(N_T,1) *diag_term_T  ;
        ones(length(CM_T_C_im1_j_T),1)*C_im1_j_T;...
        ones(length(CM_T_C_ip1_j_T),1)*C_ip1_j_T;...
        ones(length(CM_T_C_i_jm1_T),1)*C_i_jm1_T;...
        ones(length(CM_T_C_i_jp1_T),1)*C_i_jp1_T;...
        ];
    
    
    CM_T=sparse(CM_indexes_T(:,1),CM_indexes_T(:,2),CM_Values_T);
    
    
    for ii=index_Roller_con_T:N_T
        CM_T(ii,ii)=diag_term_Roller_contact_T;
    end
    
    
    %% for continous boundary x
    
    
    %     diag_term_CB_xR_T=-1*k_T*(1/delta_x_T^2)-1*ky_T*(2/delta_y_T^2) -(v_T/delta_x_T)+convection_term_Roller_T-transient_term_T;
    %
    %     T_right_T=0;
    %
    %     for ii=1:length(BC_right_T)
    %
    %         CM_T(BC_right_T(ii),BC_right_T(ii))=diag_term_CB_xR_T;  % for right side
    %
    %     end
    
    %% for continous boundary
    diag_term_CB_y_T=-2*k_T*(1/delta_x_T^2)-1*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_T-transient_term_T;
    diag_term_CB_y_T_after_roller=-2*k_T*(1/delta_x_T^2)-1*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_Roller_T-transient_term_T;
    
    for ii=1:length(BC_top_T)
        if BC_top_T(ii) <index_Roller_con_T
            CM_T(BC_top_T(ii),BC_top_T(ii))=diag_term_CB_y_T;
            CM_T(BC_bottom_T(ii),BC_bottom_T(ii))=diag_term_CB_y_T;
            
        else
            CM_T(BC_top_T(ii),BC_top_T(ii))=diag_term_CB_y_T_after_roller;
            CM_T(BC_bottom_T(ii),BC_bottom_T(ii))=diag_term_CB_y_T_after_roller;
            
            
        end
        
    end
    
    diag_term_CB_xR_T_right=-1*k_T*(1/delta_x_T^2)-2*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_Roller_T-transient_term_T;
    diag_term_CB_xR_T_left=-1*k_T*(1/delta_x_T^2)-2*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_T-transient_term_T;
    
    %     T_right_T=0;
    
    for ii=1:length(BC_right_T)
        
        CM_T(BC_right_T(ii),BC_right_T(ii))=diag_term_CB_xR_T_right;  % for right side
        CM_T(BC_Left_T(ii),BC_Left_T(ii))=diag_term_CB_xR_T_left;  % for right side
        
    end
    
    %     for four corners
    diag_term_CB_xy_T=-1*k_T*(1/delta_x_T^2)-1*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_T-transient_term_T;
    diag_term_CB_xy_T_after_roller=-1*k_T*(1/delta_x_T^2)-1*ky_T*(1/delta_y_T^2) -(v_T/delta_x_T)+convection_term_Roller_T-transient_term_T;
    %     for four corners
    for ii=1:2
        %          index_xy=indexxy(ii);
        CM_T(indexxy_T(ii),indexxy_T(ii))=diag_term_CB_xy_T_after_roller;  %  ??
    end
    for ii=3:4
        %          index_xy=indexxy(ii);
        CM_T(indexxy_T(ii),indexxy_T(ii))=diag_term_CB_xy_T;  %  ??
    end
    
    
    
    
    
    
    % This Should be modified! >> modified on 9 April 2019
    %   Temp_roller_increase_inc=T_roller;
    Temp_roller_current=Temp_roller_increase +Temp_roller_current;
    
    
    
    RHS_T(BC_Left_T)=RHS_T(BC_Left_T)-((0*k_T/delta_x_T^2)+(v_T/delta_x_T))*T_left_T;
    RHS_T(1:index_Roller_con_T-1)=RHS_T(1:index_Roller_con_T-1)+convection_term_T*T_amb_T;
    
    RHS_T(index_Roller_con_T:N_T)=RHS_T(index_Roller_con_T:N_T)+convection_term_Roller_T*Temp_roller_current;
    
    
    %      spy(CM);
    % implementing BC
    
    
    %     RHS_T(BC_Left_T)=RHS_T(BC_Left_T)-((k_T/delta_x_T^2)+(v_T/delta_x_T))*T_left_T;
    %     RHS_T(BC_bottom_T)=RHS_T(BC_bottom_T)-(ky_T/delta_y_T^2)*T_bottom_T;
    %     RHS_T(BC_top_T)=RHS_T(BC_top_T)-(ky_T/delta_y_T^2)*T_top_T;
    %     RHS_T(BC_right_T)=RHS_T(BC_right_T)-(k_T/delta_x_T^2)*T_right_T;
    
    %     RHS_T(1:index_Roller_con_T-1)=RHS_T(1:index_Roller_con_T-1)+convection_term_T*T_amb_T- (transient_term_T*T_old_T(1:index_Roller_con_T-1));
    
    %      RHS_T(index_Roller_con_T:N_T)=RHS_T(index_Roller_con_T:N_T)+convection_term_Roller_T*T_amb_T- (transient_term_T*T_old_T(index_Roller_con_T:N_T));
    
    %     RHS_T(1:index_Roller_con_T-1)=RHS_T(1:index_Roller_con_T-1)+convection_term_T*T_amb_T;
    
    
    %     RHS_T(index_Roller_con_T:N_T)=RHS_T(index_Roller_con_T:N_T)+convection_term_Roller_T*T_amb_T;
    
    RHS_T=RHS_T - (transient_term_T*T_old_T);
    
    
    
    %% ******** Solution ********
    %    Temp= RHS'/CM ;
    
    %       [L,U,P] = lu(CM);
    % Temp= U\(L\(P*RHS));
    
    [L,U,p] = lu(CM_T,'vector');
    Temp_T = U\(L\(RHS_T(p,:)));
    
    
    toc
    
    T_T=full(Temp_T);
    
    if fileID_Temp_Red_Box_T
        Temp_Red_Box_T=Measuring_Box_Cal (T_T,Measure_Box_Tape_T,delta_x_T,delta_y_T,index_middle_long_T);
        fprintf(fileID_Temp_Red_Box_T,'   %12.8f  \r\n',Temp_Red_Box_T);
    end
    
    
    %% For saving the Nip-point temperature matrix
    %        T_nip_Tape= T(BC_right);
    %       save('T_nip_Tape.mat','T_nip_Tape');
    
    
    %% Matlab Tecplot creator animation
    %     if Tecplot_check_T
    %         Tecplot_Matlab_creator_2D (tt,N_T,xnode_T,ynode_T,v_original_T,Lx_T,Ly_T,T_T,fileID_Tecplot_Tape);
    %     end
    %
    
    %%
    
    if Live_output
        figure(fig_inline_long_width);  % transient data output representation
        
        delete([p1_T  p3_T ]);
        % along long and width
        %    title(sprintf('Time = %f',(tt/inc)*Time));
        subplot(2,3,3);
        %     hold on;
        
        % title('Tmepereture along middle of the tape');
        %     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
        %     ylabel('Temp');
        p1_T=plot(linspace(0,Lx_T,length(index_middle_long_T)),T_T(index_middle_long_T),'g');
        %     legend ('Tape along Length');
        %     subplot(2,2,2)
        %     p2=plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
        %     legend ('Bottom');
        %     drawnow;
        
        
        %     figure(31);
        subplot(2,3,6);
        %     hold on;
        
        %     ylabel('Temp');
        % title('Tmepereture along width nip-point of the tape');
        %     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
        p3_T=plot(linspace(0,Ly_T,length(index_nip_point_T)),T_T(index_nip_point_T),'g');
        %     legend ('Tape along nip-point');
        %     subplot(2,2,4);
        %     p4=plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
        %     legend ('Bottom');
        drawnow;
        
        
    end
    
    
    
    %    pause(0.005);
    % %     Mov=getframe;
    %  Mov=getframe(gcf);
    
    
    
    
    %% For saving the Nip-point temperature matrix
    % >>>>>>>>>>  It should be able to save when there is a smoke, high pick in
    % Temperature.. important to estimate cooling on that specific
    % location<<<<<<<<<<<<
    
    
    %     if Pick_determination
    %
    %         if Temp_Red_Box > max_Temp_point  % Temp of substrate
    %             max_Temp_point=Temp_Red_Box;
    %
    %             % keep maximum values
    %             T_nip_Sub_pick= T(BC_Left)' ;   %as a txt file
    %             T_nip_Tape_pick= T_T(BC_right_T)' ;
    %             Current_Time_pick=Start_Time+ (tt*delta_t);
    %
    %         else Temp_Red_Box < min_Temp_point;
    %             min_Temp_point=Temp_Red_Box;
    %         end
    %
    %
    %
    %
    %         slope(1:2)=slope(2:3);
    %
    %         if Temp_Red_Box > Temp_Red_Box_OLD
    %             slope(3)=+1;
    %         else
    %             slope(3)=-1;
    %         end
    %
    %         Temp_Red_Box_OLD=Temp_Red_Box;
    %
    %
    %
    %         % relative maximum
    %         if slope (2:3)==[+1 -1]
    %
    %             if (max_Temp_point-min_Temp_point) >Tol_Temp
    %                 fprintf(fileID_71,'%f  ',T_nip_Tape_pick );
    %                 fprintf(fileID_71,'\r\n ' );
    %
    %                 fprintf(fileID_72,'%f ',T_nip_Sub_pick );
    %                 fprintf(fileID_72,'\r\n ' );
    %
    %                 fprintf(fileID_85,'%d \r\n',Current_Time_pick );
    %
    %
    %                 slope=[0 0 0];
    %                 max_Temp_point=Temp_Red_Box;
    %                 min_Temp_point=Temp_Red_Box;
    %
    %             end
    %         elseif slope (1:3)==[+1 +1 +1]  % in order to forget last written data and focus on new data
    %             max_Temp_point=Temp_Red_Box_OLD;
    %
    %
    %         end
    %
    %
    %     end
    %
    
    
    counter=counter+1;
    
    if  mod(counter,Rec_moment)==0
        T_nip_Sub= T(BC_Left)' ;   %as a txt file
        T_nip_Tape= T_T(BC_right_T)' ;
        
        
        fprintf(fileID_71,'%f  ',T_nip_Tape );
        fprintf(fileID_71,'\r\n ' );
        
        fprintf(fileID_72,'%f ',T_nip_Sub );
        fprintf(fileID_72,'\r\n ' );
        
        %         fprintf(fileID_72,'%f \r\n', );
        
        
        Current_Time=Start_Time+ (counter*delta_t);
        
        fprintf(fileID_85,'%d \r\n',Current_Time );
        
    end
    
    
    %% Sending Feedback
    
    
    %% Writting in TwinCat variables using output of computational model
    
    % Writing putput temperature
    valueToWrite(4:5) = [Temp_Red_Box_T,Temp_Red_Box];
    
    for ii=4:Active_Var
        
        
        % for writting data Value from TwinCat
        % Create AdsStream instance
        
        switch  Data_type{ii}(6)
            %             case 'I'
            %                 % for integer value
            %                 %  >>> These sentences to write value every time
            %
            %                 valueToWrite = 3+ii;
            %                 dataStreamWrite= AdsStream(2);
            %                 binWrite = AdsBinaryWriter(dataStreamWrite);
            %                 binWrite.Write(int16(valueToWrite));
            %                 adsClt.Write(hVar(ii),dataStreamWrite);
            
            
            %             adsClt.DeleteVariableHandle(hVar);
            
            case 'R'
                % for Real Datatype
                %  >>> These sentences to write value every time
                
                
                
                dataStreamWrite= AdsStream(8);
                binWrite = AdsBinaryWriter(dataStreamWrite);
                binWrite.Write(double(valueToWrite(ii)));
                adsClt.Write(hVar(ii),dataStreamWrite);
                
                %             adsClt.DeleteVariableHandle(hVar);
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    %   figure(79);
    
    figure(fig_Inline_Monitoring);
    
    
    subplot(2,2,1);
    delete (PL_Tape);
    PL_Tape=plot(t,y,'w');
    %          drawnow;
    y=[y(2:end) Temp_Red_Box_T ];
    
    subplot(2,2,2);
    y2=[y2(2:end) Temp_Red_Box ];
    delete (PL_Sub);
    PL_Sub=plot(t,y2,'w');
    %       drawnow;
    
    
    subplot(2,2,3);
    Y_in1=[Y_in1(2:end) v_original ];
    delete (PL_Vel);
    PL_Vel=plot(t,Y_in1,'w');
    %       drawnow;
    
    
    subplot(2,2,4);
    Y_in2=[Y_in2(2:end) Power_Tank ];
    delete (Pl_Pow);
    Pl_Pow=plot(t,Y_in2,'w');
    drawnow;
    
    
    
    
    %     pause(0.01);
    
    
    
    
    
    
end

adsClt.Dispose();



%% put the legend
% if Live_output
%
%     %     delete([p1 p2 p3 p4]);
%     % along long and width
%     %    title(sprintf('Time = %f',(tt/inc)*Time));
%     subplot(2,3,1);
%
%     % title('Tmepereture along middle of the tape');
%
%     legend ('Sub-Upper-length');
%     subplot(2,3,2)
%     legend ('Sub-Bottom-length');
%
%
%     %     figure(31);
%     subplot(2,3,4);
%     legend ('Sub-Upper-Nip');
%     subplot(2,3,5);
%     legend ('Sub-Bottom-Nip');
%
%
%
%     subplot(2,3,3);
%     legend ('Tape along Length');
%
%     subplot(2,3,6);
%     legend ('Tape along nip-point');
%
%
%
%
% end


%%



%% For saving the Nip-point temperature matrix
% T_nip_Sub= T(BC_Left)';   %as a txt file
% save('T_nip_Sub.mat','T_nip_Sub');  % save in a txt file in loop
% save('T_Final_Sub_3D.mat','T');

%% For saving the Nip-point temperature matrix
% T_nip_Tape= T_T(BC_right_T);
% save('T_nip_Tape.mat','T_nip_Tape');
% save('T_Final_Tape.mat','T_T');

%    fprintf(fileID_71,'%f, ',T_nip_Tape );
%         fprintf(fileID_72,'%f, ',T_nip_Sub );

% fileID_71 = fopen(sprintf('Tape_Temp_Global.txt'),'w');
% fileID_72 = fopen(sprintf('Sub3D_Temp_Global.txt'),'w');

% fclose([fileID_71 ,fileID_72, fileID_85] );
fclose(fileID_71);
fclose(fileID_72);
fclose(fileID_85);



fileID15= fopen('.\T_Final_Sub_3D.txt','w');
fprintf(fileID15,'xnode= %d, ynode=%d, znode=%d \r\n', xnode,ynode,znode );
fprintf(fileID15,'%f \r\n',T );

fileID16= fopen('.\T_Final_Tape.txt','w');
fprintf(fileID16,'xnode= %d, ynode=%d \r\n', xnode_T, ynode_T );
fprintf(fileID16,'%f \r\n',T_T );


fclose(fileID15);
fclose(fileID16);

% heating surface
T_matrix=reshape(T(N_plane:-1:1),ynode,xnode);




if Grpahic_contour
    
    T_matrix_bott=reshape(T(N:-1:N-N_plane+1),ynode,xnode);
    
    
    [X_2D,Y_2D]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
    h46=figure(46);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h46,'Visible','off');
    
    subplot(2,1,1);
    surf(X_2D,Y_2D,T_matrix);
    shading flat;
    %         view([0 -90])
    %         Frames=getframe;
    %       hold on
    %    axis equal
    title(sprintf('Upper-S, velocity= %f , Peclet number: %f ',v_original,pec));
    axis([0 Lx 0 Ly]);
    colorbar;
    
    
    pos=get(gca,'Position');
    
    x = pos(1)+[pos(3)/4 3*pos(3)/4];
    y = pos(2)+[pos(4)/2 pos(4)/2];
    annotation('textarrow',x,y,'String','Velocity','Color','w','Linewidth',2);
    
    if abs (v_2) > 1
        x = pos(1)+[pos(3)/4 pos(3)/4];
        y = pos(2)+[pos(4)/2 pos(4)/4];
        annotation('textarrow',x,y,'String','Velocity-Y','Color','r','Linewidth',2);
    end
    
    view([0 -90]);
    
    %% for Bottom surface
    
    figure(46);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h46,'Visible','off');
    %         subplot(3,2,vv);
    subplot(2,1,2);
    surf(X_2D,Y_2D,T_matrix_bott);
    shading flat
    
    %         Frames=getframe;
    %       hold on
    %    axis equal
    title(sprintf('Bottom-S, velocity= %f , Peclet number: %f , Step=%d/%d',v_original,pec));
    axis([0 Lx 0 Ly]);
    colorbar;
    pos=get(gca,'Position');
    
    x = pos(1)+[pos(3)/4 3*pos(3)/4];
    y = pos(2)+[pos(4)/2 pos(4)/2];
    annotation('textarrow',x,y,'String','Velocity','Color','w','Linewidth',2);
    
    if abs (v_2) > 1
        x = pos(1)+[pos(3)/4 pos(3)/4];
        y = pos(2)+[pos(4)/2 pos(4)/4];
        annotation('textarrow',x,y,'String','Velocity-Y','Color','r','Linewidth',2);
    end
    
    
    view([0 -90]);
    
    
    
    %% 3D
    
    h36=figure(36);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h36,'Visible','off');
    
    X=zeros(ynode*xnode,znode);
    Y=zeros(ynode*xnode,znode);
    Z=zeros(ynode*xnode,znode);
    T_show=reshape(T,ynode*xnode,znode);
    
    for kk=1:ynode*xnode
        index=((1:znode)*xnode*ynode)-kk+1;
        X(kk,1:end)=Nodes(index,1)';
        Y(kk,1:end)=Nodes(index,2)';
        Z(kk,1:end)=Nodes(index,3)';
        %      Intensity(kk,1:end)=points(index,4)';
    end
    
    
    f_obj=mesh(X,Y,Z,T_show,'MarkerFaceColor','auto',...
        'Marker','.',...
        'LineStyle','none','FaceColor','none');
    
    shading flat
    % shading interp
    %      view([0 90])
    
    hold on
    %    axis equal
    title(sprintf('velocity= %f , Peclet number: %f ',v_original,pec));
    %      axis([0 Lx 0 Ly]);
    % view([-83.5 -52.4]);
    view([42.5 -34.8000000000001]);
    colorbar;
    % disp(vv);
    
    
    
    
end



%% along long and width

if Graphic_profile
    
    h56=figure(56);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h56,'Visible','off');
    
    
    %     subplot(3,2,vv);
    subplot(2,1,1);
    hold on;
    % title('Tmepereture along middle of the tape');
    title(sprintf('Temperature along length, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    ylabel('Temperature ^{\circ}C');
    xlabel('Length (m)');
    
    plot(linspace(0,Lx,length(index_middle_long)),T(index_middle_long),'b');
    plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
    legend ('Upper','Bottom');
    
    pos=get(gca,'Position');
    
    x = pos(1)+[pos(3)/4 3*pos(3)/4];
    y = pos(2)+[pos(4)/2 pos(4)/2];
    annotation('textarrow',x,y,'String','Velocity','Color','k','Linewidth',2)
    
    
    
    h56=figure(56);
    
    set(h56,'Visible','off');
    subplot(2,1,2);
    hold on;
    ylabel('Temperature ^{\circ}C');
    xlabel('Width (m)');
    % title('Tmepereture along width nip-point of the tape');
    title(sprintf('Temperature along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    plot(linspace(0,Ly,length(index_nip_point)),T(index_nip_point),'b');
    plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
    legend ('Upper','Bottom');
    
    
    
    T_nip_point=T(index_nip_point);
    Nip_point_Temp_Sub=T_nip_point(floor(length(index_nip_point)/2));
    assignin('base','Nip_point_Temp_Sub_3D',Nip_point_Temp_Sub);
    
end
%% Tape



if Graphic_Profile_T
    
    h21=figure(21);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h21,'Visible','off');
    
    subplot(2,2,3);
    hold on;
    % title('Tmepereture along middle of the tape');
    title(sprintf('Tape-Temperature along length, V= %f , Pec: %f, h=%d ',v_original_T,pec_T,h_T));
    title(sprintf('Tape-Temperature along length, V=%f , Pec=%f, h1,2=%d %d ',v_original_T,pec_T,h1_T,h2_T));
    
    ylabel('Temperature ^{\circ}C');
    xlabel('Length (m)');
    
    H_long=linspace(0,Lx_T,length(T_T(index_middle_long_T))); % H_long > horizental
    plot(H_long,T_T(index_middle_long_T));
    
    pos=get(gca,'Position');
    
    x = pos(1)+[pos(3)/4 pos(3)/2];
    y = pos(2)+[pos(4)/2 pos(4)/2];
    annotation('textarrow',x,y,'String','Velocity','Color','k','Linewidth',2)
    
    %     plot(T(index_middle_long));
    
    h21=figure(21);
    set(h21,'Visible','off');
    
    subplot(2,2,4);
    hold on;
    ylabel('Temperature ^{\circ}C');
    xlabel('Width (m)');
    % title('Tmepereture along width nip-point of the tape');
    title(sprintf('Tape-Temperature along width nip-point, V=%f , Pec=%f, h1,2=%d %d ',v_original_T,pec_T,h1_T,h2_T));
    
    T_nip_point_T=T_T(index_nip_point_T);
    H_Nip_T=linspace(0,Ly_T,length(T_nip_point_T)); % H_Nip > horizental
    plot(H_Nip_T,T_nip_point_T);
    Nip_point_Temp_Tape=T_nip_point_T(floor(length(index_nip_point_T)/2));
    
    assignin('base','Nip_point_Temp_Tape',Nip_point_Temp_Tape);
    
    %     plot(T(index_nip_point));
    
end

%% outPut 3D
if Graphic_contour_T
    h100=figure(100);
    javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
    set(h100,'Visible','off');
    
    
    T_matrix2=reshape(T_T,ynode_T,xnode_T);
    hold on;
    subplot(2,1,2)
    
    [X_2D,Y_2D]=meshgrid(linspace (0,Lx_T,xnode_T),linspace(0,Ly_T,ynode_T) );
    
    surf(X_2D,Y_2D,T_matrix2);
    title('Temperature distribution of Tape');
    colorbar;
    shading flat;
    
    
    pos=get(gca,'Position');
    x = pos(1)+[pos(3)/4 pos(3)/2];
    y = pos(2)+[pos(4)/2 pos(4)/2];
    annotation('textarrow',x,y,'String','Velocity','Color','w','Linewidth',2);
    
end

%%

len_T=length(T_T);
T_disp_T=zeros(len_T,1);


counter=0;



% N=xnode*ynode;
for ii=1:xnode_T
    for jj=ynode_T:-1:1
        counter=counter+1;
        index=N_T-(jj)+1-((ii-1)*ynode_T);
        T_disp_T(index)=T_T(counter);
        
    end
end



T_matrix_T=reshape(T_disp_T,ynode_T,xnode_T);

T_matrix={T_matrix, T_matrix_T}  ;   % for Substrate 3D and Tape



if Tecplot_check
    
    fclose(fileID_Tecplot_Sub3D);
end

if Tecplot_check_T
    
    fclose(fileID_Tecplot_Tape);
end


