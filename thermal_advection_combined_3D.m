%%The problem of 3D advection-diffusion case


function T_matrix=thermal_advection_combined_3D(Lx,Ly,Lz,xnode,ynode,znode,E_points,v_original,materials,Temp_Incoming,h,Temp_Conv_inside,...
       Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box,Tecplot_check,Graphic_profile,Grpahic_contour)


% $ relation between optical and thermal nodes > 
fid13 = fopen('.\Supp_files\edge_factors.txt');
out = textscan(fid13,'%s','delimiter',',');
    edge_factor_Sub=str2num(out{1,1}{4});
      fclose(fid13);
      
    
% edge_factor_Sub=1.5;





%  fileID_Tecplot_Sub3D = fopen(sprintf('Temp_3D-FDM_V.plt'),'w');
           
T_amb=20;
% format long;
% h=500;%1e5;%6; % W/m^2

% Thickness should be large enough to neglect effect of thickness, to be
% z-direction independent

% Lx=40e-3;
% Ly=20e-3;
% Lz=1*1.5e-4;

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
% Temp_Conv_inside=20;

% xnode=36;
% ynode=25;
% znode=6;
% Previous delta x,y,z
% delta_x=Lx/xnode;
% delta_y=Ly/ynode;
% delta_z=Lz/znode;
Lx=abs(Lx);

delta_x=Lx/(xnode-1); % it was Lx/(xnode-1);
delta_y=Ly/(ynode-1);
delta_z=Lz/(znode-0);

N=ynode*xnode*znode;
N_plane=ynode*xnode;

% For Tecplot output
advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode)

%% How the heat flux is implemented should be investigated !?

% Total_power=450;   %W



T_BC=sparse(N,1);

fileID1 = fopen('.\nodes.txt');

node = textscan(fileID1,' %*d %f %f %f ','Delimiter',',') ;
Nodes=cell2mat(node);

T_right=0;
T_left=Temp_Incoming;
T_top=0;
T_bottom=0;
% T_top_Z=20;
% T_bot_Z=20;

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











 g_dot=1/((xnode-0)*(ynode-0)*(znode/znode-0));%*A; %20e-2;%/(ro*cp);
% g_dot=g_dot*(delta_x*delta_y*delta_z)^2;
g_dot= g_dot/(delta_x*delta_y*delta_z);
% g_dot=g_dot*(1/delta_z);

g_dot=ones(N_plane,1)*g_dot;
g_dot([BC_L_Xleft_Zbot,BC_L_Xright_Zbot,BC_L_Ybot_Zbot,BC_L_Ytop_Zbot])=g_dot([BC_L_Xleft_Zbot,BC_L_Xright_Zbot,BC_L_Ybot_Zbot,BC_L_Ytop_Zbot])*edge_factor_Sub;
g_dot([ BC_P_Xleft_Ybot_Z_bot,    BC_P_Xleft_Ytop_Z_bot,BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ytop_Z_bot])=...
g_dot([ BC_P_Xleft_Ybot_Z_bot,    BC_P_Xleft_Ytop_Z_bot,BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ytop_Z_bot])*(edge_factor_Sub^2);


g_dot_plane=g_dot; %(1:N_plane);

%%
counter=0;
for kk=xnode-1:-1:0
    counter=counter+1;
        index_middle_long(counter)=ceil(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=1; % number before the last row % 2D and 3D are diiferent from each other
index_nip_point=(ynode*(row_nip-1)+1):((row_nip)*(ynode));


index_middle_long_Bott=index_middle_long +(N-N_plane);
index_nip_point_Bott=index_nip_point+(N-N_plane);

%%




    
%     tic

    v0=ro*cp*v_original;%/(delta_z^1);25*25*
    
    v=v0*abs(cosd(Rel_v_theta));
    v_2=v0*abs(sind(Rel_v_theta));
    
    
    
    %      fileID6 = fopen(sprintf('Temp_3D-FDM_V_0to%f.plt',v_original),'a');

        CM=sparse(N);
        %         CM=zeros(N);
        
        for kk=1:N
            CM(kk,kk)=1;
            
        end
        
        
        % T=sparse(1,N);
        
        % T=zeros(1,N);
        RHS=zeros(N,1);
%         RHS=sparse(N,1);


        RHS(1:N_plane)=-E_points.*g_dot_plane*N_plane;  
     
        % RHS=zeros(N,1);
        
        pec=(v_original*delta_x)/alpha;
        
        %         convection_term=-h*A*(delta_x*delta_y*delta_z)^2; % non-modified
        
        convection_term=-h*A/(xnode*ynode);%*((delta_x*delta_y*delta_z)^2);
        %     convection_term=convection_term/(A*delta_z);% because delta_z is constant
        %
        convection_term=convection_term/(delta_x*delta_y*delta_z);% because delta_z is constant
        %           convection_term=convection_term*Lz/ (delta_z);% because delta_z is constant
        
        %         diag_term=-2*k* ( (delta_y*delta_z)^2 +(delta_x*delta_z)^2 + (delta_x*delta_y)^2 )-(v*delta_x*(delta_y*delta_z)^2);
%         diag_term=-2*k* ( (1/delta_x^2) +(1/delta_y^2)+ (1/delta_z^2) )-(v*(1/delta_x))-(v_2*(1/delta_y));
        
         diag_term=-2*k* (1/delta_x^2) - 2*ky*(1/delta_y^2)- 2*kz* (1/delta_z^2) -(v*(1/delta_x))-(v_2*(1/delta_y));
        
        CM= CM*diag_term;
        
        
 %%
 
 
         
        C_ip1_j_k=(k*(1/delta_x)^2);
        C_im1_j_k=(k*(1/delta_x)^2)+(v*(1/delta_x));
        C_i_jp1_k=(ky*(1/delta_y)^2);
        C_i_jm1_k=(ky*(1/delta_y)^2)+(v_2*(1/delta_y));
        
        C_i_j_kp1=(kz*(1/delta_z)^2);
        C_i_j_km1=(kz*(1/delta_z)^2);
        
        
        
        for kk=1:N
            
            
            if mod(kk,ynode)~=0
                if kk <N
                    CM(kk,kk+1)=C_i_jp1_k;
                    CM(kk+1,kk)=C_i_jm1_k;
                end
                
            end
            
            if kk <N-ynode+1
                if mod(kk,N_plane)~=0
                    
                    CM(kk,kk+ynode)=C_ip1_j_k;
                    CM(kk+ynode,kk)=C_im1_j_k;
                else
                    CM(kk-ynode+1:kk,(kk-ynode+1:kk)+ynode)=0;
                    CM((kk-ynode+1:kk)+ynode,kk-ynode+1:kk)=0;
                end
            end
            
            if kk <N-(ynode*xnode)+1
                CM(kk,kk+(ynode*xnode))=C_i_j_kp1;
                
                CM(kk+(ynode*xnode),kk)=C_i_j_km1;
                
                
            end
            
        end
        
   
        
        %%  %% for continous boundary z
%         T_top_Z=0;
%         T_bot_Z=0;
%         diag_term_Cz_R=-k*2*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));
%         
%         for ii=1:length(BC_bot_Z)
%             CM(BC_bot_Z(ii),BC_bot_Z(ii))=diag_term_Cz_R;  % only for the right side
%             
%         end
%         
%         for ii=1:length(BC_top_Z)
%             CM(BC_top_Z(ii),BC_top_Z(ii))=diag_term_Cz_R;  % only for the right side
%             
%         end
        
        
   
T_top_Z=0;
T_bot_Z=0;
%     diag_term_Cz_R=-1*k* ( 2*(1/delta_x)^2 +2*(1/delta_y)^2+ (1/delta_z)^2 )-(v*(1/delta_x))-transient_term;

diag_term_Cz_R=-k*2*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));


for ii=1:length(BC_bot_Z)
    CM(BC_bot_Z(ii),BC_bot_Z(ii))=diag_term_Cz_R;  % only for the right side
    
    CM(BC_top_Z(ii),BC_top_Z(ii))=diag_term_Cz_R;  % only for the right side
    
end



diag_term_Cx_R=-k*1*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));
diag_term_Cy_R=-k*2*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));

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

diag_term_Cxz=-k*1*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));

for ii=1:length(BC_L_Xleft_Zbot)
    CM(BC_L_Xleft_Zbot(ii),BC_L_Xleft_Zbot(ii))=diag_term_Cxz;
    CM(BC_L_Xleft_Ztop(ii),BC_L_Xleft_Ztop(ii))=diag_term_Cxz;
    CM(BC_L_Xright_Zbot(ii),BC_L_Xright_Zbot(ii))=diag_term_Cxz;
    CM(BC_L_Xright_Ztop(ii),BC_L_Xright_Ztop(ii))=diag_term_Cxz;
    
end


diag_term_Cyz=-k*2*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));

for ii=1:length(BC_L_Ybot_Zbot)
    CM(BC_L_Ybot_Zbot(ii),BC_L_Ybot_Zbot(ii))=diag_term_Cyz;
    CM(BC_L_Ytop_Zbot(ii),BC_L_Ytop_Zbot(ii))=diag_term_Cyz;
    CM(BC_L_Ybot_Ztop(ii),BC_L_Ybot_Ztop(ii))=diag_term_Cyz;
    CM(BC_L_Ytop_Ztop(ii),BC_L_Ytop_Ztop(ii))=diag_term_Cyz;
    
end

diag_term_Cxy=-k*1*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*2* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));

for ii=1:length(BC_L_Xleft_Ybot)
    CM(BC_L_Xleft_Ybot(ii),BC_L_Xleft_Ybot(ii))=diag_term_Cxy;
    CM(BC_L_Xleft_Ytop(ii),BC_L_Xleft_Ytop(ii))=diag_term_Cxy;
    CM(BC_L_Xright_Ybot(ii),BC_L_Xright_Ybot(ii))=diag_term_Cxy;
    CM(BC_L_Xright_Ytop(ii),BC_L_Xright_Ytop(ii))=diag_term_Cxy;
    
end


% points BC


diag_term_Cxyz=-k*1*(1/delta_x)^2 -ky*1*(1/delta_y)^2 -kz*1* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y));


CM(BC_P_Xleft_Ybot_Z_bot,BC_P_Xleft_Ybot_Z_bot)=diag_term_Cxyz;
CM(BC_P_Xleft_Ytop_Z_bot,BC_P_Xleft_Ytop_Z_bot)=diag_term_Cxyz;
CM(BC_P_Xleft_Ybot_Z_top,BC_P_Xleft_Ybot_Z_top)=diag_term_Cxyz;
CM(BC_P_Xleft_Ytop_Z_top,BC_P_Xleft_Ytop_Z_top)=diag_term_Cxyz;


CM(BC_P_Xright_Ybot_Z_bot,BC_P_Xright_Ybot_Z_bot)=diag_term_Cxyz;
CM(BC_P_Xright_Ytop_Z_bot,BC_P_Xright_Ytop_Z_bot)=diag_term_Cxyz;
CM(BC_P_Xright_Ybot_Z_top,BC_P_Xright_Ybot_Z_top)=diag_term_Cxyz;
CM(BC_P_Xright_Ytop_Z_top,BC_P_Xright_Ytop_Z_top)=diag_term_Cxyz;


    
        
        
        %%
        
        
        % convection in the bottom surface
        %to be more efficient
        % Create separate index and value arrays.
        % Call sparse to assemble the index and value arrays.
   
          % convection air upprt and bottom surface
    
    Conv_inside=convection_term*Temp_Conv_inside;
        
    Conv_amb=convection_term*T_amb;
    RHS(1:N_plane)=RHS(1:N_plane)+ Conv_amb;
    % Create separate index and value arrays.
    % Call sparse to assemble the index and value arrays.

        % >> modified on 12 July 2022
for hh=1:N_plane  
        CM(hh,hh)=CM(hh,hh)+convection_term;
%         RHS(hh)=RHS(hh)+ Conv_inside;
    end

    
    for hh=(znode-1)*N_plane+1:N    % for convection of Top_z surface
        CM(hh,hh)=CM(hh,hh)+convection_term;
        RHS(hh)=RHS(hh)+ Conv_inside;
    end
        
        
        % Note: The effect of continoius boundary is not that much in case of
        % high velocity
        

        
        
        %      spy(CM);
        %Now coeffiecient matrix has been constructed
        
        %        RHS(BC_right)=RHS(BC_right)-((k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
        %        RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
        %
        %         RHS(BC_bottom)=RHS(BC_bottom)-((ky*(1/delta_y)^2)+(v_2*(1/delta_y)))*T_bottom;
        %         RHS(BC_top)=RHS(BC_top)-(ky*(1/delta_y)^2)*T_top;
        %
        %         RHS(BC_top_Z)=RHS(BC_top_Z)-(kz*(1/delta_z)^2)*T_top_Z;
        %         RHS(BC_bot_Z)=RHS(BC_bot_Z)-(kz*(1/delta_z)^2)*T_bot_Z;
        
        RHS(BC_right)=RHS(BC_right)-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
        
        % 2*2 lines
        RHS([BC_L_Xright_Zbot BC_L_Xright_Ztop])=RHS([BC_L_Xright_Zbot BC_L_Xright_Ztop])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
        RHS([BC_L_Xright_Ybot, BC_L_Xright_Ytop])=RHS([BC_L_Xright_Ybot, BC_L_Xright_Ytop])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
        
        %          4 corner points
        RHS([BC_P_Xright_Ybot_Z_bot,  BC_P_Xright_Ytop_Z_bot, BC_P_Xright_Ybot_Z_top, BC_P_Xright_Ytop_Z_top])=RHS([BC_P_Xright_Ybot_Z_bot,  BC_P_Xright_Ytop_Z_bot, BC_P_Xright_Ybot_Z_top, BC_P_Xright_Ytop_Z_top])-((0*k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
        
        
        RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
        
        RHS(BC_bottom)=RHS(BC_bottom)-((ky*(1/delta_y)^2)+(v_2*(1/delta_y)))*T_bottom;
        RHS(BC_top)=RHS(BC_top)-(ky*(1/delta_y)^2)*T_top;
        
        RHS(BC_top_Z)=RHS(BC_top_Z)-(kz*(1/delta_z)^2)*T_top_Z;
        RHS(BC_bot_Z)=RHS(BC_bot_Z)-(kz*(1/delta_z)^2)*T_bot_Z;
        
        
        
        
        
        T= RHS'/CM ;
        
        %    [L,U,P] = lu(CM);
        % T= U\(L\(P*RHS));
        
        % [L,U,p] = lu(CM,'vector');
        % T = U\(L\(RHS(p,:)));
        
%         toc
        
        T=full(T) ;
        
        %% for extracting the measured box
        if fileID_Temp_Red_Box
    Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box_Sub,delta_x,delta_y,index_middle_long);
        fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
        end
        
        %% For saving the Nip-point temperature matrix
       T_nip_Sub= T(BC_Left)';
      save('T_nip_Sub.mat','T_nip_Sub');
        %%
        
        
        
        
        
        % heating surface
        T_matrix=reshape(T(N_plane:-1:1),ynode,xnode);
        
%         %% Graphical output
%         

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
        %         Frames(tt)=getframe;
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
        set(h46,'Visible','off');
%         subplot(3,2,vv);
           subplot(2,1,2);
        surf(X_2D,Y_2D,T_matrix_bott);
        shading flat
       
        %         Frames(tt)=getframe;
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
        %%
        
        

        
      
        
        
        
        
        %% Tecplot
        if Tecplot_check
            
        fileID6 = fopen('.\Temp_3D_Sub_Steady.plt','w');
    
        
        tt=1;
        
        fprintf(fileID6,' "TITLE = 3D-Temperature  by Amin zaami" \n');
        fprintf(fileID6,'VARIABLES = "X", "Y", "Z", "Temp" \n');
        fprintf(fileID6,'zone  N=    %d E=   %d DATAPACKING = POINT, ZONETYPE = FEBRICK  \n',N, (xnode-1)*(ynode-1)*(znode-1) );
        fprintf(fileID6,'STRANDID=1, SOLUTIONTIME=  %d \n',tt );
        
        eleNx=xnode;
        eleNy=ynode;
        eleNz=znode;
        
        
        xnode_P=linspace(0,Lx,eleNx);
        ynode_P=linspace(0,Ly,eleNy);
        znode_P=linspace(0,Lz,eleNz);
        counter=0;
        
        
        for kk=1:znode
            for ii=1: xnode
                
                for jj=1: ynode
                    counter=counter+1;
                    
                    fprintf(fileID6,' %f %f  %f %f \r\n', xnode_P(ii), ynode_P(jj), znode_P(kk), T(counter));
                end
                
            end
        end
        
        
        
        %% element creation
        
        shift=(eleNx)*(eleNy);
        num_ele_plane=(eleNx-1)*(eleNy-1);  %number of element in a plane
        
        node_plane=zeros(4,1);
        node_plane_up=zeros(4,1);
        
        
        counter_e=0;
        n_th=0;
        for ii=1: (eleNx-1)*(eleNy-1)*(eleNz-1)    % number of element
            
            if rem(ii,eleNy-1)==1
                counter_e=counter_e+1;
            end
            
            node_plane=[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
            node_plane_up=shift+[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
            
            fprintf(fileID6,' %d %d %d %d %d %d %d %d \r\n', [node_plane node_plane_up]);
            
            if mod(ii,num_ele_plane)==0
                n_th=n_th+1;
                counter_e=n_th*shift-1;
                
            end
            
            counter_e=counter_e+1;
            
        end
        
        
        
        
  
    fclose(fileID6);

        end

% movie(Frames,20)