%%The problem of 3D advection-diffusion case
%%%%%%%%%% ******  problem with convective boundary condition  ******


%% Note >>> convection has a problem with z-node by changing it

function T_matrix=FDM_3D_transient_combined(Lx,Ly,Lz,xnode,ynode,znode,E_points,v_original,materials,Temp_Incoming,h,T_amb,...
       Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box,Tecplot_check,Graphic_profile,Grpahic_contour)


% load initial value
% load T_initial;

% if video_ID % from user
% Video = VideoWriter('animation-3D.avi','Uncompressed AVI');
% video.FrameRate = 60;
%
% open(Video);

 fileID6 = fopen('Temp_3D-FDM_V.plt','w');


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



Time=10; % from user
inc=20;% from user
delta_t=Time/inc;





delta_x=Lx/(xnode-0); % it was Lx/(xnode-1);
delta_y=Ly/(ynode-0);
delta_z=Lz/(znode-0);

N=ynode*xnode*znode;
N_plane=ynode*xnode;

advection_diffusion_3D_nodes(Lx,Ly,Lz,xnode,ynode,znode)

%% How the heat flux is implemented should be investigated !?




g_dot=1/((xnode-0)*(ynode-0));
g_dot= g_dot/(delta_x*delta_y*delta_z);


% T_BC=sparse(N,1);

fileID1 = fopen('nodes.txt');

node = textscan(fileID1,' %*d %f %f %f ','Delimiter',',') ;
Nodes=cell2mat(node);


% T_left=[ T_nip_Tape; T_nip_Sub];
T_right=20;
T_left=Temp_Incoming;
T_top=20;
T_bottom=20;

BC_Left=[];
BC_top=[];
BC_bottom=[];
BC_right=[];

for ii=0:znode-1
    BC_Left=[BC_Left (1:ynode) + (N_plane *ii)];
    BC_top=[BC_top (ynode:ynode:N_plane)+(N_plane *ii)];
    BC_bottom=[BC_bottom (1:ynode:N_plane-1)+(N_plane *ii)];
    BC_right=[BC_right (N_plane:-1:N_plane-ynode+1)+(N_plane *ii) ];
    
end

BC_bot_Z=1:N_plane;
BC_top_Z=N_plane*(znode-1)+1:N;



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

%%


transient_term=ro*cp/delta_t;

% Temp=zeros(1,N);

Temp=T_amb*ones(N,1);

% CM_orig=sparse(1:N,1:N,ones(1,N),N,N);

% step=1;


    v0=ro*cp*v_original;%/(delta_z^1);25*25*
    
    v=v0*abs(cosd(Rel_v_theta));
    v_2=v0*abs(sind(Rel_v_theta));


pec=(v_original*delta_x)/alpha;



%%

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


%%


    convection_term=-h*A/(xnode*ynode);%*((delta_x*delta_y*delta_z)^2);
    
    convection_term=convection_term/(delta_x*delta_y*delta_z);% because delta_z is constant



figure(61);  % transient data output representation

  hold on;
    p1=[];
    p2=[];
    p3=[];
    p4=[];

for tt=1:inc
    
      tic;
      
    T_old=Temp';
    
  

    RHS=zeros(N,1);
    %         RHS=sparse(N,1);
%      RHS(1:N_plane)=-g_dot;

     
%         RHS=sparse(N,1);
        RHS(1:N_plane)=-E_points*g_dot*N_plane;
    
    
    
%     diag_term=-2*k* ( (1/delta_x^2) +(1/delta_y^2)+ (1/delta_z^2) )-(v*(1/delta_x))-transient_term;
    
    diag_term=-2*k* (1/delta_x^2) - 2*ky*(1/delta_y^2)- 2*kz* (1/delta_z^2) -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
    
%     CM= CM*diag_term;
%     C_ip1_j_k=(k*(1/delta_x)^2);
%     C_im1_j_k=(k*(1/delta_x)^2)+(v*(1/delta_x));
%     C_i_jp1_k=(k*(1/delta_y)^2);
%     C_i_jm1_k=(k*(1/delta_y)^2);
%     
%     C_i_j_kp1=(k*(1/delta_z)^2);
%     C_i_j_km1=(k*(1/delta_z)^2);
    
    
    
    
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
    T_top_Z=0;
    T_bot_Z=0;
%     diag_term_Cz_R=-1*k* ( 2*(1/delta_x)^2 +2*(1/delta_y)^2+ (1/delta_z)^2 )-(v*(1/delta_x))-transient_term;
  
        diag_term_Cz_R=-k*2*(1/delta_x)^2 -ky*2*(1/delta_y)^2 -kz* (1/delta_z)^2 -(v*(1/delta_x))-(v_2*(1/delta_y))-transient_term;
        
    
    for ii=1:length(BC_bot_Z)
        CM(BC_bot_Z(ii),BC_bot_Z(ii))=diag_term_Cz_R;  % only for the right side
        
%     end
%     
%     for ii=1:length(BC_top_Z)
        CM(BC_top_Z(ii),BC_top_Z(ii))=diag_term_Cz_R;  % only for the right side
        
    end
    %
    
    %%
    
    
    % convection in the bottom surface
    %to be more efficient
    % Create separate index and value arrays.
    % Call sparse to assemble the index and value arrays.
    
    for hh=(znode-1)*N_plane+1:N    % for convection of Top_z surface
        CM(hh,hh)=CM(hh,hh)+convection_term;
        RHS(hh)=RHS(hh)+ (convection_term*T_amb);
    end
    
    
    for ii=1:length(CM_ind_0)
        CM(CM_ind_0(ii,1),CM_ind_0(ii,2))=0;
    end
    
    
    % Note: The effect of continoius boundary is not that much in case of
    % high velocity
    
    
%     C_ip1_j_k=(k*(1/delta_x)^2);
%     C_im1_j_k=(k*(1/delta_x)^2)+(v*(1/delta_x));
%     C_i_jp1_k=(k*(1/delta_y)^2);
%     C_i_jm1_k=(k*(1/delta_y)^2);
%     
%     C_i_j_kp1=(k*(1/delta_z)^2);
%     C_i_j_km1=(k*(1/delta_z)^2);
    
%     CM_ind_C_i_jp1_k=0;
%     CM_ind_C_i_jm1_k=0;
%     
%     for kk=1:N
%         
%         
%         if mod(kk,ynode)~=0
%             if kk <N
%                 CM(kk,kk+1)=C_i_jp1_k;
%                 CM(kk+1,kk)=C_i_jm1_k;
%                 
% 
% %                 
%             end
%             
%         end
%         
%         if kk <N-ynode+1
%             if mod(kk,N_plane)~=0
%                 
%                 CM(kk,kk+ynode)=C_ip1_j_k;
%                 CM(kk+ynode,kk)=C_im1_j_k;
% 
%                 
%             else
%                 CM(kk-ynode+1:kk,(kk-ynode+1:kk)+ynode)=0;
%                 CM((kk-ynode+1:kk)+ynode,kk-ynode+1:kk)=0;
%                 
% 
%                 
%             end
%         end
%         
%         if kk <N-(ynode*xnode)+1
%             CM(kk,kk+(ynode*xnode))=C_i_j_kp1;
%             CM(kk+(ynode*xnode),kk)=C_i_j_km1;
%             
%             
% 
%             
%         end
%         
%     end
    
%     indexes


    
  

    
    %      spy(CM);
    %Now coeffiecient matrix has been constructed
    
    
    
    %         RHS(BC_Left)=RHS(BC_Left)-((k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
    %                 RHS(BC_right)=RHS(BC_right)-(k*(1/delta_x)^2)*T_right;
    
    % T_left should be matrix of the incloming conditions Temperature T(BC_Left)
    
%     RHS(BC_right)=RHS(BC_right)-((k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
%     RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
%     
%     
%     
%     RHS(BC_bottom)=RHS(BC_bottom)-(k*(1/delta_y)^2)*T_bottom;
%     RHS(BC_top)=RHS(BC_top)-(k*(1/delta_y)^2)*T_top;
%     
%     RHS(BC_top_Z)=RHS(BC_top_Z)-(k*(1/delta_z)^2)*T_top_Z;
%     RHS(BC_bot_Z)=RHS(BC_bot_Z)-(k*(1/delta_z)^2)*T_bot_Z;
    
    
    
         RHS(BC_right)=RHS(BC_right)-((k*(1/delta_x)^2)+(v*(1/delta_x)))*T_left;
       RHS(BC_Left)=RHS(BC_Left)-(k*(1/delta_x)^2)*T_right;
        
        RHS(BC_bottom)=RHS(BC_bottom)-((ky*(1/delta_y)^2)+(v_2*(1/delta_y)))*T_bottom;
        RHS(BC_top)=RHS(BC_top)-(ky*(1/delta_y)^2)*T_top;

        RHS(BC_top_Z)=RHS(BC_top_Z)-(kz*(1/delta_z)^2)*T_top_Z;
        RHS(BC_bot_Z)=RHS(BC_bot_Z)-(kz*(1/delta_z)^2)*T_bot_Z;
    
    
    
    
    RHS=RHS - (transient_term*T_old);
    
    Temp= RHS'/CM ;
    
    %    [L,U,P] = lu(CM);
    % T= U\(L\(P*RHS));
    
    % [L,U,p] = lu(CM,'vector');
    % T = U\(L\(RHS(p,:)));
    
   
    
    T=full(Temp) ;
     toc;
    
        %% for extracting the measured box
        if fileID_Temp_Red_Box
    Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box_Sub,delta_x,delta_y,index_middle_long);
        fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
        end
     
        
                
    
    
     
    %% Matlab Tecplot creator animation
     if Tecplot_check
       Tecplot_Matlab_creator (tt,N,xnode,ynode,znode,v_original,Lx,Ly,Lz,T,fileID6);
     end
    %         %%
    %
    %
    %
    %         figure(3);
    %
    %         X=zeros(ynode*xnode,znode);
    %         Y=zeros(ynode*xnode,znode);
    %         Z=zeros(ynode*xnode,znode);
    %         T_show=reshape(T,ynode*xnode,znode);
    %
    %         for kk=1:ynode*xnode
    %             index=((1:znode)*xnode*ynode)-kk+1;
    %             X(kk,1:end)=Nodes(index,1)';
    %             Y(kk,1:end)=Nodes(index,2)';
    %             Z(kk,1:end)=Nodes(index,3)';
    %             %      Intensity(kk,1:end)=points(index,4)';
    %         end
    %
    %
    % %         subplot(3,2,vv);
    %         f_obj=mesh(X,Y,Z,T_show,'MarkerFaceColor','auto',...
    %             'Marker','.',...
    %             'LineStyle','none','FaceColor','none');
    %
    %         axis equal;
    %         shading flat;
    %         % shading interp
    %         %      view([0 90])
    %
    % %         hold on
    %
    %         title(sprintf('Velocity= %f , Peclet number: %f, Time=%f ',v_original,pec,tt*delta_t));
    %         %      axis([0 Lx 0 Ly]);
    %         % view([-83.5 -52.4]);
    %         view([42.5 -34.8000000000001]);
    %         colorbar;
    %         % disp(vv);
    %
    % %            axis equal;
    %             Mov(tt)=getframe(gcf);
    %            pause(0.01);
    
    
    
    %%
  
    
    delete([p1 p2 p3 p4]);
    % along long and width
   title(sprintf('Time = %f',(tt/inc)*Time));
    subplot(2,2,1);
    %     hold on;
    
    % title('Tmepereture along middle of the tape');
%     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    ylabel('Temp');
    p1=plot(linspace(0,Lx,length(index_middle_long)),T(index_middle_long),'b');
    legend ('Upper');
    subplot(2,2,2)
    p2=plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
    legend ('Bottom');
%     drawnow;
    
    
%     figure(31);
    subplot(2,2,3);
    %     hold on;
    
    ylabel('Temp');
    % title('Tmepereture along width nip-point of the tape');
%     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    p3=plot(linspace(0,Ly,length(index_nip_point)),T(index_nip_point),'b');
    legend ('Upper');
    subplot(2,2,4);
    p4=plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
    legend ('Bottom');
    drawnow;
    
    
    
end

          %% For saving the Nip-point temperature matrix
       T_nip_Sub= T(BC_Left)';
      save('T_nip_Sub.mat','T_nip_Sub');
  

        % heating surface
        T_matrix=reshape(T(N_plane:-1:1),ynode,xnode);




if Grpahic_contour

        T_matrix_bott=reshape(T(N:-1:N-N_plane+1),ynode,xnode);
        
        
        [X_2D,Y_2D]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
        h46=figure(46);
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
        































% if video_ID
%         writeVideo(Video,Mov);
%         close(Video);
% end

%%
% T_matrix=reshape(T(N_plane:-1:1),ynode,xnode);
% T_matrix_bott=reshape(T(N:-1:N-N_plane+1),ynode,xnode);
% 
% 
% [X_2D,Y_2D]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
% figure(4);
% % subplot(3,2,vv);
% %    subplot(2,2,vv);
% surf(X_2D,Y_2D,T_matrix);
% shading flat
% view([0 90])
% %         Frames(tt)=getframe;
% %       hold on
% %    axis equal
% title(sprintf('Upper-S, velocity= %f , Peclet number: %f ',v_original,pec,tt));
% axis([0 Lx 0 Ly]);
% colorbar;

%% for Bottom surface

% figure(5);
% % subplot(3,2,vv);
% %    subplot(2,2,vv);
% surf(X_2D,Y_2D,T_matrix_bott);
% shading flat
% view([0 90])
% %         Frames(tt)=getframe;
% %       hold on
% %    axis equal
% title(sprintf('Bottom-S, velocity= %f , Peclet number: %f ',v_original,pec,tt));
% axis([0 Lx 0 Ly]);
% colorbar;

%         % along long and width
%           figure(31);
%     subplot(3,2,vv);
%     hold on;
%     % title('Tmepereture along middle of the tape');
%     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
%     ylabel('Temp');
%     plot(linspace(0,Lx,length(index_middle_long)),T(index_middle_long),'b');
%      plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
%      legend ('Upper','Bottom');
%
%
%
%     figure(32);
%     subplot(3,2,vv);
%     hold on;
%     ylabel('Temp');
%     % title('Tmepereture along width nip-point of the tape');
%     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
%     plot(linspace(0,Ly,length(index_nip_point)),T(index_nip_point),'b');
%       plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
%             legend ('Upper','Bottom');



%%









%% Tecplot
%
%
%             fileID6 = fopen(sprintf('Temp_3D-FDM_V_0to%f.plt',v_original),'w');
%
%     %      fileID6 = fopen(sprintf('Temp_3D-FDM_V_0to%f.plt',v_original),'a');
%
%
%
%
%
%         fprintf(fileID6,' "TITLE = 3D-Temperature  by Amin zaami" \n');
%         fprintf(fileID6,'VARIABLES = "X", "Y", "Z", "Temp" \n');
%         fprintf(fileID6,'zone  N=    %d E=   %d DATAPACKING = POINT, ZONETYPE = FEBRICK  \n',N, (xnode-1)*(ynode-1)*(znode-1) );
%         fprintf(fileID6,'STRANDID=1, SOLUTIONTIME=  %d \n',tt );
%
%         eleNx=xnode;
%         eleNy=ynode;
%         eleNz=znode;
%
%
%         xnode_P=linspace(0,Lx,eleNx);
%         ynode_P=linspace(0,Ly,eleNy);
%         znode_P=linspace(0,Lz,eleNz);
%         counter=0;
%
%
%         for kk=1:znode
%             for ii=1: xnode
%
%                 for jj=1: ynode
%                     counter=counter+1;
%
%                     fprintf(fileID6,' %f %f  %f %f \r\n', xnode_P(ii), ynode_P(jj), znode_P(kk), T(counter));
%                 end
%
%             end
%         end
%
%
%
%         %% element creation
%
%         shift=(eleNx)*(eleNy);
%         num_ele_plane=(eleNx-1)*(eleNy-1);  %number of element in a plane
%
%         node_plane=zeros(4,1);
%         node_plane_up=zeros(4,1);
%
%
%         counter_e=0;
%         n_th=0;
%         for ii=1: (eleNx-1)*(eleNy-1)*(eleNz-1)    % number of element
%
%             if rem(ii,eleNy-1)==1
%                 counter_e=counter_e+1;
%             end
%
%             node_plane=[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
%             node_plane_up=shift+[counter_e, counter_e+1,eleNy+counter_e+1,eleNy+counter_e];
%
%             fprintf(fileID6,' %d %d %d %d %d %d %d %d \r\n', [node_plane node_plane_up]);
%
%             if mod(ii,num_ele_plane)==0
%                 n_th=n_th+1;
%                 counter_e=n_th*shift-1;
%
%             end
%
%             counter_e=counter_e+1;
%
%         end
%
%
%
%
%
%     fclose(fileID6);





%% Use BC of

% Save
% T(BC_Left)  is the end of length...
%     should be used as the initial length as T(BC_right) which also includes Tape conditions


% movie(Frames,20)