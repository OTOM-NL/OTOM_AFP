
%% For different velocity - transient




function T_matrix=FDM_2D_transient_combined_Tape(Lx,Ly,xnode,ynode,E_points,v_original,materials,Temp_Incoming,h,thickness,...
               Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box,Graphic_Profile,Graphic_contour,Tecplot_check)


% Video = VideoWriter('animation-2D.avi','Uncompressed AVI');
% video.FrameRate = 60;
% % 
% open(Video);


Time=10;% from user
inc=50;% from user
delta_t=Time/inc;


k=  materials(1);
    %density
    ro= materials(2);
    %sprecific heat capacity
    cp= materials(3);
    
    
     
    if length(materials)>3
        ky=materials(4);
        kz=materials(5);
    else
        ky=k;
        kz=k;
        
    end



    %thermal diffusivity
alpha=(k/ro*cp);
A=Lx*Ly;

delta_x=Lx/xnode;
delta_y=Ly/ynode;

N=ynode*xnode;

% Heat flux is constant here
g_dot=1/thickness;
g_dot=g_dot/(delta_x*delta_y);  % for the lement size
g_dot=g_dot/(xnode*ynode);  % total energy devided per number of nodes for each node


%   T_amb=25;  % Temperature for convection from surface, not edge!!
T_right=0;
T_left=Temp_Incoming;  % incoming velocity
T_top=0;
T_bottom=0;

 T_amb=20; % out of plane temperature


BC_Left=1:ynode;
BC_top=ynode:ynode:N;
BC_bottom=1:ynode:N-1;
BC_right=N:-1:N-ynode+1 ;

BC_corners=[1 ynode N-ynode+1 N ];

for kk=1:xnode-1
        index_middle_long(kk)=floor(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=0; % number before the last row
index_nip_point=(ynode*(xnode-row_nip-1)+1):((xnode-row_nip)*(ynode));



transient_term=ro*cp/delta_t;
Temp=T_amb*ones(N,1);


% %      figure(1);
%      pos=get(gcf,'Position');

     
     
     
figure(62);  % transient data output representation

  hold on;
    p1=[];
    p2=[];
    p3=[];
    p4=[];

% for vv=2:2
    
    for tt=1:inc
    
    T_old=Temp;
    
        tic
    %     T_old=zeros(1,N);
    CM=sparse(N);
    %         CM=zeros(N);
    
    for kk=1:N
        CM(kk,kk)=1;
     end
    
    %  CM=eye(N);
    %  Temp=sparse(1,N);
    % Temp=zeros(1,N);
    % RHS=zeros(N,1);
 
    % The sign of velocity is important because the program uses upwind method in the direction of the velocity the computational points are calculated

    
      RHS=sparse(N,1);
        RHS(1:end)=-E_points*g_dot *N; % /N because of intensity input power come as a matrix
     
 v=ro*cp*v_original;
  pec=(v_original*delta_x)/alpha;
    
    
   
    %%
    pec=(v_original*delta_x)/alpha;
    
    % convection_term=-h;%*A;
    % convection for all nodes
    
%     convection_term=-h*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
%         convection_term=convection_term/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3
% %     convection_term=convection_term/(A*thickness);% because delta_z is constant
    
    h1=h(1);
     h2=h(2);

  convection_term=-h1*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term=convection_term/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3
        
       convection_term_Roller=-h2*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term_Roller=convection_term_Roller/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3

           

  diag_term=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)+convection_term-transient_term;
 
 diag_term_Roller_contact=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)+convection_term_Roller-transient_term;

index_Roller_con=1+(Row_number_Tape_Roller_tangent*ynode);


    
    CM= CM*diag_term;

for ii=index_Roller_con:N

    CM(ii,ii)=diag_term_Roller_contact;
    
end 

    %continous in y-direction
    %% for continous boundary
         diag_term_CB_y=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2)-(v/delta_x)+convection_term-transient_term;
        
    
        for ii=1:length(BC_top)
            
            CM(BC_top(ii),BC_top(ii))=diag_term_CB_y;
            CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y;
            
        end
        
    

    
    
    %% for continous boundary x
      diag_term_CB_xR_nip=-1*k*(1/delta_x^2)-1*ky*(2/delta_y^2) -(v/delta_x)+convection_term_Roller-transient_term;
    diag_term_CB_xR_left=-1*k*(1/delta_x^2)-1*ky*(2/delta_y^2) -(v/delta_x)+convection_term-transient_term;
        

    
    for ii=1:length(BC_right)
        
        CM(BC_right(ii),BC_right(ii))=diag_term_CB_xR_nip;  % for right side
          CM(BC_Left(ii),BC_Left(ii))=diag_term_CB_xR_left;  % for right side
            
    end
        
    
    
        diag_term_CB_Corners=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term-transient_term;
        
        for ii=1:length(BC_corners)
            CM(BC_corners(ii),BC_corners(ii))=diag_term_CB_Corners;
            
        end
        
    
    
    %define temp terms

    
       C_ip1_j=(k/delta_x^2);
    C_im1_j=(k/delta_x^2)+(v/delta_x);
    C_i_jp1=(ky/delta_y^2);
    C_i_jm1=(ky/delta_y^2);
    
    for kk=1:N
        
        if mod(kk,ynode)~=0
            if kk <N
                CM(kk,kk+1)=C_i_jp1;
                CM(kk+1,kk)=C_i_jm1;
                
%                 CM_C(kk,kk+1)=C_i_jp1;
%                 CM_C(kk+1,kk)=C_i_jm1;
            end
            
            %                RHS(kk,1)=-g_dot/delta_x/delta_y;
        end
        
        if kk <N-ynode+1
            CM(kk,kk+ynode)=C_ip1_j;
%             CM_C(kk,kk+ynode)=C_ip1_j;
            
            CM(kk+ynode,kk)=C_im1_j;
            
            
        end
        
    end
    
  
    
    
    %      spy(CM);
    % % % implementing BC
    
%     RHS(BC_Left)=RHS(BC_Left)-((k/delta_x^2)+(v/delta_x))*T_left;
%     RHS(BC_bottom)=RHS(BC_bottom)-(k/delta_y^2)*T_bottom;
%     RHS(BC_top)=RHS(BC_top)-(k/delta_y^2)*T_top;
%     RHS(BC_right)=RHS(BC_right)-(k/delta_x^2)*T_right;
    
      RHS(BC_Left)=RHS(BC_Left)-((k/delta_x^2)+(v/delta_x))*T_left;
    RHS(BC_bottom)=RHS(BC_bottom)-(ky/delta_y^2)*T_bottom;
    RHS(BC_top)=RHS(BC_top)-(ky/delta_y^2)*T_top;
    RHS(BC_right)=RHS(BC_right)-(k/delta_x^2)*T_right;
    
    RHS(1:index_Roller_con-1)=RHS(1:index_Roller_con-1)+convection_term*T_amb- (transient_term*T_old(1:index_Roller_con-1));
    
    
     RHS(index_Roller_con:N)=RHS(index_Roller_con:N)+convection_term_Roller*T_amb- (transient_term*T_old(index_Roller_con:N));
    
    
    
    
    
    
    
%     RHS=RHS+convection_term*T_amb - (transient_term*T_old);
    
    
    %% ******** Solution ********
    %    Temp= RHS'/CM ;
    
    %       [L,U,P] = lu(CM);
    % Temp= U\(L\(P*RHS));
    
    [L,U,p] = lu(CM,'vector');
    Temp = U\(L\(RHS(p,:)));
    
 
    toc
    
    T=full(Temp);
    
    if fileID_Temp_Red_Box
    Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box_Tape,delta_x,delta_y,index_middle_long);     
          fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
    end
    
    
        %% For saving the Nip-point temperature matrix
%        T_nip_Tape= T(BC_right);
%       save('T_nip_Tape.mat','T_nip_Tape');
    
    
       %% Matlab Tecplot creator animation
     if Tecplot_check
       Tecplot_Matlab_creator_2D (tt,N,xnode,ynode,v_original,Lx,Ly,T);
     end
    %        
    
       %%
  
    
    delete([p1  p3 ]);
    % along long and width
   title(sprintf('Time = %f',(tt/inc)*Time));
    subplot(1,2,1);
    %     hold on;
    
    % title('Tmepereture along middle of the tape');
%     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    ylabel('Temp');
    p1=plot(linspace(0,Lx,length(index_middle_long)),T(index_middle_long),'b');
    legend ('along Length');
%     subplot(2,2,2)
%     p2=plot(linspace(0,Lx,length(index_middle_long_Bott)),T(index_middle_long_Bott),'r');
%     legend ('Bottom');
%     drawnow;
    
    
%     figure(31);
    subplot(1,2,2);
    %     hold on;
    
    ylabel('Temp');
    % title('Tmepereture along width nip-point of the tape');
%     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    p3=plot(linspace(0,Ly,length(index_nip_point)),T(index_nip_point),'b');
    legend ('along nip-point');
%     subplot(2,2,4);
%     p4=plot(linspace(0,Ly,length(index_nip_point_Bott)),T(index_nip_point_Bott),'r');
%     legend ('Bottom');
    drawnow;
    
    
    
    
    
    
    
    
%         T_matrix=reshape(Temp,ynode,xnode);
 

%     subplot(3,2,vv);
%     
%     [X,Y]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
%     surf(X,Y,T_matrix);
%     
%     %       mesh(T_matrix);
%     %      shading flat
%     shading interp
%     view([0 90]);
%     
% %     hold on
%     axis equal;
%     title(sprintf('Velocity= %f , Peclet: %f, h=%d, Time=%f',v_original,pec,h,tt*delta_t));
%     axis([0 Lx 0 Ly]);
%        colorbar;
% %        caxis([0 500]);
%     
%    pause(0.005);
% %     Mov=getframe;
%  Mov(tt)=getframe(gcf);

 
    
    end

 % %      figure(1);
%      pos=get(gcf,'Position');


    
if Graphic_Profile
    
   h21=figure(21);
   set(h21,'Visible','off');
   
    subplot(2,2,3);
    hold on;
    % title('Tmepereture along middle of the tape');
    title(sprintf('Tape-Temperature along length, V= %f , Pec: %f, h=%d ',v_original,pec,h));
       title(sprintf('Tape-Temperature along length, V=%f , Pec=%f, h1,2=%d %d ',v_original,pec,h1,h2));
    
      ylabel('Temperature ^{\circ}C');
     xlabel('Length (m)');
    
        H_long=linspace(0,Lx,length(T(index_middle_long))); % H_long > horizental 
        plot(H_long,T(index_middle_long));
        
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
    title(sprintf('Tape-Temperature along width nip-point, V=%f , Pec=%f, h1,2=%d %d ',v_original,pec,h1,h2));
    
    T_nip_point=T(index_nip_point);
     H_Nip=linspace(0,Ly,length(T_nip_point)); % H_Nip > horizental 
    plot(H_Nip,T_nip_point);
    Nip_point_Temp_Tape=T_nip_point(floor(length(index_nip_point)/2));
    
    assignin('base','Nip_point_Temp_Tape',Nip_point_Temp_Tape);
    
%     plot(T(index_nip_point));
    
end
  
%% outPut 3D
if Graphic_contour
h100=figure(100);
  set(h100,'Visible','off');
  
  
    T_matrix2=reshape(T,ynode,xnode);
   hold on;
subplot(2,1,2)

[X_2D,Y_2D]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
    
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

    len_T=length(T);
    T_disp=zeros(len_T,1);
    
    
    counter=0;
    

  
% N=xnode*ynode;
    for ii=1:xnode
        for jj=ynode:-1:1
            counter=counter+1;
            index=N-(jj)+1-((ii-1)*ynode);
            T_disp(index)=T(counter);
            
        end
    end
    

    
    T_matrix=reshape(T_disp,ynode,xnode);


%  figure(11);
%     subplot(3,2,vv);
%     hold on;
%     % title('Tmepereture along middle of the tape');
%     title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
%     ylabel('Temp');
%     plot(linspace(0,Lx,length(index_middle_long)),Temp(index_middle_long));
%     
%     figure(12);
%     subplot(3,2,vv);
%     hold on;
%     ylabel('Temp');
%     % title('Tmepereture along width nip-point of the tape');
%     title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
%     plot(linspace(0,Ly,length(index_nip_point)),Temp(index_nip_point));


% figure;
% axis off;
% movie(Mov);
% 
% writeVideo(Video,Mov);
% close(Video);
