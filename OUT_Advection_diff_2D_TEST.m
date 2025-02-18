
%% For different velocity
%%The problem of 2D case
% initial temperature withheat flux boundary condition
%Purpose: to investigate effect of velocity in solution resluts
% ******  problem with convective boundary condition  ******
% Challenge  among velocity, convection, heat flux, and diffusion, which all
% depends on the enough time for interaction with the surfaces

clc;
clear;
close all;
h=1000; % W/m^2

%thermal conductivity tape
k=0.72; %e-3;
%density
ro=1560; %e-9;
%sprecific heat capacity
cp=1425; %e0;
%thermal diffusivity
alpha=(k/ro*cp);
%%
Lx=40e-3;
Ly=20e-3;
thickness=1.5e-4;
A=Lx*Ly;

xnode=24;
ynode=36;
delta_x=Lx/xnode;
delta_y=Ly/ynode;

N=ynode*xnode;


% Heat flux is constant here
g_dot=400/thickness;
g_dot=g_dot/(delta_x*delta_y);  % for the lement size
g_dot=g_dot/(xnode*ynode);  % total energy devided per number of nodes for each node

T_BC=sparse(N,1);

T_right=0;
T_left=0;
T_top=0;
T_bottom=0;

BC_Left=1:ynode;
BC_top=ynode:ynode:N;
BC_bottom=1:ynode:N-1;
BC_right=N:-1:N-ynode+1 ;

for kk=1:xnode-1
        index_middle_long(kk)=floor(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=1; % number before the last row
index_nip_point=(ynode*(xnode-row_nip-1)+1):((xnode-row_nip)*(ynode));


for vv=1:6
    
    %     tic
    %     T_old=zeros(1,N);
    CM=sparse(N);
    %         CM=zeros(N);
    
    for kk=1:N
        CM(kk,kk)=1;
     end
    
    %  CM=eye(N);
    %  T=sparse(1,N);
    % T=zeros(1,N);
    % RHS=zeros(N,1);
    RHS=sparse(N,1);
    
    RHS(1:end)=-g_dot;
    
    % The sign of velocity is important because the program uses upwind method in the direction of the velocity the computational points are calculated
    v_original=((vv-1)/10);
    v=ro*cp*v_original;%/ro/delta_x/delta_y;
    T_amb=20; % out of plane temperature
    %%
    pec=(v_original*delta_x)/alpha;
    
    % convection_term=-h;%*A;
    % convection for all nodes
    
    convection_term=-h*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term=convection_term/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3
%     convection_term=convection_term/(A*thickness);% because delta_z is constant
    
   diag_term=-2*k*((1/delta_x^2)+(1/delta_y^2))-(v/delta_x)+convection_term;
    diag_term_C=-2*k*((1/delta_x^2)+(1/delta_y^2));
    
    CM_C=CM;
    
    CM= CM*diag_term;
    CM_C= CM*diag_term_C;  % CM constant
    
    
    
    
    %% for continous boundary x
    diag_term_CB_xR=-1*k*((1/delta_x^2)+(2/delta_y^2))-(v/delta_x)+convection_term;
    T_right=0;
    
    for ii=1:length(BC_right)
        
        CM(BC_right(ii),BC_right(ii))=diag_term_CB_xR;  % for right side
        
    end
        
    %continous in y-direction
    %% for continous boundary
     diag_term_CB_y=-1*k*((2/delta_x^2)+(1/delta_y^2))-(v/delta_x)+convection_term;
    
      for ii=2:xnode-1
    index_y=N-(ynode*(ii-1));
         CM(index_y,index_y)=diag_term_CB_y;  %  ??
         CM(index_y-ynode+1,index_y-ynode+1)=diag_term_CB_y;
    
      end
    
%     for four corners
     indexxy=[N,N-ynode+1; 1 ynode];
     diag_term_CB_xy=-1*k*((1/delta_x^2)+(1/delta_y^2))-(v/delta_x)+convection_term;
     for ii=1:4
         index_xy=indexxy(ii);
         CM(index_xy,index_xy)=diag_term_CB_xy;  %  ??
    
     end
    
    
    %define temp terms
    C_ip1_j=(k/delta_x^2);
    C_im1_j=(k/delta_x^2)+(v/delta_x);
    C_i_jp1=(k/delta_y^2);
    C_i_jm1=(k/delta_y^2);
    
    for kk=1:N
        
        if mod(kk,ynode)~=0
            if kk <N
                CM(kk,kk+1)=C_i_jp1;
                CM(kk+1,kk)=C_i_jm1;
                
                CM_C(kk,kk+1)=C_i_jp1;
                CM_C(kk+1,kk)=C_i_jm1;
            end
            
            %                RHS(kk,1)=-g_dot/delta_x/delta_y;
        end
        
        if kk <N-ynode+1
            CM(kk,kk+ynode)=C_ip1_j;
            CM_C(kk,kk+ynode)=C_ip1_j;
            
            CM(kk+ynode,kk)=C_im1_j;
            
            
        end
        
    end
    
    tic
    
    Advection_term_diag=-(v/delta_x);
    
    for kk=1:N
        CM_C(kk,kk)=CM_C(kk,kk)+Advection_term_diag;
    end
    
    
    % jj=1+ynode:N;  % for the coloumn
    % ii=1:N-ynode;  % for the row
    % index=N*(jj-1)+ii;
    % CM(index)=C_ip1_j;
    %
    %
    % % %      Velocity term
    
    
    
    ii=1+ynode:N;  % for the coloumn
    jj=1:N-ynode;  % for the row
    
    index=N*(jj-1)+ii;
    CM_C(index)=C_im1_j;
    
    
    % for i=1:N-ynode
    %
    %         CM_C(ii(i),jj(i))=C_im1_j;
    %
    % end
    
    
    
    % for ii=1:N-ynode
    % CM_C(index(ii))=C_im1_j;
    % end
    
    
    
    % %
    % %
    % % tt=(ynode*(0:xnode-1));
    % % ii=2:ynode;
    % % jj=1:ynode-1;
    % %
    % % for cc=1:length(tt)
    % %
    % % A=ii+tt(cc);
    % % B=jj+tt(cc);
    % %
    % %
    % % index1=N*(B-1)+A;
    % % index2=N*(A-1)+B;
    % %
    % %
    % %
    % % index3=N*(B(2)-1)+A(2);
    % % index4=N*(A(2)-1)+B(2);
    % %
    % %
    % %  CM(index1)=C_i_jm1;
    % %   CM(index3)=C_i_jm1;
    % %     CM(index2)=C_i_jp1;
    % %   CM(index4)=C_i_jp1;
    % %
    % %
    % % end
 
    
    
    %      spy(CM);
    % % % implementing BC
    
    RHS(BC_Left)=RHS(BC_Left)-((k/delta_x^2)+(v/delta_x))*T_left;
    RHS(BC_bottom)=RHS(BC_bottom)-(k/delta_y^2)*T_bottom;
    RHS(BC_top)=RHS(BC_top)-(k/delta_y^2)*T_top;
    RHS(BC_right)=RHS(BC_right)-(k/delta_x^2)*T_right;
    
    RHS=RHS+convection_term*T_amb;
    
    
    %% ******** Solution ********
    %    T= RHS'/CM ;
    
    %       [L,U,P] = lu(CM);
    % T= U\(L\(P*RHS));
    
    [L,U,p] = lu(CM,'vector');
    T = U\(L\(RHS(p,:)));
    
    toc
    
        T_matrix=reshape(T,ynode,xnode);
     figure(1);
    subplot(3,2,vv);
    
    [X,Y]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
    surf(X,Y,T_matrix);
    
    %       mesh(T_matrix);
    %      shading flat
    shading interp
    view([0 90])
    
    hold on
    axis equal
    title(sprintf('Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    axis([0 Lx 0 Ly]);
        colorbar;
    
    figure(2);
    subplot(3,2,vv);
    hold on;
    % title('Tmepereture along middle of the tape');
    title(sprintf('Tmepereture along long, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    ylabel('Temp');
    plot(T(index_middle_long));
    
    figure(3);
    subplot(3,2,vv);
    hold on;
    ylabel('Temp');
    % title('Tmepereture along width nip-point of the tape');
    title(sprintf('Tmepereture along width nip-point, Velocity= %f , Peclet: %f, h=%d ',v_original,pec,h));
    plot(T(index_nip_point));
    
    
end








