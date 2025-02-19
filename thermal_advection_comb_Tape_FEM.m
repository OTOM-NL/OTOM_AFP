



function T_matrix=thermal_advection_comb_Tape_FEM(Lx,Ly,xnode,ynode,E_points,v_original,materials,Temp_Right_T_Roller,h,thickness,...
               Row_number_Tape_Roller_tangent,Measure_Box_Tape,fileID_Temp_Red_Box,Graphic_Profile,Graphic_contour,Tecplot_check_T)


%% initialization 
% $ relation between optical and thermal nodes > 
% fid13 = fopen('.\Supp_files\edge_factors.txt');

fid13 = fopen(strcat(cd,'/Supp_files/edge_factors.txt'));

out = textscan(fid13,'%s','delimiter',',');
    edge_factor_T=str2num(out{1,1}{2});
      fclose(fid13);
      
% edge_factor_T=1.85 ;




% No increase in T_Roller
% h_conv_T = h

% tic ;

fileID_Tecplot_Tape = fopen('.\Temp_2D_Tape_Steady.plt','w');
            

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


% T_BC=sparse(N,1);



  T_amb=20;  % Temperature for convection from surface, not edge!!
T_right=0;
T_left=Temp_Right_T_Roller(1);  % incoming velocity
T_top=0;
T_bottom=0;
        
BC_Left=1:ynode;
BC_top=ynode:ynode:N;
BC_bottom=1:ynode:N-1;
BC_right=N:-1:N-ynode+1 ;

% BC_corners=[1 ynode N-ynode+1 N ];
indexxy=[N,N-ynode+1, 1 ynode];


for kk=1:xnode-1
        index_middle_long(kk)=floor(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=1; % number before the last row
index_nip_point=(ynode*(xnode-row_nip-1)+1):((xnode-row_nip)*(ynode));




% Heat flux is constant here
g_dot=1/thickness;
g_dot=g_dot/(delta_x*delta_y);  % for the lement size
g_dot=g_dot/(xnode*ynode);  % total energy devided per number of nodes for each node



g_dot=ones(N,1)*g_dot;
g_dot([BC_bottom BC_top])=g_dot([BC_bottom BC_top])*edge_factor_T;
g_dot([BC_right BC_Left])=g_dot([BC_right BC_Left])*edge_factor_T;





 CM=sparse(N);
    %         CM=zeros(N);
    
    for kk=1:N
        CM(kk,kk)=1;
     end
  RHS=sparse(N,1);
        RHS(1:end)=-E_points.*g_dot *N; % /N because of intensity input power come as a matrix
     
 v=ro*cp*v_original;
  pec=(v_original*delta_x)/alpha;
  
  
  % h(1) for the flat part, and h(2) for the tangent part which have contact
  % with Roller
  
     h1=h(1);
     h2=h(2);

  convection_term=-h1*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term=convection_term/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3

       convection_term_Roller=-h2*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term_Roller=convection_term_Roller/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3

           
        
        
        
%  diag_term=-2*k*((1/delta_x^2)+(1/delta_y^2))-(v/delta_x)+convection_term;
 
  diag_term=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)+convection_term;
 
 
diag_term_Roller_contact=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)+convection_term_Roller;

index_Roller_con=1+(Row_number_Tape_Roller_tangent*ynode);



        CM= CM*diag_term;





    %continous in y-direction
    %% for continous boundary
%          diag_term_CB_y=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2)-(v/delta_x)+convection_term;
%         
%     
%         for ii=1:length(BC_top)
%             
%             CM(BC_top(ii),BC_top(ii))=diag_term_CB_y;
%             CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y;
%             
%         end
%         
    
    %% for continous boundary x
%       diag_term_CB_xR_nip=-1*k*(1/delta_x^2)-1*ky*(2/delta_y^2) -(v/delta_x)+convection_term_Roller;
%     diag_term_CB_xR_left=-1*k*(1/delta_x^2)-1*ky*(2/delta_y^2) -(v/delta_x)+convection_term;
%         
% 
%     
%     for ii=1:length(BC_right)
%         
%         CM(BC_right(ii),BC_right(ii))=diag_term_CB_xR_nip;  % for right side
%           CM(BC_Left(ii),BC_Left(ii))=diag_term_CB_xR_left;  % for right side
%             
%     end
%         
%     
%     
%         diag_term_CB_Corners=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term;
%         
%         for ii=1:length(BC_corners)
%             CM(BC_corners(ii),BC_corners(ii))=diag_term_CB_Corners;
%             
%         end
    
    
    
    
    

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
       
        end
        
        if kk <N-ynode+1
            CM(kk,kk+ynode)=C_ip1_j;
%             CM_C(kk,kk+ynode)=C_ip1_j;
            
            CM(kk+ynode,kk)=C_im1_j;
            
            
        end
        
    end
    
    
    
    for ii=index_Roller_con:N

    CM(ii,ii)=diag_term_Roller_contact;
    
    end
    
    
    
    

        
    %% for continous boundary
diag_term_CB_y=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term;
diag_term_CB_y_after_roller=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term_Roller;

for ii=1:length(BC_top)
    if BC_top(ii) <index_Roller_con
        CM(BC_top(ii),BC_top(ii))=diag_term_CB_y;
        CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y;
        
    else
        CM(BC_top(ii),BC_top(ii))=diag_term_CB_y_after_roller;
        CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y_after_roller;
        
        
    end
    
end

diag_term_CB_xR_right=-1*k*(1/delta_x^2)-2*ky*(1/delta_y^2) -(v/delta_x)+convection_term_Roller;
diag_term_CB_xR_left=-1*k*(1/delta_x^2)-2*ky*(1/delta_y^2) -(v/delta_x)+convection_term;

%     T_right=0;

for ii=1:length(BC_right)
    
    CM(BC_right(ii),BC_right(ii))=diag_term_CB_xR_right;  % for right side
    CM(BC_Left(ii),BC_Left(ii))=diag_term_CB_xR_left;  % for right side
    
end

%     for four corners
diag_term_CB_xy=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term;
diag_term_CB_xy_after_roller=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term_Roller;
%     for four corners
for ii=1:2
    %          index_xy=indexxy(ii);
    CM(indexxy(ii),indexxy(ii))=diag_term_CB_xy_after_roller;  %  ??
end
for ii=3:4
    %          index_xy=indexxy(ii);
    CM(indexxy(ii),indexxy(ii))=diag_term_CB_xy;  %  ??
end


    
    
    
    
    
    
    
    % Should be changed
    T_Roller=Temp_Right_T_Roller(2);
    
    
%      Advection_term_diag=-(v/delta_x);
%     
%     for kk=1:N
%         CM_C(kk,kk)=CM_C(kk,kk)+Advection_term_diag;
%     end
%     
    
    
    RHS(BC_Left)=RHS(BC_Left)-(0*(k/delta_x^2)+(v/delta_x))*T_left;
%     RHS(BC_bottom)=RHS(BC_bottom)-(ky/delta_y^2)*T_bottom;
%     RHS(BC_top)=RHS(BC_top)-(ky/delta_y^2)*T_top;
%     RHS(BC_right)=RHS(BC_right)-(k/delta_x^2)*T_right;
    
    RHS(1:index_Roller_con-1)=RHS(1:index_Roller_con-1)+convection_term*T_amb;
    
    
     RHS(index_Roller_con:N)=RHS(index_Roller_con:N)+convection_term_Roller*T_Roller;
    
    
    
    
    
    
       
    [L,U,p] = lu(CM,'vector');
    T = U\(L\(RHS(p,:)));
    
    
    if fileID_Temp_Red_Box
    Temp_Red_Box=Measuring_Box_Cal (full(T),Measure_Box_Tape,delta_x,delta_y,index_middle_long);     
          fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
    end
    
       %% For saving the Nip-point temperature matrix
       T_nip_Tape= T(BC_right);
      save('T_nip_Tape.mat','T_nip_Tape');
      %%
    if Tecplot_check_T
       Tecplot_Matlab_creator_2D (1,N,xnode,ynode,v_original,Lx,Ly,full(T),fileID_Tecplot_Tape);
     end
    %      
      
      
%      Graphical output
%     
% %     toc
    
if Graphic_Profile
    
   h21=figure(21);
   javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

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
    javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

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
  
%% outPut
if Graphic_contour
h100=figure(100);
mouse3d;
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

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
    
%     T_matrix=T_matrix';

  if Tecplot_check_T
  
            fclose(fileID_Tecplot_Tape);
   end

