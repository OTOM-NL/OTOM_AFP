


function T_matrix=thermal_advection_comb_FEM(Lx,Ly,xnode,ynode,E_points,v_original,materials,Temp_Incoming,h,thickness,Rel_v_theta,Measure_Box_Sub,...
                    fileID_Temp_Red_Box,Tecplot_check,Graphic_profile,Grpahic_contour)



%% initialization 
fileID_Tecplot= fopen('.\Temp_2D_Sub_Steady.plt','w');

% tic ;

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

T_BC=sparse(N,1);

%%
T_amb=20;  % Temperature for convection from surface, not edge!!

%%



T_right=0;
T_left=Temp_Incoming;  % incoming velocity
T_top=0;
T_bottom=0;
        
BC_Left=1:ynode;
BC_top=ynode:ynode:N;
BC_bottom=1:ynode:N-1;
BC_right=N:-1:N-ynode+1 ;

indexxy=[N,N-ynode+1, 1 ynode];


for kk=1:xnode-1
        index_middle_long(kk)=ceil(ynode/2)+(ynode*kk);
end
% Define a nip-point line
row_nip=0; % number before the last row
index_nip_point=(ynode*(xnode-row_nip-1)+1):((xnode-row_nip)*(ynode));






% Heat flux is constant here
g_dot=1/thickness;
g_dot=g_dot/(delta_x*delta_y);  % for the lement size
g_dot=g_dot/(xnode*ynode);  % total energy devided per number of nodes for each node

g_dot=ones(N,1)*g_dot;
g_dot([BC_bottom BC_top])=g_dot([BC_bottom BC_top])*2;
g_dot([BC_right BC_Left])=g_dot([BC_right BC_Left])*2;








 CM=sparse(N);
    %         CM=zeros(N);
    
    for kk=1:N
        CM(kk,kk)=1;
     end
  RHS=sparse(N,1);
        RHS(1:end)=-E_points.*g_dot*N;

 v0=ro*cp*v_original;
 
    v=v0*abs(cosd(Rel_v_theta));
    v_2=v0*abs(sind(Rel_v_theta));
 
 
  pec=(v_original*delta_x)/alpha;

%    h1=h(1);
%     h2=h(2);
  
  
  convection_term=-h*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
        convection_term=convection_term/(delta_x*delta_y*thickness);   % (delta_x*delta_y*thickness); volume of the element, now convection_term is W/m^3

        
%          convection_term_Roller=-h*A/(xnode*ynode);  % A/(xnode*ynode); is the area of each element !
%     convection_term_Roller=convection_term_Roller/(delta_x*delta_y*thickness); 
    
    
%  diag_term=-2*k*((1/delta_x^2)+(1/delta_y^2))-(v/delta_x)+convection_term;
   diag_term=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)-(v_2/delta_y)+convection_term;
  
%     diag_term_Roller_contact=-2*k*(1/delta_x^2)-2*ky*(1/delta_y^2)-(v/delta_x)+convection_term_Roller;
    
%     index_Roller_con=1+(Row_numberape_Roller_tangent*ynode);
 

        CM= CM*diag_term;

     

  %define temp terms
  
       C_ip1_j=(k/delta_x^2);
    C_im1_j=(k/delta_x^2)+(v/delta_x);
    C_i_jp1=(ky/delta_y^2);
    C_i_jm1=(ky/delta_y^2)+(v_2/delta_y);
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
      
        
    %% for continous boundary
diag_term_CB_y=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term;
% diag_term_CB_y_after_roller=-2*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term_Roller;

for ii=1:length(BC_top)
%     if BC_top(ii) <index_Roller_con
        CM(BC_top(ii),BC_top(ii))=diag_term_CB_y;
        CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y;
        
%     else
%         CM(BC_top(ii),BC_top(ii))=diag_term_CB_y_after_roller;
%         CM(BC_bottom(ii),BC_bottom(ii))=diag_term_CB_y_after_roller;
%         
  
    
end

% diag_term_CB_xR_right=-1*k*(1/delta_x^2)-2*ky*(1/delta_y^2) -(v/delta_x);
diag_term_CB_xR=-1*k*(1/delta_x^2)-2*ky*(1/delta_y^2) -(v/delta_x)+convection_term;

%     T_right=0;

for ii=1:length(BC_right)
    
    CM(BC_right(ii),BC_right(ii))=diag_term_CB_xR;  % for right side
    CM(BC_Left(ii),BC_Left(ii))=diag_term_CB_xR;  % for right side
    
end

%     for four corners
diag_term_CB_xy=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term;
% diag_term_CB_xy_after_roller=-1*k*(1/delta_x^2)-1*ky*(1/delta_y^2) -(v/delta_x)+convection_term_Roller;
%     for four corners
for ii=1:4
    %          index_xy=indexxy(ii);
    CM(indexxy(ii),indexxy(ii))=diag_term_CB_xy;  %  ??
% end
% for ii=3:4
    %          index_xy=indexxy(ii);
    CM(indexxy(ii),indexxy(ii))=diag_term_CB_xy;  %  ??
end


    
    
    
%      Advection_term_diag=-(v/delta_x);
%     
%     for kk=1:N
%         CM_C(kk,kk)=CM_C(kk,kk)+Advection_term_diag;
%     end
%     
    
    
    RHS(BC_Left)=RHS(BC_Left)-((0*k/delta_x^2)+(v/delta_x))*T_left;
%     RHS(BC_bottom)=RHS(BC_bottom)-((ky/delta_y^2)+(v_2/delta_y))*T_bottom;
%     RHS(BC_top)=RHS(BC_top)-(ky/delta_y^2)*T_top;
%     RHS(BC_right)=RHS(BC_right)-(k/delta_x^2)*T_right;
  
    
    
    RHS=RHS+convection_term*T_amb;
    
       
    [L,U,p] = lu(CM,'vector');
    T = U\(L\(RHS(p,:)));
    
%     toc
    
    
      if fileID_Temp_Red_Box
    Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box_Sub,delta_x,delta_y,index_middle_long);
 fprintf(fileID_Temp_Red_Box,'   %12.8f  \r\n',Temp_Red_Box);
      end
    
      
      
        if Tecplot_check
       Tecplot_Matlab_creator_2D (1,N,xnode,ynode,v_original,Lx,Ly,full(T),fileID_Tecplot);
     end
      
    
    
    if Graphic_profile
    
   h21=figure(21);
   mouse3d;
   javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

     set(h21,'Visible','off');
     
    subplot(2,2,1);
    hold on;
    % title('Tmepereture along middle of the tape');
    title(sprintf('Sub-Temperature along length, V= %f , Pec: %f, h=%d ',v_original,pec,h));
      ylabel('Temperature ^{\circ}C');
     xlabel('Length (m)');
    
      H_long=linspace(0,Lx,length(T(index_middle_long))); % H_long > horizental 
        plot(H_long,T(index_middle_long));
        
             pos=get(gca,'Position');
        
        x = pos(1)+[pos(3)/4 3*pos(3)/4];
y = pos(2)+[pos(4)/2 pos(4)/2];
annotation('textarrow',x,y,'String','Velocity','Color','k','Linewidth',2);


    
    h21=figure(21);
    
      set(h21,'Visible','off');
      
    subplot(2,2,2);
    hold on;
      ylabel('Temperature ^{\circ}C');
     xlabel('Width (m)');
    % title('Tmepereture along width nip-point of the tape');
    title(sprintf('Sub-Temperature along width nip-point, V=%f , Pec=%f, h=%d ',v_original,pec,h));
    
%     H_Nip=linspace(0,Ly,length(T(index_nip_point))); %  H_Nip > horizental 
%     plot(H_Nip,T(index_nip_point));
    
       T_nip_point=T(index_nip_point);
     H_Nip=linspace(0,Ly,length(T_nip_point)); % H_Nip > horizental 
    plot(H_Nip,T_nip_point);
    Nip_point_Temp_Sub=T_nip_point(floor(length(index_nip_point)/2));
    
    assignin('base','Nip_point_Temp_Sub',Nip_point_Temp_Sub);
    
    end
  
    
    
%% outPut contour
if Grpahic_contour
h100=figure(100);
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));


  set(h100,'Visible','off');
  
hold on;
subplot(2,1,1);

     
    T_matrix2=reshape(T,ynode,xnode);
    
    [X_2D,Y_2D]=meshgrid(linspace (0,Lx,xnode),linspace(0,Ly,ynode) );
       
    surf(X_2D,Y_2D,T_matrix2);
       colorbar;
     title('Temperature distribution of Substrate');
     shading flat;
     view([0 -90]);
     
          pos=get(gca,'Position');
                x = pos(1)+[pos(3)/4 pos(3)/1];
y = pos(2)+[pos(4)/2 pos(4)/2];
annotation('textarrow',x,y,'String','Velocity','Color','w','Linewidth',2);

if abs (v_2) > 1
 x = pos(1)+[pos(3)/4 pos(3)/4];
y = pos(2)+[pos(4)/2 pos(4)/4];
annotation('textarrow',x,y,'String','Velocity-Y','Color','r','Linewidth',2);
end


end




    len=length(T);
    T_disp=zeros(len,1);
    
    
    counter=0;
    
%     for ii=1:ynode
%         for jj=1:xnode
%             counter=counter+1;
%             index=(ii*xnode) -jj+1;
%             T_disp(index)=T(counter);
%             
%         end
%     end
% %     
  
% N=xnode*ynode;
    for ii=xnode:-1:1
        for jj=1:ynode
            counter=counter+1;
            index=(jj)+((ii-1)*ynode);
            T_disp(index)=T(counter);
            
        end
    end
    

    
    T_matrix=reshape(T_disp,ynode,xnode);
    
%     T_matrix=T_matrix';



  if Tecplot_check
  
            fclose(fileID_Tecplot);
   end



