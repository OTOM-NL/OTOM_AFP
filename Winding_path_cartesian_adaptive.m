
%% A kinematic path on the ellipsoidal shape
% Problem: only works if z0=0 for upper dome

function [h_path,Points_center_path,Lag_points, Boarders,starting_index]=Winding_path_cartesian_adaptive(R_cyl,z_cyl_end,L_prim,W,th_y,No_dev_L,No_dev,Nip_Mov,nip_point_M);

% h_center_path >> for controlling on the main path
% 
points=[];
Boarders=[];
starting_index=[];


if length(R_cyl)==3
   
    c1=R_cyl(2);  % > Bottom dome
     c2=R_cyl(3); % > upper dome
     
else
       c1=R_cyl(1);
     c2=R_cyl(1);
%       R_cyl=R_cyl(1);
    
end

% delete previous point or line
 Delete_curves;


% number_of_Div=20;



x0=nip_point_M(1);
y0=nip_point_M(2);
z0=nip_point_M(3);



L=z_cyl_end; %L_prim;
% W=pi/20;


% m=0.05;
% th_y=acotd(m)+180;

% only Cover a semi-sphere >> Not general approach
if mod(th_y,90)==0
    th_y=th_y+0.1;
%     warningdlg('For numerical stability it adds 0.1 deg')
end

m=abs(cotd(th_y-180));

% W_m=1/m;


a=R_cyl(1);
b=R_cyl(1);
% c1=0.1;
% c2=0.12;
% phi_limit=atan((a/(c1*m))*sin(pi/2));

% Center of bottom dome
x0_E1=0;
y0_E1=0;
z0_E1=0;

% Center of upper dome
x0_E2=0;
y0_E2=0;
z0_E2=z_cyl_end;


% define the interval to find the next point
% The step should be big for smaller c1, c2 of domes
step=pi/1.01;
end_limit=step;
start_limit=0;


% delta_x from the thermal model
delta_x=L_prim/No_dev_L;

% distance value between point to be achieved
delta_s0=delta_x;

  % nonlinear finding step among points
  % calculate all points, but later choose a point for start
  
counter=0;

while (start_limit) <= pi+pi/20
    
    
  theta0=atan2(y0,x0);
    
        theta_orig=linspace(start_limit,end_limit,2);
        
%         limit=0*pi/20 ;%1*asin(1*z0/c1);
        
       Rot=0*pi/100;
        % theta0  for Rotation
        theta=theta_orig+theta0+Rot;%- (m*limit);
        
        phi=atan((a/(c2*m)).*(sin(theta_orig)));%+asin(1*z0/c1);
        
        index_start=1;

        yp=a.*cos(phi(index_start:end)).*sin(theta(index_start:end));
        xp=b.*cos(phi(index_start:end)).*cos(theta(index_start:end));
        zp=c2.*sin(phi(index_start:end)) +z0_E2 ; %.*ones(size(theta));
        
        x=xp(2);
        x0=xp(1);
        y=yp(2);
        y0=yp(1);
        z=zp(2);
        z0=zp(1);
        
       
        % derivatives
%     xp_der=-(b*sin(th)*(a^2 + c2^2*m^2))/(c2^2*m^2*((a^2*sin(th)^2 + c2^2*m^2)/(c2^2*m^2))^(3/2));
%     yp_der=(a*cos(th))/((- a^2*cos(th)^2 + a^2 + c2^2*m^2)/(c2^2*m^2))^(3/2);
%     zp_der=(a*cos(th))/(m*((- a^2*cos(th)^2 + a^2 + c2^2*m^2)/(c2^2*m^2))^(3/2));
%               
        % Newton-Raphson method >> does not work for periodic function
%      end_limit= end_limit- ( (delta_s*(delta_s-delta_s0))/...
%             (xp_der*(x-x0) + yp_der*(y-y0) + zp_der*(z-z0)));
        

fun = @(th) F_path(x0,y0,z0,a,c2,m,b,z0_E2,delta_s0,theta0+Rot,th);    % function of x alone
end_limit = fzero(fun,start_limit+[0 step]);
        

% find the result, and go for next point
    counter=counter+1;
    % we do not know the precis length of theta_orig_adapted
    theta_orig_adapted(counter)=end_limit;
    % changes for next round
           start_limit=end_limit;
           % add end_limit with a certain value
    end_limit=end_limit + step;

end

material metal; % dull; shiny;
hold on;

% [~ ,index_start]=max((theta_orig_adapted>theta0)>0);


 theta_orig=atan2(nip_point_M(2),nip_point_M(1));
 
% theta_orig_adapted=theta_orig_adapted+theta0;
% Data adapted representation
 phi=atan((a/(c2*m)).*(1*(sin(theta_orig_adapted)-0)+ 0*(m*z0/a)));%+asin(1*z0/c1);
    
        yp=a.*cos(phi(index_start:end)).*sin(theta_orig+theta_orig_adapted(index_start:end));
        xp=b.*cos(phi(index_start:end)).*cos(theta_orig+theta_orig_adapted(index_start:end));
        zp=c2.*sin(phi(index_start:end)) +z0_E2 ; %.*ones(size(theta));
             
%   plot3(xp,yp,zp,'r*');


% Now we have evenly districuted elucidian point ,called the Euclidean norm




%  figure
% plot(theta)
% hold on
% plot(phi,'r')




T_v=zeros(length(xp)-1,3);
n_xyz=T_v;
for ii=1:length(xp)-1
    tx=xp(ii+1)-xp(ii);
    ty=yp(ii+1)-yp(ii);
    tz=zp(ii+1)-zp(ii);
    
    % tangent vector for each point
    T_v(ii,:)=[tx ty tz]/(norm([tx ty tz]));
    
    % normal vector for xp,yp,zp
    
    
    if zp(ii) <=z0_E1
        
        %gradiant calculations
        nx1=(xp(ii)-x0_E1)/a^2;
        ny1=(yp(ii)-y0_E1)/b^2;
        nz1=(zp(ii)-z0_E1)/c1^2;
        
        Norm_nxyz= norm([nx1,ny1,nz1]);
        nx1=nx1/Norm_nxyz;
        ny1=ny1/Norm_nxyz;
        nz1=nz1/Norm_nxyz;
        n_xyz(ii,:)=[nx1 ny1 nz1];
        
        Center_point=[0 0 z0_E1];
       
        % normal of which surface or part is considered ?
        
        % Because we want to capture transition between cylinder and dome
    elseif  zp(ii) >=z0_E2 /2
        
        %gradiant calculations
        nx1=(xp(ii)-x0_E2)/a^2;
        ny1=(yp(ii)-y0_E2)/b^2;
        nz1=(zp(ii)-z0_E2)/c2^2;
        
        Norm_nxyz= norm([nx1,ny1,nz1]);
        nx1=nx1/Norm_nxyz;
        ny1=ny1/Norm_nxyz;
        nz1=nz1/Norm_nxyz;
        n_xyz(ii,:)=[nx1 ny1 nz1];
               Center_point=[0 0 z0_E2];
      
        
    end
    
    
end

% No_dev=5;

tape_width=W;
% vector of tape width
cp_v=cross(T_v,n_xyz);

% is a 3D matrix for storing points along long
xyz_int=zeros(3,length(cp_v),No_dev-1);
% t=linspace(0,0.5,2);

% R_cyl=[a   c1 c2];

% figure;
% hold on;


%%
% Delta_y=(1/(No_dev-1))*tape_width;

%%
% Now only works for z0_E2




% for Two side of main path
for jj=1:2
    
    
    % number of nodes for each side of width is > No_dev-1
    for ww=1:No_dev-1
        val_w=(ww/(No_dev-1))*tape_width*(-1)^jj;
        %     val_w=(ww/No_dev-1)*val(jj);
        
        
        for ii=1:length(cp_v)
            line=cp_v(ii,1:3)*val_w ; %+ ([xp(ii) yp(ii) zp(ii) ]);
            
            point=[Center_point;xp(ii) yp(ii) zp(ii) ];
            
            point_w=line+point;
            
            % find intersection of point_wwith sphere, for width of the tape
            % plot3(point_w(:,1),point_w(:,2),point_w(:,3),'c');
            
            [xyz_int(1:3,ii,ww)]=Line_Mandrel_intersection_4subW(point_w(1,:)',[xp(ii) yp(ii) point(2,3)-point(1,3)]',R_cyl,z0_E2);
            
            % save points
            % thermal_points(1:3,(1:No_row_nodes) + (pp-1)*No_row_nodes )=xyz_int
            
            % plot3(point_w(:,1),point_w(:,2),point_w(:,3),'c');
            %  plot3(xyz_int(1,:),xyz_int(2,:),xyz_int(3,:),'c*');
            
        end
        
    end
    
    
    if jj==1
        % points rows in upper side
        up_point=xyz_int;
        index_up=find(xyz_int(3,:,end)< z0_E2);
    elseif jj==2
        bot_point=xyz_int;
        index_bot=find(xyz_int(3,:,end)< z0_E2);
    end
    
    % plot3(xyz_int(1,:),xyz_int(2,:),xyz_int(3,:),'c*');
    
end


[Mutual,st1,st2] = intersect(index_up, index_bot);

% the final index range that after it, we should use helical path

if ~isempty (Mutual)
    final_index0=1:min(Mutual);
    
else
    final_index0=1:length(up_point);
    
end


% plot3(xp(final_index),yp(final_index),zp(final_index),'r.-');
% X_width_point=(1,final_index,ii)

X_width_point=[];%zeros(1,(No_dev-1)*final_index);
Y_width_point=[]; %zeros(1,(No_dev-1)*final_index);
Z_width_point=[];%zeros(1,(No_dev-1)*final_index);

for ii=1:No_dev-1
    X_width_point=[X_width_point up_point(1,final_index0,ii) bot_point(1,final_index0,ii) ];
      Y_width_point=[Y_width_point up_point(2,final_index0,ii) bot_point(2,final_index0,ii) ];
        Z_width_point=[Z_width_point up_point(3,final_index0,ii) bot_point(3,final_index0,ii) ];
        
%     plot3(up_point(1,final_index,ii),up_point(2,final_index,ii),up_point(3,final_index,ii),'c*');
%     plot3(bot_point(1,final_index,ii),bot_point(2,final_index,ii),bot_point(3,final_index,ii),'c*');
end

% for performance we plot all one-time
%  plot3(X_width_point,Y_width_point,Z_width_point,'c.');
%     plot3(bot_point(1,final_index,ii),bot_point(2,final_index,ii),bot_point(3,final_index,ii),'c*');

% Start position for helical winding
% xp(min(Mutual)),yp(min(Mutual)),zp(min(Mutual))

Main_path_cyl=zeros(3,0);
thermal_points=zeros(3,0);

if ~isempty (Mutual)
    start_cyl=(up_point(1:3,final_index0(end))+bot_point(1:3,final_index0(end)))/2;
    
%     plot3(start_cyl(1),start_cyl(2),start_cyl(3),'b*');
    
    
    
    % R=a;
%     L=a*1;
%     th_y=acotd(m)+180;
%     No_dev_L=20;
    
    Nip_Mov=0;
    Length_on_cylinder=1*L;
    
    
 
    
    
    fun_cyl = @(No_dev_L) F_path_cyl(Length_on_cylinder,start_cyl,delta_s0,th_y,R_cyl,No_dev_L);    % function of x alone
No_dev_L_adapted = fzero(fun_cyl,0+[2 1e3 ]);
    
    No_dev_L_adapted=ceil(No_dev_L_adapted);
    
    
    [Main_path_cyl,thermal_points,Ps]=helical_3D_points_4subW(R_cyl,z0_E2,Length_on_cylinder,tape_width,th_y,No_dev_L_adapted,No_dev,Nip_Mov,start_cyl);
%     plot3(Ps(1,:),Ps(2,:),Ps(3,:),'y*')
final_index=final_index0(1:end-1);
else
     final_index=final_index0;
end
% for performance we plot all one-time
%  plot3(X_width_point,Y_width_point,Z_width_point,'c.');
h_Lag_points=plot3([X_width_point thermal_points(1,:)],[Y_width_point thermal_points(2,:)],[Z_width_point thermal_points(3,:)],'g.');


% They have one mutual point, that is why starts at 2nd index
h_center_path=plot3([xp(final_index) Main_path_cyl(1,1:end)],[yp(final_index) Main_path_cyl(2,1:end)],[zp(final_index) Main_path_cyl(3,1:end)],'r.-');

Points_center_path=[xp(final_index) Main_path_cyl(1,1:end);...
    yp(final_index) Main_path_cyl(2,1:end);...
    zp(final_index) Main_path_cyl(3,1:end)];


h_path=[h_center_path,h_Lag_points];


% re-arranging thermal nodes along length direction
W_step=2*No_dev-1;
[s1,s2]=size(thermal_points);
% No_node_along_length v>NAL
NAL=s2/W_step;

cyl_point=zeros(3,NAL,W_step);
for ii=1:W_step
cyl_point(:,:,ii)=thermal_points(:,(1+(ii-1)):W_step:s2) ;
end

% length(up_point(1,final_index,1)) thermal_points(1,:)
%% save lagrangian points all-to-gether with a specific pattern to access
Length_path=length([xp(final_index) Main_path_cyl(1,1:end)]);
Lag_points=zeros(3,Length_path,W_step);

counter=1;
for ii=No_dev-1:-1:1
    Lag_points(1:3,:,counter)=[up_point(1,final_index0,ii) cyl_point(1,:,counter) ;...
       up_point(2,final_index0,ii) cyl_point(2,:,counter) ;...
       up_point(3,final_index0,ii) cyl_point(3,:,counter) ];
    
counter=counter+1;
        
end
    % in the middle red points are coming

 Lag_points(1:3,:,counter)=[xp(final_index) Main_path_cyl(1,1:end);...
     yp(final_index) Main_path_cyl(2,1:end); ...
     zp(final_index) Main_path_cyl(3,1:end)];

 
 counter=counter+1;
 
for ii=1:No_dev-1
   Lag_points(1:3,:,counter)=[bot_point(1,final_index0,ii) cyl_point(1,:,counter) ;...
       bot_point(2,final_index0,ii) cyl_point(2,:,counter) ;...
       bot_point(3,final_index0,ii) cyl_point(3,:,counter) ];
    
counter=counter+1;
end

% Now Lagrangian nodes are along length of curved path on a 3D matrix,
% 3rd dimention is about the width nodes

% bot_point(1,final_index,ii)

