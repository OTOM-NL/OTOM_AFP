
function [Rot_Roller_axis,beta_macro,S_normal_active_G,kr1_active_G]=Rot_Matrix_Finder_local_BRDF_Tape(BRDF_Tape,n_Macro,laser_direction)

Fib_Rot= BRDF_Tape(1);
div= BRDF_Tape(2);
sig_T= BRDF_Tape(3);
sig_F= BRDF_Tape(4);
Amp= BRDF_Tape(5);
threshold= BRDF_Tape(6);



%  div=20; % Sample data from USER
%
%    sig_T=0.1; % x from USER
%     sig_F=0.5; % Y from USER
%
%   Amp=35; %from USER
%      threshold=20; %from USER- a value which is worth to consider

%% To avoid Numerical issues
    n_Macro(n_Macro==0)=1e-7;
     laser_direction(laser_direction==0)=1e-7;


%% direction for determining phi, respect to fibers
project_ray=-(laser_direction-dot(laser_direction,n_Macro)*n_Macro' );

project_ray=project_ray./norm(project_ray);

% t=linspace(0,0.5,2);

% figure(1);
%   Pr_v(1,:)=project_ray(1)*t+xyz_int(1);
%                 Pr_v (2,:)=project_ray(2)*t+xyz_int(2);
%                      Pr_v (3,:)=project_ray(3)*t+xyz_int(3);
%                                   plot3(Pr_v(1,:),Pr_v(2,:),Pr_v(3,:),'b--');

%                                   Fib_Rot=45; % in degree

% phi-refference angle of the tape, direction with z-axis in local axis
% phi=acosd(project_ray);


Fib_or=0+Fib_Rot; % repect to the tangent to the helical direction


% Vec_Pr_v=[Pr_v(1,1)-Pr_v(1,2),Pr_v(2,1)-Pr_v(2,2),Pr_v(3,1)-Pr_v(3,2)];
Vec_Pr_v=-project_ray;

% if R_cyl(1)>0
% for cylinder and dome

% [Rot_Roller_axis,phi_Fiber]=Rot_Matrix_Fiber_BRDF(xyz_int,Fib_or,R_cyl,z_cyl_end,Vec_Pr_v,n_Macro');

% else % for tape placement
% fiber direction - z-axis
Vec_tang=[0 ;0 ;1]; % refference axis to measure with


vector_cross_Pr_v=cross(Vec_Pr_v,Vec_tang);


% should be Rotated
theta_normal_cross_product_Pr_v=acosd(dot(vector_cross_Pr_v,n_Macro)/(norm(vector_cross_Pr_v)*norm(n_Macro)));


phi_2vec= acosd(dot(Vec_Pr_v,Vec_tang)/(norm(Vec_Pr_v)*norm(Vec_tang)));






%   phi_0=atand(project_ray(2)/project_ray(1));
% =0-(sign(cosd(theta_normal_cross_product_Pr_v))*phi_0)+


phi_Fiber=-sign(cosd(theta_normal_cross_product_Pr_v))*phi_2vec+ 0;%0*90; %+Fib_Rot+0;

%    phi_Fiber=-phi_2vec+ 180;


%     Rot_z
%     Rot_Roller_axis=([cosd(-Fib_or)  -sind(-Fib_or) 0 ;  sind(-Fib_or)  cosd(-Fib_or) 0 ; 0 0 1]);
%      phi_Fiber=0-(sign(cosd(theta_normal_cross_product_Pr_v))*acosd((Pr_v*Vec_tang')/(norm(Vec_tang)*norm(Pr_v))));

% end



% Rot_Roller_axis_Arb  Arbitrary

v1=[0 ;0  ;1];   % the initial laser head on xy plane
%  v2=[0 1 0]';  % normal of the surface
v2=n_Macro';
Rot_Roller_axis_Arb=(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix
%

% CCW_rot=[0 -1 ; 1 0];

% Vec_tang=CCW_rot*n_Macro;

% axes X of BRDF goes to Z-axis in local tape-roller axis

Vec_Xaxis_roller=[1 ;0 ;0];

Vec_X_rotated=Rot_Roller_axis_Arb*Vec_Xaxis_roller;

theta_2Refconfig=acosd(dot(Vec_X_rotated,Vec_tang)/(norm(Vec_tang)*norm(Vec_X_rotated)));



theta_normal_cross_product=cross(Vec_X_rotated,Vec_tang);




theta_normal_cross_product_CW=acosd(dot(theta_normal_cross_product,n_Macro)/(norm(theta_normal_cross_product)*norm(n_Macro)));


Fib_or=-sign(cosd(theta_normal_cross_product_CW))* theta_2Refconfig+Fib_or+180;

%  end


% Make a Rotation that transform axis to final axis
% Rot_y=[cosd(Fib_or) 0 sind(Fib_or) ; 0 1 0; -sind(Fib_or) 0 cosd(Fib_or)];

Rot_z=([cosd(-Fib_or)  -sind(-Fib_or) 0 ;  sind(-Fib_or)  cosd(-Fib_or) 0 ; 0 0 1]);





Rot_Roller_axis=Rot_Roller_axis_Arb*Rot_z;








%%
% Use this rotation to transform Icoming Ray to local axis of BRDF
% intersection point goes to [0 0 0]

%  Incom_Ray_BRDF=Rot_plane*laser_direction;


%  n=100;  % number of surface devision

% div=20; % Sample data from USER


beta_macro=90-(1*acosd(dot(laser_direction,-n_Macro)/(norm(laser_direction)*norm(n_Macro))));

% Angle of incidence , Elevation angle
% theta=90-beta_macro;  % got from macro calculation , Angle incidence respect to macro normal
theta=beta_macro;

% phi=0;   % got from macro calculation , fiber orientation respect to incoming ray

% BRDF figure
%     figure(15);
%     cla;

fib_or=(0)*pi/180;  % fibre orientation from USER, 0 along helical path,cross-section

% making sample points along the fiber
point_th=linspace(0,pi/2,div);
point_phi=linspace(-fib_or,-fib_or,div); %half round
point_phi2=point_phi+pi;

%     [xo,yo,zo]=sph2cart([point_phi point_phi2(2:end)],[point_th point_th(end-1:-1:1) ],1); % R=1, directions are important


%     figure(15);
%     hold on;
%
%     plot3(xo,yo,zo,'g*');

%     sig_T=0.1; % x from USER
%     sig_F=0.5; % Y from USER
%
%
%     Amp=35; %from USER



%% decide which points are worth to consider
% integral of the points  shold be done later

%% >>>>>>>>>>>>  Sarting from here 4 OTOM  <<<<<<<<

%      [xo,yo,zo]=sph2cart([point_phi point_phi2(2:end)],[point_th point_th(end-1:-1:1) ],1); % R=1, directions are important

ALLpoint_th=[point_th point_th(end-1:-1:1) ];
ALLpoint_phi=[point_phi point_phi2(2:end)]+fib_or;

No_all=2*div-1; % number of all sample points along fiber
%     threshold=20; %from USER- a value which is worth to consider

counter=0; % number of points to be considered
Selected_point=zeros(No_all,2);
weight_microfacet=zeros(No_all,1);

for jj=1:No_all  % correspoding weight values/ density of microfacet
    p_point =Amp*exp(   -(1)*tan(ALLpoint_th(jj)+pi/2).^2 * ( cos(ALLpoint_phi(jj)).^2 /(2*sig_F^2) + sin (ALLpoint_phi(jj)).^2 /(2*sig_T^2))   );
    
    if  p_point  >threshold  % enough possible energy to transfer?
        counter=counter+1;
        weight_microfacet(counter)=p_point;
        Selected_point(counter,1:2)=[ALLpoint_th(jj),ALLpoint_phi(jj)];
    end
    
end

% Remove th rest which were not filled
weight_microfacet=weight_microfacet(1:counter);
Selected_point=Selected_point(1:counter,1:2);

% phi comes first, bring to local cartesian
[x_sel,y_sel,z_sel]=sph2cart(Selected_point(:,2)-fib_or,Selected_point(:,1),1);


%     figure(15);
%     hold on;
%show selected points
%     plot3(x_sel,y_sel,z_sel,'rd','MarkerSize',10);




%% a sample ray  As an example
% reflection calculations
%     figure(15);

% move intersection point to the [0 0 0]
% Transform macro normaL surface to [0 0 1]



% Azimuthal angle, relative angle of incidence into the fibers
phi=Fib_Rot+phi_Fiber;%-Fib_or;

[xp,yp,zp]=sph2cart(phi*pi/180,theta*pi/180,1);

%     xp = 1 * cos(theta) * cos(phi);
%     yp = 1* cos(theta) * sin(phi );
%     zp = 1 *sin(theta);

% plot source in local axis
%       plot3(xp,yp,zp,'r*');


Start_point=[xp yp zp];
n_ray=-Start_point; % direction > Incoming Ray direction
n_ray=n_ray/norm(n_ray); %normalizartion
%      n_ray=-Incom_Ray_BRDF';  % in local coordinate

% plot incoming ray
%    vector=[xp yp zp; 0 0 0];
%   plot3(vector(:,1),vector(:,2),vector(:,3),'r');



% angle of incident
%     beta=zeros(1,length(x_sel));
counter_kr1=0;

%        kr1_active=[];


% the angle between incoming ray and normal, should be surface normal of macro
% model
% How much ebergy is absorbed > micro-model does not involved
%             S_normal_macro=[0 0 1]; % Needs attention
% should be calculated before
%             beta_macro=1*acosd(dot(n_ray,S_normal_macro)/(norm(n_ray)*norm(S_normal_macro)));


S_normal_active=zeros(length(x_sel),3);
kr1_active=zeros(length(x_sel),4);


for kk=1:length(x_sel)
    
    
    % local the point of intersection
    %         x_int=0;
    %         y_int=0;
    %         z_int=0;
    
    % ray direction in local coordinate > Transformation required
    %          Rx=n_ray(1);
    %         Ry=n_ray(2);
    %         Rz= n_ray(3);
    
    
    % to show reflected ray
    %         t=linspace(0,1,2);
    
    % normal of micro facet surface, originating from [0 0 0]
    S_normal=[x_sel(kk),y_sel(kk),z_sel(kk)];
    % Ref_line=zeros(1,3);
    
    
    % direction of the reflection
    kr1=zeros(1,3);
    
    
    % if S_normal(1,1)
    
    kr1(1,1:3)=2*(-n_ray*S_normal(1,:)')*S_normal(1,:)-(-n_ray);
    %kr1 direction should be outward of the surface !
    
    
    % check whther it is visible from viewer and objective, upper
    % hemisphere reflection and incoming ray
    
    if kr1(3)>0
        
        % count how many are considered
        counter_kr1=counter_kr1+1;
        
        
        %plot microfacet involved normals, can be reported as results
        %as well
        
        %Save in .txt file > Transform to Global coordinate later
        %             vector_normal=[0 0 0;S_normal];
        
        
        %           plot3(vector_normal(:,1),vector_normal(:,2),vector_normal(:,3),'y--');
        
        %Save txt Data > Use as normals
        % Local normal lines
        
        S_normal_active(counter_kr1,1:3)=S_normal;
        
        kr1_active(counter_kr1,1:3)=kr1;  % just for output in general program
        kr1_active (counter_kr1,4)=weight_microfacet(kk); % bring related energy density related to microfacet
        
        
        % reflection line for representation
        %             Ref_line_X(1,:)=kr1(1,1)*t+x_int(1);
        %            Ref_line_Y (1,:)=kr1(1,2)*t+y_int(1);
        %            Ref_line_Z (1,:)=kr1(1,3)*t+z_int(1);
        
        %             plot3(Ref_line_X(1,:),Ref_line_Y(1,:),Ref_line_Z(1,:),'k--')
        
        % text(Ref_line_X(1,2),Ref_line_Y(1,2),Ref_line_Z(1,2) , sprintf('Ref'));
        %Save txt Data
        
        % Local reflections
        
        %             Ref_line=[Ref_line_X;Ref_line_Y;Ref_line_Z];
        %              kr1_active (1:counter_kr1,4) > to indicate energy
        
        % Active microfacet directions and corresponding energy
        % transfer capability, related to microfacet density
        % orinetations
        
        
        
        %              beta(kk)
    end
    
    
end





% Remove th rest which were not filled
S_normal_active=S_normal_active(1:counter_kr1,1:3);
kr1_active=kr1_active(1:counter_kr1,1:4);


%% Transform to the Global axis



% Make a Rotation that transform axis to final axis
% Rot_y=[cosd(th_y) 0 sind(th_y) ; 0 1 0; -sind(th_y) 0 cosd(th_y)];

% The rotation matrix from initial Config to Final Config of coordinate
% axis






%  v1=[0  ;0 ;1];   % the initial laser head on xy plane
% %  v2=[0 1 0]';  % normal of the surface
%  v2=-n_Macro';
%  Rot_plane=transpose(fcn_RotationFromTwoVectors(v1, v2)); % find the rotation matrix

%  z_rot=50;
% Rot_z=([cosd(-z_rot)  -sind(-z_rot) 0 ;  sind(-z_rot)  cosd(-z_rot) 0 ; 0 0 1]);



% Rot_Roller_axis2=(Rot_plane)*Rot_z;

% Transorm to local to local axis

% Send out the global BRDF microfacet normals and its reflections
kr1_active_G=Rot_Roller_axis*kr1_active(:,1:3)';

% to show intensity of a normals, percentage relation
kr1_active_G(4,:)=kr1_active(:,4)'./sum(kr1_active(:,4));
S_normal_active_G=Rot_Roller_axis*S_normal_active(:,1:3)';



%  figure(1);
%         t=linspace(0,1,2);
%
%
%   X(1,:)=n_Macro(1)*t+xyz_int(1);
%                 Y (1,:)=n_Macro(2)*t+xyz_int(2);
%                      Z (1,:)=n_Macro(3)*t+xyz_int(3);
%                                   plot3(X,Y,Z,'k');
%
%                      Tr_BRDF=Rot_Roller_axis*v1;
%
%                        X(1,:)=Tr_BRDF(1)*t+xyz_int(1);
%                 Y (1,:)=Tr_BRDF(2)*t+xyz_int(2);
%                      Z (1,:)=Tr_BRDF(3)*t+xyz_int(3);
%
%
%                 plot3(X,Y,Z,'c--');

% xyz_int,fib_or,S_normal


