
function  [points_New, Boarders,starting_index]=sub_on_ellips(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv)


 %% List of inputs
% z_cyl_end distance from the axis center
% L_prim length of the tape
% w width of the Tape
% sui orientation of the substrate
% No_dev_L devision of the substrate along length
% No_dev devision along half width
% Nip_Mov index before the end because of the deformation
% tv movement from the axis center

% a,b, c are the radius of the ellipsoidic shape



starting_index=floor(No_dev_L*(Nip_Mov/L_prim)); 


points=zeros(3,No_dev_L* (2*No_dev-1));
points_New=points;


if length(R_cyl)==3
   
    c1=R_cyl(2);
     c2=R_cyl(3);
else
       c1=R_cyl(1);
     c2=R_cyl(1);
%       R_cyl=R_cyl(1);
    
end



%%
a=R_cyl(1);
b=R_cyl(1);
% c1=3;

% c2=1;

wind_angle=sui*pi/180;

wind_angle=wind_angle- (.001);  % to avoid numerical instability



%  dcm_obj = datacursormode(gcf);
% c_info = getCursorInfo(dcm_obj);
% 
%  P_cursor=c_info.Position;

x0=tv(1);
y0=tv(2);

if tv(3)> z_cyl_end
    % 2nd Dome  
    
z0=tv(3)-z_cyl_end;
c=c2;
else
%     1st dome

 z0=tv(3);   
c=c1;
end

% x0=P_cursor(1);
% y0=P_cursor(2);
% z0=P_cursor(3);

L=L_prim;
W=w;

 theta0=atan2(x0,y0);
    phi0=asin(z0/c);


wind_angle_lateral=wind_angle+pi/2;

[th_limit12,status]=ellipsoid_path_integration(a,b,c,wind_angle_lateral,W,theta0,phi0);







number_of_Div=No_dev;
xyz=zeros(3,2*number_of_Div-1);

for ii=1:2
    
th_limit=real(th_limit12(ii)*((-1)^(ii-1)));   % for other side of Tape

% th_limit=pi/200;

theta0=atan2(x0,y0);
% phi=asin(z0/c)

theta=linspace(theta0,theta0+th_limit,number_of_Div);



theta_range=linspace(0,th_limit,number_of_Div);
% theta=linspace(pi/2,1*pi,160);
phi=asin(z0/c)+ tan(wind_angle_lateral)*theta_range;
% phi=pi/6;


% phi0=asin(z0/c);

% first trajectory
% xp=a.*cos(0).*sin(theta);
% yp=b.*cos(0).*cos(theta);
% zp=c.*sin(phi); %.*ones(size(theta));


% second trajectory
xp=a.*cos(phi).*sin(theta);
yp=b.*cos(phi).*cos(theta);
zp=c.*sin(phi); %.*ones(size(theta));


% make 2 semi-line in opposite direction

if ii==1
xyz(1:3,((ii-1)*number_of_Div)+(1:number_of_Div))=[xp(end:-1:1);yp(end:-1:1);zp(end:-1:1)];

elseif ii==2
 xyz(1:3,((ii-1)*number_of_Div)+(1:number_of_Div-1))=[xp(2:end);yp(2:end);zp(2:end)];

end


% plot3(xp (zp>0),yp (zp>0),zp(zp>0),'r.-');
% plot3(xp (zp<0),yp (zp<0),zp(zp<0),'b.-');

% if tv(3)> z_cyl_end
% plot3(xp,yp,z_cyl_end+zp,'w.-');
% else
%   plot3(xp,yp,zp,'w.-');  
% hold on;
% end

end

% tp prevent repition first index is ignored
% xyz=xyz(1:3,2:end);


% 
% [th_limit12,status]=ellipsoid_path_integration(a,b,c,wind_angle,L,theta0,phi0);
% 
%% Index numbering in previoius version was adjusted based on end to begining along the width !!
xyz(:,:)=xyz(:,end:-1:1);
%%

% for ii=1:2

    
%     xyz(1:3,((ii-1)*number_of_Div)+(1:number_of_Div))
xp=xyz(1,:);
yp=xyz(2,:);
zp=xyz(3,:);

% Borader nodes
% xp=xyz(1,[number_of_Div-1, end]);
% yp=xyz(2,[number_of_Div-1, end]);
% zp=xyz(3,[number_of_Div-1, end]);

% boarder = xyz(1:3,[number_of_Div-1, end])



% number_of_Div=25;

counter=0;
Boarders=zeros(3,No_dev_L*2);


for pp=1:length(xp)
    
    
    
    theta0=atan2(xp(pp),yp(pp));
    phi0=asin(zp(pp)/c);
    
    
    
[th_limit12,status]=ellipsoid_path_integration(a,b,c,wind_angle,L,theta0,phi0);


    
th_limit=th_limit12(1);%*((-1)^(ii-1));   % for other side of Tape

% th_limit=pi/200;

theta0=atan2(xp(pp),yp(pp));
% phi=asin(z0/c)



theta=linspace(theta0,theta0+th_limit,No_dev_L);



theta_range=linspace(0,th_limit,No_dev_L);
% theta=linspace(pi/2,1*pi,160);
phi=asin(zp(pp)/c)+ tan(wind_angle)*theta_range;
% phi=pi/6;


% first trajectory
% xp=a.*cos(0).*sin(theta);
% yp=b.*cos(0).*cos(theta);
% zp=c.*sin(phi); %.*ones(size(theta));


% second trajectory
xTraj=a.*cos(phi).*sin(theta);
yTraj=b.*cos(phi).*cos(theta);
zTraj=c.*sin(phi); %.*ones(size(theta));



% xyz(1:3,pp*length(xTraj))=[xTraj;yTraj;zTraj];


% plot3(xp (zp>0),yp (zp>0),zp(zp>0),'r.-');

% plot3(xp (zp<0),yp (zp<0),zp(zp<0),'b.-');

% just for representation !!
% if tv(3)> z_cyl_end
% plot3(xTraj,yTraj,z_cyl_end+zTraj,'w.-');
% else
%     plot3(xTraj,yTraj,0+zTraj,'w.-');
% end
% 
% hold on;


points(1:3,1+(pp-1)*No_dev_L:(pp)*No_dev_L )=[xTraj;yTraj;zTraj];




if pp==1 
    
    counter=counter+1;
Boarders(1:3,1+(counter-1)*No_dev_L:(counter)*No_dev_L)=[xTraj;yTraj;zTraj];

elseif  pp== length(xp)
    counter=counter+1;
    Boarders(1:3,1+(counter-1)*No_dev_L:(counter)*No_dev_L)=[xTraj(end:-1:1);yTraj(end:-1:1);zTraj(end:-1:1)];
    
% boarder = xyz(1:3,[number_of_Div-1, end])
end



end

if tv(3)> z_cyl_end
Boarders(3,:)=Boarders(3,:)+z_cyl_end;
points(3,:)=points(3,:)+z_cyl_end;

end


%% To make in a good order- width order
% because other part of codes, thermal models, works in this order

% points_New=points;

No_dev_W=2*No_dev-1;




for ii=1:No_dev_L
points_New(:,(1:No_dev_W)+ (ii-1)*No_dev_W)=points(:,((1:No_dev_W)-1)*No_dev_L + ii);

end


points_New=points_New(:,No_dev_W*starting_index+1:end);

% end

