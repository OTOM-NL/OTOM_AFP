
function [points_T_G,z,dl,th_v]=Tape_points (deg_T,R_tape,W_tape,L_flat,tv,Rot_Roller_axis,W_R,...
    Width_Nodes,Angular_Nodes,L_flat_Nodes,theta_ind)

% we assume th or theta is inthe degree
% Q: what is the best way to distribute the points
% node_space=0.55;
% Angular_space=3; % in degree

% dl=0:L_flat_space:L_flat;
dl=linspace(0,L_flat, L_flat_Nodes);
% dl(2)

% useless, repetition of one point
dl(1)=[];

deg=(deg_T);

if deg_T< theta_ind
    deg=(theta_ind);

end


% -3 as Tolerance
Tolerance=5;

% Having the same step size as flat part of the Tape

delta_theta=(dl(1)/R_tape)*(180/pi);
 th_v=theta_ind-Tolerance:delta_theta:deg;


% if deg_T-theta_ind >15  % set an arbitrary criteria
%     th_v=linspace(theta_ind-Tolerance,deg,Angular_Nodes);
% %     th_v=theta_ind-1:(dl(2)/R_tape):deg;  % the size of both sections should be the same
% else
%     th_v=theta_ind-Tolerance:1:deg;
% end


% th_v=1:Angular_space:th+1;

% th_v=linspace(theta_ind,deg,Angular_Nodes);
% th_v=theta_ind-1:1:deg;
th_v=th_v*(pi/180);
% z=0:node_space:W_tape;  % it should be in the middle
z=linspace(0,W_tape, Width_Nodes);
z=z-W_R/2;

% W_R width of the Roller
% z=z+W_R/4;   % it should be in the middle

z=z+(W_R-W_tape)/2;


% Making x,y of the curved part !
x=-R_tape*sin(th_v); % vector based
y=R_tape*(1-cos(th_v)); %vector based

% L_flat_space=node_space;

% % dl=0:L_flat_space:L_flat;
% dl=linspace(0,L_flat, L_flat_Nodes);


% Findidng strat point fot the Flat part!
X0_new=x(end)-dl*cosd(deg_T); %vector based
Y0_new=y(end)+dl*sind(deg_T); %vector based



% Filling for the curved part
points_T=zeros((length(th_v)+length(dl))*length(z),4);
for jj=1:length(th_v)
% for ii=1:length(z)
row_index=(1:length(z))+ (jj-1)*length(z);
% points_T(row_index,1:2)  =[x(jj) y(jj) ;
    points_T(row_index,1)  =x(jj) ;
    points_T(row_index,2)  = y(jj);
points_T(row_index,3)  = z;
% end
end

for jj=1:length(dl)
    row_index=(1:length(z))+ (jj-1+length(th_v))*length(z);
% points_T(row_index,1:2)  =[X0_new(jj) Y0_new(jj) ];
points_T(row_index,1)  =X0_new(jj) ;
    points_T(row_index,2)  = Y0_new(jj);
points_T(row_index,3)  = z;
    
end




% plot3(points_T(:,1),points_T(:,2),points_T(:,3),'r*');
% axis equal

% Now we have points for the thermal domain of the tape
%Transformation to the global axis systems
 points_T_G=Transformation_L2G (tv,Rot_Roller_axis,points_T(:,1:3)');


% plot3(points_T_G(1,:),points_T_G(2,:),points_T_G(3,:),'md');
% axis equal

% The thermal points should be saved for the calculation from here
