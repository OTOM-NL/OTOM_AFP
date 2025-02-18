function  [points_New, Boarders,starting_index]=sub_on_flat(R_cyl,z_cyl_end,L_prim,w,sui,No_dev_L,No_dev,Nip_Mov,tv)



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


No_dev_W=2*No_dev-1;

points=zeros(3,No_dev_L* No_dev_W);
points_New=points;



x0=tv(1);
y0=tv(2);

z0=0;

 x_limit=w* cosd(sui+90);
 y_limit=w* sind(sui+90);


 x_start=linspace (x0-x_limit ,x0+x_limit, No_dev_W );
 y_start=linspace (y0-y_limit ,y0+y_limit, No_dev_W );


 for ii=1:length(x_start)

 x_sub=x_start(ii)+ linspace (0,L_prim,No_dev_L).* cosd(sui);
 y_sub=y_start(ii)+linspace (0,L_prim,No_dev_L).* sind(sui);
 z_sub=x_sub*0;
 
  points(1,(1:No_dev_L)+ (ii-1)*No_dev_L)=x_sub;
 points(2,(1:No_dev_L)+ (ii-1)*No_dev_L)=y_sub;
 points(3,(1:No_dev_L)+ (ii-1)*No_dev_L)=z_sub;
 
 end
 
 

 

counter=0;
Boarders=zeros(3,No_dev_L*2);


% first row and last row in good order ?!
Boarders=points(1:3,[1:No_dev_L,end:-1: (No_dev_L*(No_dev_W-1)+1)]);

% plot3(Boarders(1,:),Boarders(2,:),Boarders(3,:),'g*')
for ii=1:No_dev_L
points_New(:,(1:No_dev_W)+ (ii-1)*No_dev_W)=points(:,((1:No_dev_W)-1)*No_dev_L + ii);

end


points_New=points_New(:,No_dev_W*starting_index+1:end);

% end

