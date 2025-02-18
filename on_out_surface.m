



function theta_all=on_out_surface (p,Ps)
%% This function compare specific point whether it is in the boarder of the points or Not in xz plane, 
% the sum of the angle for each two point form the boarder should be roughly 360


 [m,n]=size(Ps);
if n<2
Ps=[Ps Ps];
end


%%

% p=[8;31.45]; % x,z
theta_all=0;

% plot3(p(1),10,p(2),'c*');

for ii=1:length(Ps)

a=Ps(1:2,ii);

if mod(ii,length(Ps)) ==0
    ii=0;
end


b=Ps(1:2,ii+1);


line1=a-p;
line2=b-p;

theta=acosd(dot(line1,line2)/(norm(line1)*norm(line2)));
theta_all=theta_all+abs(theta);


end

theta_all=real(theta_all); 

%  plot(Ps(1,:),Ps(2,:),'r.')
% plot(p(1),p(2),'c*')
