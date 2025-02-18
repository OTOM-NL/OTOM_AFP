
% Temp increase of gas in the Tank 40 L pressure vessel
% estimation 
function Temp_roller_increase=Roller_temp_estimation(ro,r_roller,w_roller,cp,time,Q)

%%  Should be modified 

% ro=ro;
V=(r_roller^2)*w_roller*pi;
m=ro*V;  %
% cp=1*1e3;

% h=20;
% A=0.06 *0.025;
% delta_t=45-20;
% time=130;


% Q=h*A*delta_t*time;
% Q=265; % joule per second
Q=Q*time;
Temp_roller_increase=Q/(m*cp);


