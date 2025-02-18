function R=fcn_RotationFromTwoVectors(v1, v2)

%% Refrence
% @MISC {1119301,
%     TITLE = {Calculate Rotation Matrix to align Vector A to Vector B in 3d?},
%     AUTHOR = {Shiyu (https://math.stackexchange.com/users/7156/shiyu)},
%     HOWPUBLISHED = {Mathematics Stack Exchange},
%     NOTE = {URL:https://math.stackexchange.com/q/1119301 (version: 2015-05-08)},
%     EPRINT = {https://math.stackexchange.com/q/1119301},
%     URL = {https://math.stackexchange.com/q/1119301}
% }
%%
% R*v1=v2
% v1 and v2 should be column vectors and 3x1

% 1. rotation vector
w=cross(v1,v2);
if norm(w) ~=0
w=w/norm(w);
end

w_hat=fcn_GetSkew(w);
% 2. rotation angle
cos_tht=v1'*v2/norm(v1)/norm(v2);
tht=acos(cos_tht);
% 3. rotation matrix, using Rodrigues' formula
R=eye(size(v1,1))+w_hat*sin(tht)+w_hat^2*(1-cos(tht));

function x_skew=fcn_GetSkew(x)
x_skew=[0 -x(3) x(2);
 x(3) 0 -x(1);
 -x(2) x(1) 0];