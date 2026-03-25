
%normal and reflection for Tape -Ray 3D
function [kr1_active_G,beta_macro]=tape_ray_reflection_3D_BRDF (R,S1,S2,xyz_int_L,x0,y0,laser_direction,tv,Rot_Roller_axis,fileID_Normals,fileID_Refs,xyz_int_G,BRDF_Tape)


% BRDF_Tape=[Fib_Rot,div,sig_T,sig_F,Amp,threshold]


% Fib_Rot= 90; % in degree
