
% OTO  > Optical Thermal Optimization of LATW

function Nip_point_Temp_T=OTO_Luncher_Opt(xin,ID_parameter)

%% Help section 
% xin >> is the parameter to be optimized
% ID  >> decides which parameter should be optimized from the list box 

%%

% clc;
% clear;

%% Geometrical
 th_y=-180;
 W_tape = 6; % width of the tape
thick_T=0.4e-3; % thickness of the tape
thick_sub=2*thick_T;
R_tape=5;   % without thickness, it will be added in general 3D tape cylinder
L_flat=R_tape*3;
deg_tape=360/8;
sui=0 ; %pi/8;  % tape_direction in radian, effect of rotation and sui
L_prim=17;   % tape long
% th_1=30; % starting angle in degree

w=8; % width of the tape
z1=21; % position of the nip-point z
R_cyl= 15; %80; %input('The radius of cylinder =');
z_cyl_end= 50; %input('Enetr the end point of cylinder =');
tv3=z_cyl_end/2;  % for the tape
W_R=1.1*W_tape;  % width of the Roller
Rxyz=[-0.8239 -0.1972 -0.3455];
% Laser_head_Ax=2;
% Laser_head_Ay=4;
% Laser_head_nx=10;
% Laser_head_ny=8;
Laser_head=[2 4 10 8];
L_xyz0=[50 27 42];
H_indentation=0;

Laser_head_Rot=0;

%% Process parameters 
materials_Tape=[0.72e-3;1560e-9;1425];
materials_sub=[0.72e-3;1560e-9;1425];
Velocity=0.5 ;
Total_energy=1e2;
ID=2;   % laser distribution pattern 
absorbtion=0.8;
Temp_Right_sub=20;
Temp_Right_T=20;

%% computational parameters
step_size_angle=0.04;  % in radian
node_num_Z=8;
node_space=12;  % change to the number of devision
Angular_space=3;
L_flat_space=3*node_space;  % change to the number of devision

%% DEsicion section

switch  ID_parameter
    case   1
        th_y=xin;
    case 2
        Velocity=xin;
    case 3
        Total_energy=xin;
 
end



%%

Node_z_3D_thermal=4;

%%  >>>>>>>>>  compute section   <<<<<<<<<<
%%  >>>>>>>>>  compute section   <<<<<<<<<<


[Nip_point_Temp_T]=thermal_domain_generator_opt(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,step_size_angle,w,thick_sub,node_num_Z,R_cyl,z_cyl_end,...
    tv3,W_R,materials_Tape,materials_sub,Velocity,Total_energy,...
     Rxyz,ID,Laser_head,L_xyz0,absorbtion,...
     node_space,Angular_space,L_flat_space,...
     Temp_Right_sub,Temp_Right_T,H_indentation,Node_z_3D_thermal,Laser_head_Rot);