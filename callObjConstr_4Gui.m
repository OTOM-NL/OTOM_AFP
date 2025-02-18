



function [x, Fval,exitFlag,Output] = callObjConstr_4Gui(handles)

% Lb=[-180];
% Ub=[-110];
% handles  >> is a struc data
%% Inputs
%% return struc data into numeric and Run


%% Geometrical Parameters
th_y=str2double(handles.Geometrical_parameters{1});
W_tape = str2double(handles.Geometrical_parameters{2}); % width of the tape
thick_T=str2double(handles.Geometrical_parameters{3}); % thickness of the tape
R_tape=str2double(handles.Geometrical_parameters{4});   % without thickness, it will be added in general 3D tape cylinder
L_flat=str2double(handles.Geometrical_parameters{5});
deg_tape=str2double(handles.Geometrical_parameters{6});
sui=str2double(handles.Geometrical_parameters{7}) ; %pi/8;  % tape_direction in radian, effect of rotation and sui
L_prim=str2double(handles.Geometrical_parameters{8});   % tape long
% th_1=str2double(handles.Geometrical_parameters{9}); % starting angle in degree
w=str2double(handles.Geometrical_parameters{9}); % width of the substrate tape
thick_sub=str2double(handles.Geometrical_parameters{10}); % thiness of substrate
R_cyl= str2num(handles.Geometrical_parameters{11}); %80; %input('The radius of cylinder =');
z_cyl_end= str2double(handles.Geometrical_parameters{12}); %input('Enetr the end point of cylinder =');
Roller_Pos_TV=str2num(handles.Geometrical_parameters{13});  % for the tape
W_R=str2double(handles.Geometrical_parameters{14});  % width of the Roller
Rxyz=str2num(handles.Geometrical_parameters{15});
Laser_head=str2num(handles.Geometrical_parameters{16});
L_xyz0=str2num(handles.Geometrical_parameters{17});
   Laser_head_Rot=str2double(handles.Geometrical_parameters{18});
    


%% Process Parameters
materials_Tape=str2num(handles.Process_parameters{1});
materials_sub=str2num(handles.Process_parameters{2});
Velocity=str2double(handles.Process_parameters{3});
Total_energy=str2double(handles.Process_parameters{4});
ID=str2num(handles.Process_parameters{5});   % laser distribution pattern
absorbtion_waste=str2num(handles.Process_parameters{6});
Temp_Right_sub=str2num(handles.Process_parameters{7}); % Includes also mandrel Temperature
Temp_Right_T=str2num(handles.Process_parameters{8}); 
h_conv_Sub=str2double(handles.Process_parameters{9});
h_conv_T=str2num(handles.Process_parameters{10});
Roller_Force=str2double(handles.Process_parameters{11});


assignin('base','Consolidation_Force',Roller_Force);

if Roller_Force ~=0
    
    if isempty(handles.fit_Func)
        
        uiwait(warndlg('No Force-displacement data, Please load the data.'));
        H_indentation=0;
        %         Roller_def_Callback(hObject, eventdata, handles);
    else
        
        fitresult= handles.fit_Func;
        H_indentation=fitresult(Roller_Force);
        figure(50);
        plot(Roller_Force,H_indentation,'ko');
        text (Roller_Force,H_indentation,' Maximum normal deformation');
        javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
        
    end
    
else
    H_indentation=0;
    
end

if H_indentation > R_tape
    h=errordlg('Normal displacement is more than Roller radius!');
    H_indentation=R_tape;
end



%% Computational parameters

step_size_angle=str2double(handles.Computational_parameters{1});  % in radian
No_dev=str2double(handles.Computational_parameters{2});
node_space=str2double(handles.Computational_parameters{3});  % change to the number of devision
Angular_space=str2double(handles.Computational_parameters{4});
L_flat_space=str2double(handles.Computational_parameters{5});  % change to the number of devision
  % change to the number of devision
%%

%% Temperature of the mandrel
if length (Temp_Right_sub) ==2
    T_amb_mandrel=Temp_Right_sub(2);
else
T_amb_mandrel=30;
end
%%  Check graphical outputs
if isempty(handles.checkbox1)
    Graphic_chekbox=ones(1,9);
else
Graphic_chekbox=[handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,...
    handles.checkbox5,handles.checkbox6,handles.checkbox7,handles.checkbox8,...
    handles.checkbox9];
end

Node_z_3D_thermal=handles.Node_z_3D_thermal;
dos('animator_19_frame_camtes.exe -i  &');
%%

Temp_desired_tape=str2double(get(handles.OP_edit2,'String'));   % desired temperature of the Tape
Temp_desired_sub=str2double(get(handles.OP_edit1,'String'));  % desired temperature of the substrate

ID_parameter=get(handles.OP_popupmenu1,'Value') ;  % show which parameter will be optimized

Exact_Temp_checkBox=get(handles.OP_checkbox1,'Value') ;  
Uniform_Temp_checkBox=get(handles.OP_checkbox2,'Value') ;  

%%

% use or not use divergence of the laser
Laser_div_OnOff=get(handles.checkbox_laser_divergence,'Value');

if Laser_div_OnOff

 fid23 = fopen('.\Supp_files\Laser_characteristic.txt','r');
 out = textscan(fid23,'%s ','delimiter',',');
 fclose(fid23);
 
Divergence_factor=str2double(out{1}{2});

else
    Divergence_factor=0;
    
end

%%
% poolobj = gcp('nocreate');
% 
% if isempty(poolobj)
% %     poolsize = 0;
%     parpool('local',2);
% else
% %     poolsize = poolobj.NumWorkers;
% end






oldData=[];

%%
Lb=str2num(handles.Range{1}); % for multiple input variable, it is a combination of all variable
Ub=str2num(handles.Range{2});


% [x, fval]=ga(@myObjective,1,[] ,[], [], [],Lb,Ub,@ellipseparabola);


A = []; b = []; % No linear inequality constraints
Aeq = []; beq = []; % No linear equality constraints
options = optimoptions('ga','FunctionTolerance',1e0,'ConstraintTolerance',1e-4,'MaxStallGenerations',5,...
'MaxStallTime',3,'UseParallel','never');

% 'PlotFcns',@gaplotpareto,'UseParallel','always'

nvar=length(Lb);

[x,Fval,exitFlag,Output]=gamultiobj(@myObjective,nvar,A, ...
    b,Aeq,beq,Lb,Ub,options);
fprintf('The number of points on the Pareto front was: %d\n', size(x,1));
fprintf('The number of generations was : %d\n', Output.generations);
figure(20);
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

plot3(x,Fval(:,1),Fval(:,2),'r*');
s=sprintf('Number of design= %d, Pareto design= %d',Output.funccount,size(x,1));
title(['Multi-obj-optimization, ',s]);
xlabel('Var');
ylabel('Obj1 (mean difference)');
zlabel('Obj2 (STD)');
grid on;



% [x fval] = fmincon(@myObjective,x0,[],[],[],[],[],[], ...
%     @ellipseparabola,options);

function f = myObjective(xin)
% Temp=thermal_domain_generator_opt(xin);
% [Nip_point_Temp_T]=thermal_domain_generator_opt

% >>>> where is xin

switch  ID_parameter
    case   1
        th_y=xin;
    case 2
        Velocity=xin;
    case 3
        Total_energy=xin;
        case 4
        L_xyz0=xin;
        
  case 5
                %Laser_direction (Rx,Ry,Rz)
     Rxyz  =xin;

    case 6
        % Laser Head Size (Ax,Ay)
      Laser_head(1:2)  =xin;
         
      case 7
                %Laser-ID parameter  (Par_x, Par_Y)
        ID (2:3)  =xin;
      
    case 8
        % Mandrel radius R_cyl
     R_cyl (1)   =xin;
    
        
end
disp(xin);
 

 
 [Temp_sub_nip,Temp_tape_nip ]=thermal_domain_generator_opt(th_y,W_tape,thick_T,R_tape,L_flat,deg_tape,...
    sui,L_prim,step_size_angle,w,thick_sub,No_dev,R_cyl,z_cyl_end,...
    Roller_Pos_TV,W_R,materials_Tape,materials_sub,Velocity,Total_energy,...
    Rxyz,ID,Laser_head,L_xyz0,absorbtion_waste,...
    node_space,Angular_space,L_flat_space,...
    Temp_Right_sub,Temp_Right_T,h_conv_Sub,h_conv_T,H_indentation,Node_z_3D_thermal,T_amb_mandrel,...
             Graphic_chekbox,handles.Input_reader_mode,Laser_head_Rot,...
             handles.BRDF_mode,Divergence_factor);
 
 
 
 
%  Temp_sub_nip,Temp_tape_nip 
 
f=zeros(1,2);
counter=0;
 if Exact_Temp_checkBox
counter=counter+1;
% m=length(Temp);
     f(counter)=mean(abs(Temp_tape_nip-Temp_desired_tape));
     counter=counter+1;
       f(counter)=mean(abs(Temp_sub_nip-Temp_desired_sub));
 end
 
  if Uniform_Temp_checkBox
counter=counter+1;
     f(counter)=mean ([std(Temp_sub_nip),std(Temp_tape_nip)]);
 end
 
% f=[norm(abs(Temp-Temp_desired)), std(Temp)];



newRow=[xin, f(:)'];
newData = [oldData; newRow];
set(handles.OP_uitable1, 'Data', newData);
oldData = get(handles.OP_uitable1,'Data');
disp(f);
h1=figure(1);



% Roller_Pos_TV(1)
% Roller_Pos_TV(2)
% Roller_Pos_TV(3)

 axis([-1*R_cyl(1)+Roller_Pos_TV(1), 1*R_cyl(1)+Roller_Pos_TV(1), -1*R_cyl(1)+Roller_Pos_TV(2), 1*R_cyl(1)+Roller_Pos_TV(2),  -(R_cyl(1))+Roller_Pos_TV(3), (R_cyl(1))+Roller_Pos_TV(3)]);
% axis tight;
% view([150 -60]);
 pause(0.01);

%% This drawnow causes the matlab to stop easily during execution
drawnow;
 axcap = screencapture(h1.CurrentAxes);
set(gca,'CameraViewAngle',2);
% h23=figure(23);
% imshow(h23,axcap);


imshow(axcap, 'Parent', handles.axes1);

 



%  set(h1,'Visible','off');

% f = 3*(xin(1)-xin(2))^4 + 4*(xin(1)+xin(3))^2/(1+sum(xin.^2)) ...
%     + cosh(xin(1)-1) + tanh(xin(2)+xin(3));
end

function [c,ceq] = ellipseparabola(x)
% c(1) = (x(1)^2)/9 + (x(2)^2)/4 - 1;
% c(2) = x(1)^2 - x(2) - 1;
c=[];
ceq = [];

end
end