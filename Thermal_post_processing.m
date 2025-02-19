



function Thermal_post_processing(jobname,R_cyl,z_cyl_end,nip_point_M_all,Rot_Roller_axis_all,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,H_indentation,tv_def_all,CV_mesh,UOT_pathfile,...
    Laser_head,Points_in_domain_all_sub,Points_in_domain_all_Tape,Total_energy,Stlname,Scale_factor, ...
    Laser_modul_No,varargin)

ss=1;

warning('off');

Tape_Rack_steps=varargin{1};

%     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'r');   % file includes the xyz + ID
%     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'r');     % file includes the xyz + ID
%     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'r');
%
%
%     %%  instead of plotting in each loop, plot all afterwards
%     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'r');
%     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'r');
%     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'r');
%     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'r');
%     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'r');
%     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'r');
%
%     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'r');
%



str_name=split(Stlname,"/");


fig_HO=figure('Name',strcat('UOT Analysis <>',str_name{end}),'NumberTitle','off');
mouse3d;
% % fig_HO=figure('Name','Optical UOT Analysis','NumberTitle','off');
cameratoolbar(fig_HO, 'Show');
javaFrame    = get(gcf,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

set(gcf,'color','w');
xlabel('x');
ylabel('y');
zlabel('z');

% create structure of handles
% Fig_handles = guihandles(fig_HO);
% % Add some additional data as a new field called numberOfErrors
% % Fig_handles.numberOfErrors = 0;
% % Save the structure
% guidata(fig_HO,Fig_handles)




Fig_handles = guidata(gcbo);


% plot Mandrel
% First is cylinder , second is bottom dome part, third is upper Dome part
% Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);


if isempty(Stlname)
    Mandrel_plot=plot3D_cylinder(R_cyl,z_cyl_end);

else
 [~,~,~,Mandrel_plot]=stldemo(Stlname);   
 
 
 Mandrel_plot.Vertices=Mandrel_plot.Vertices*Scale_factor;
 
 

end


     set(Mandrel_plot,...
    'FaceLighting','phong',...
'AmbientStrength',.4,'DiffuseStrength',.8,...
'SpecularStrength',.2,'SpecularExponent',25, ...
'SpecularColorReflectance',0.5,'BackFaceLighting','reverselit');

     camlight('headlight');

%      camlight('left');
% camlight('right');
% lightangle(-35,50);
% material(Mandrel_plot,'shiny');


assignin('base','Mandrel_plot',Mandrel_plot);

if ~isempty(Tape_Rack_steps)

  Tape_Rack=copyobj(Tape_Rack_steps.TapeRack_patch,fig_HO.CurrentAxes);
%     Ver_stl_TRack=Tape_Rack_steps.Vertices{ss};
%      Face_stl_TapeRack=Tape_Rack_steps.Faces;

else
    Tape_Rack=[];

end



set(fig_HO,'Visible','on');

% plot3(tv(1),tv(2),tv(3),'g*','Markersize',5)


%    Roller_Pos_TV=nip_point_M_all(ss,:,:);


%     tv(1)=Roller_Pos_TV(1)-(-tv_def(1));
%     tv(2)=Roller_Pos_TV(2)-tv_def(2);
%     tv(3)=Roller_Pos_TV(3)+tv_def(3);
%
%      Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
% output is real center of roller without deformation, end of arrow of
% the local axis
%   [n1,n2,n3,x0,y0,z0,Fig_handles.h_R,Fig_handles.h_T,Center_roller]=Tape_plot_3D_Tshow(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);

% assignin('base','h_R',Fig_handles.h_R);
% assignin('base','h_T',Fig_handles.h_T);




%     XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
%     XYZ_int=cell2mat(XYZ_int);
%     Fig_handles.Laser_int_points=plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
%
%
%
%     XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     XYZ_Rays=cell2mat(XYZ_Rays);
%    Fig_handles.Graphics_laser= plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',2 );
%
%
%
%     Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
%     Fig_handles.Normal_lines1=plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:');
%
%
%
%
%     Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
%     Fig_handles.Graphics_reflection1=plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--');
%




%     Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
%    Fig_handles.Normal_lines2= plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
%
%
%
%
%     Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
%     Fig_handles.Graphics_reflection2=plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--');
%
%
%
%     Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
%     Ref_int_xyz=cell2mat(Ref_int_xyz);
%     Fig_handles.reflection_int=plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');


%% Changing the graphical object, Visibility or deleting

%     Graphics_laser=findobj(h,'color','g','LineStyle',':');
%     assignin('base','Graphics_laser',Fig_handles.Graphics_laser);
%     % set(Graphics_laser,'Visible','off');
% %     Graphics_reflection1=findobj(h,'color','r');
% %     Graphics_reflection2=findobj(h,'color','b');
%     Fig_handles.Graphics_reflection=[Fig_handles.Graphics_reflection1;Fig_handles.Graphics_reflection2];
%     assignin('base','Graphics_reflection',Fig_handles.Graphics_reflection);
%
% %     reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
%     assignin('base','reflection_int',Fig_handles.reflection_int);
%
% %     Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
%     assignin('base','Laser_int_points',Fig_handles.Laser_int_points);
%
% %     Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
% Fig_handles.Normal_lines=[Fig_handles.Normal_lines1;Fig_handles.Normal_lines2];
%     assignin('base','Normal_lines',Fig_handles.Normal_lines);

% Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
% assignin('base','Thermal_points',Thermal_points);




% set(Graphics_reflection,'Visible','off');

Laser_head_Ax=Laser_head(1);
Laser_head_Ay=Laser_head(2);

Laser_head_nx=Laser_head(3);
Laser_head_ny=Laser_head(4);
ID=[0 0 0];

Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

%% Create a push button

Fig_handles.Graphics_laser=[];
Fig_handles.Graphics_reflection=[];
Fig_handles.h_R=[];
Fig_handles.h_T=[];
Fig_handles.reflection_int=[];
Fig_handles.Laser_int_points=[];
Fig_handles.Normal_lines=[];
Fig_handles.Sub=[];
Fig_handles.P_sub_in=[];




guidata(gcbo,Fig_handles);

rb1 = uicontrol('Style', 'radiobutton','Position',[20 50 90 40],'Value',true,...
    'String','Rays','Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Fig_handles.Graphics_laser),'Backgroundcolor','w');

btn1 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser rays',...
    'Position', [20 20 90 20],...
    'Callback', 'delete(Graphics_laser)','Backgroundcolor','w');


rb2 = uicontrol('Style', 'radiobutton','Position',[150 50 90 40],'Value',true,...
    'String','Reflections','Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Fig_handles.Graphics_reflection),'Backgroundcolor','w');

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Clear Reflections',...
    'Position', [150 20 90 20],...
    'Callback', 'delete(Graphics_reflection)','Backgroundcolor','w');



rb3 = uicontrol('Style', 'radiobutton','Position',[280 50 90 40],'Value',true,...
    'String','Mandrel','Callback', @(rb3,event) Visible_plotButtonPushed(rb3,Mandrel_plot),'Backgroundcolor','w');

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Clear Mandrel',...
    'Position', [280 20 90 20],...
    'Callback', 'delete(Mandrel_plot)','Backgroundcolor','w');


rb4 = uicontrol('Style', 'radiobutton','Position',[400 50 90 40],'Value',true,...
    'String','Roller','Callback', @(rb4,event) Visible_plotButtonPushed(rb4,Fig_handles.h_R),'Backgroundcolor','w');

btn4 = uicontrol('Style', 'pushbutton', 'String', 'Clear Roller',...
    'Position', [400 20 90 20],...
    'Callback', 'delete(h_R)','Backgroundcolor','w');


rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,...
    'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]),'Backgroundcolor','w');

btn5 = uicontrol('Style', 'pushbutton', 'String', 'Clear Tape-Sub',...
    'Position', [520 20 90 20],...
    'Callback', 'delete(Tape_sub)','Backgroundcolor','w');




rb6 = uicontrol('Style', 'radiobutton','Position',[640 50 90 40],'Value',true,...
    'String','1st reflections','Callback', @(rb6,event) Visible_plotButtonPushed(rb6,Fig_handles.reflection_int),'Backgroundcolor','w');
btn6 = uicontrol('Style', 'pushbutton', 'String', 'Clear reflection_int',...
    'Position', [640 20 90 20],...
    'Callback', 'delete(reflection_int)','Backgroundcolor','w');

%%
rb10 = uicontrol('Style', 'radiobutton','Position',[1000 50 90 40],'Value',true,...
    'String','Point_Sub inside','Callback', @(rb10,event) Visible_plotButtonPushed(rb10,Fig_handles.P_sub_in),'Backgroundcolor','w');
btn10 = uicontrol('Style', 'pushbutton', 'String', 'Clear P_sub_in',...
    'Position', [1000 20 90 20],...
    'Callback', 'delete(P_sub_in)','Backgroundcolor','w');





rb7 = uicontrol('Style', 'radiobutton','Position',[760 50 90 40],'Value',true,...
    'String','Laser_int','Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Fig_handles.Laser_int_points),'Backgroundcolor','w');

btn7 = uicontrol('Style', 'pushbutton', 'String', 'Clear Laser_int_points',...
    'Position', [760 20 90 20],...
    'Callback', 'delete(Laser_int_points)','Backgroundcolor','w');


%
% rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
%     'String','Thermal_P','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Thermal_points));

rb8 = uicontrol('Style', 'radiobutton','Position',[900 50 90 40],'Value',true,...
    'String','Normal lines','Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Fig_handles.Normal_lines),'Backgroundcolor','w');


%     string=[]
btn8 = uicontrol('Style', 'pushbutton', 'String', 'Clear Normal_lines',...
    'Position', [900 20 90 20],...
    'Callback', 'delete(Normal_lines)','Backgroundcolor','w');




rb9 = uicontrol('Style', 'radiobutton','Position',[1040 50 90 40],'Value',true,...
    'String','Tape Rack','Callback', @(rb9,event) Visible_plotButtonPushed(rb9,Tape_Rack),'Backgroundcolor','w');

if ~isempty(Tape_Rack_steps)
  rb9.Enable='off';
end

text_orientation = uicontrol('Style','text','Units','normalized',...
    'String',sprintf('Total Absorbed Energy: \n Tape=%f \n Mandrel=%f \n Roller=%f \n Roller=%f',0 ,0,0,0),'FontUnit','normalized','Fontsize',0.10,...
    'Position',[0.75 0.75 0.25 0.25],'Visible','on','Backgroundcolor','w','HorizontalAlignment','right');



% text_orientation = annotation('textbox','Units','normalized',...
%     'String',sprintf('Total Absorbed Energy: \n Tape=%f \n Mandrel=%f \n Roller=%f',0 ,0,0),'Fontsize',0.15,'FontUnit','normalized',...
%     'Position',[0.8 0.8 0.2 0.15],'Visible','on','Backgroundcolor','w','HorizontalAlignment','right','FaceAlpha',0.5,'EdgeColor','none');



% jFigPeer = get(handle(gcf),'JavaFrame');
% jWindow = jFigPeer.fFigureClient.getWindow;
% com.sun.awt.AWTUtilities.setWindowOpacity(jWindow,0.7)

sld_transparency= uicontrol('Style', 'slider',...
    'Min',0,'Max',1,'Value',1,...
    'Units','normalized'...
    ,'Position', [0.09 0.98 0.08 0.02],'Callback',  @(sld_transparency,event) change_transp(sld_transparency));


pop_items = uicontrol('Style', 'popupmenu','String', {1:length(nip_point_M_all)},...
    'Units','normalized'...
    ,'Position', [0.01 0.98 0.08 0.02],...
    'Callback',  @(pop_items,event) show_step_outputs(text_orientation,pop_items,jobname,nip_point_M_all,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def_all,...
    rb1,rb2,rb3,rb4,rb6,rb7,rb8,rb10,...
    sld_transparency,Tape_Rack,Tape_Rack_steps));


%     'Units','normalized'...
%     ,'Position', [0.09 0.98 0.08 0.02],

Animation_items = uicontrol('Style', 'pushbutton','String', 'Animation',...
         'Units','normalized'...
     ,'Position', [0.01 0.95 0.08 0.02],...
    'Callback',  @(Animation_items,event) animation_creation);

% 'Position', [1400 20 140 20 ],...

text_box_ID = uicontrol('Style', 'edit','String', '0 0 0',...
    'Units','normalized'...
    ,'Position', [0.01 0.92 0.08 0.02]    );


btn_save_absorbed_data = uicontrol('Style', 'pushbutton', 'String', 'Save absorbed energy',...
    'Position', [1200 20 140 20],...
    'Callback',  @(btn_save_absorbed_data,event) save_absorbed_data,'Backgroundcolor','r');


    btn_optimized_laser = uicontrol('Style', 'radiobutton', 'String', 'Use optimized Laser',...
        'Position', [1200 50 140 20],...
        'Backgroundcolor','w');

% text_orientation = uicontrol('Style','pushbutton','Units','normalized',...
%     'String','Save absorbed energy','Fontsize',10,...
%     'Position',[0.8 0.6 0.2 0.15],'Visible','on','Backgroundcolor','w','HorizontalAlignment','right');

% 'Callback', @surfzlim


[cdata,map] = imread('complete.jpg');

h2=msgbox('Please select a step number from Pop-up menu at top left corner!',...
    'Notice','custom',cdata,map);
javaFrame    = get(h2,'JavaFrame');
iconFilePath = 'OTOM-icon.png';
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
set(h2,'color','w');



function animation_creation
    
Video = VideoWriter(strcat(UOT_pathfile,'Kinematics_Optical'),'MPEG-4');
% Video = VideoWriter('.\Kinematics_Optical','MPEG-4');
% video.FrameRate = 10;
Video.FrameRate=20; Video.Quality = 99;
open(Video);

h1=gca;


counter=0;
camera_center=get(gca,'CameraPosition');

for nn=1:1:length(nip_point_M_all)
pop_items.Value=nn;



show_step_outputs(text_orientation,pop_items,jobname,nip_point_M_all,...
    N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def_all,...
    rb1,rb2,rb3,rb4,rb6,rb7,rb8,rb10,...
    sld_transparency,Tape_Rack,Tape_Rack_steps)
	drawnow;
% 	pause(0.1);
     Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
     
         Roller_Pos_TV=reshape(nip_point_M_all(ss,:,:),[1,3]);
      
        set(gca,'CameraViewAngle',R_cyl(1)*45);
  camroll(0.05);
         set(gca,'CameraTarget',Roller_Pos_TV);
%           set(gca,'CameraPosition',camera_center+(camera_center*Rot_Roller_axis)*(nn/length(nip_point_M_all)));
           
            
           
%         CameraPosition
%         CameraTarget
%         CameraUpVector
%        figure(59);
% pause(0.1);
counter=counter+1;
Mov(counter)=getframe(gcf);
end
    
writeVideo(Video,Mov);

   fclose all;
end



    function change_transp (sld_transparency)
        data=get(gca,'ch');
        
        for i=1:7
            data(i).Color(4) = get(sld_transparency,'value');
            
        end
        
    end


    function Visible_plotButtonPushed(rb1,Graphics_laser)
        %         x = linspace(0,2*pi,100);
        %         y = sin(x);
        %         plot(ax,x,y)
        if rb1.Value ==false
            set(Graphics_laser,'Visible','off');
        else
            set(Graphics_laser,'Visible','on');
        end
    end


% save absorbed  energy data
    function save_absorbed_data
        
       A= get(  btn_optimized_laser,'Value');
       
       if A
           
           
           
           
             fileID_selected_laser = fopen(strcat(UOT_pathfile, sprintf('Selected_Laser_ID.txt')),'r');
             ID_power= textscan(fileID_selected_laser,' %f %f %f %f %f %f %f %f  ','Delimiter',',','HeaderLines',0) ;
            ID_power=cell2mat(ID_power);
            
              [m,n]= size(ID_power);
              
              
                  prompt = {'Starting step:'};
dlg_title = 'Step of UOT Pareto generation from  optimization ';
 defaultans = {'50'};
 Answer= inputdlg(prompt,dlg_title,[1 50],defaultans);
              
              
        Start_Step= str2num(cell2mat(Answer));
        
       else
        
        ID=str2num( get(text_box_ID,'String'));
        Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

         Power_Actual_all=zeros(Laser_head_ny,Laser_head_nx,Laser_modul_No);
for ii=1:Laser_modul_No
   
    Power_Actual_all(:,:,ii)=Power_Actual;
end



     
        [m,n]= size(tv_def_all);
        Start_Step=1;
        
        
       end
        
        
        
        fileID_absorbed_objects = fopen(strcat(UOT_pathfile, sprintf('Absorbed_power_objects.txt')),'w');
        
        
        fprintf(fileID_absorbed_objects,'Mandrel:    Roller:    Substrate:   Tape: \r\n' );
        
       
        counter=0;
        
        for ss=(1:m)+Start_Step-1
            %         tv_def=tv_def_all(ss,:,:);
            
            counter=counter+1;
              if A
                    ID=ID_power(counter,1:end-1);
             Total_energy=ID_power(counter,end);
        Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

         Power_Actual_all=zeros(Laser_head_ny,Laser_head_nx,Laser_modul_No);
for ii=1:Laser_modul_No
   
    Power_Actual_all(:,:,ii)=Power_Actual;
end


        
              end
        
        
        
      
        
            
            fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'r');
            fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'r');
            fileID3 = fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'r');
            
               %TapeRack
fileID21 = fopen(strcat(UOT_pathfile, sprintf('TapeRack_ints%d.txt',ss)),'r');
            
            fileID4 = fopen(strcat(UOT_pathfile, sprintf('int_xyz%d.txt',ss)),'r');
            fileID5 = fopen(strcat(UOT_pathfile, sprintf('Laser_Rays%d.txt',ss)),'r');
            fileID6 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Mandrel%d.txt',ss)),'r');
            fileID7 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Mandrel%d.txt',ss)),'r');
            fileID8 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Tape%d.txt',ss)),'r');
            fileID9 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Tape%d.txt',ss)),'r');
            fileID11 = fopen(strcat(UOT_pathfile, sprintf('Ref_int_xyz%d.txt',ss)),'r');
            
            
            
            Mandrel_E= textscan(fileID1,' %*f %*f %*f %f %f ','Delimiter',',','HeaderLines',1) ;
            Mandrel_E=cell2mat(Mandrel_E);
            
            Mandrel_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Mandrel_E(:,2)).*Mandrel_E(:,1);
            
            
            Tape_E = textscan(fileID2,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
            Tape_E=cell2mat(Tape_E);
            
            Tape_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Tape_E(:,2)).*Tape_E(:,1);
            
            
            %         Rolelr_energy = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
            %         Rolelr_energy=cell2mat(Rolelr_energy);
            
            Roller_E = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
            Roller_E=cell2mat(Roller_E);
            
            Roller_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Roller_E(:,2)).*Roller_E(:,1);
            
               % Tape Rack
                  TRack_E = textscan(fileID21,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        TRack_E=cell2mat(TRack_E);
        
          TRack_energy=(Total_energy/(Laser_head_nx*Laser_head_ny*Laser_modul_No))*Power_Actual_all(TRack_E(:,2)).*TRack_E(:,1);
          

            %%
            P_sub_in=Points_in_domain_all_sub{ss};
            P_Tape_in=Points_in_domain_all_Tape{ss};
            
            
            Absorbed_sub= (Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(P_sub_in(:,5)).*P_sub_in(:,4);
            Absorbed_Tape= (Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(P_Tape_in(:,5)).* P_Tape_in(:,4);
            
            
            fprintf(fileID_absorbed_objects,'  %f %f  %f %f %f \r\n',sum(Mandrel_energy),sum(Roller_energy),...
                sum(Absorbed_sub), sum(Absorbed_Tape) ,sum(TRack_energy) );
            
          
            
            
            %         set (text_orientation ,'String',sprintf('Total Absorbed Energy: \n Tape=%f \n Mandrel=%f \n Roller=%f \n Substrate=%f \n Tape=%f',sum(Tape_energy),sum(Mandrel_energy),sum(Roller_energy),...
            %              sum(Absorbed_sub), sum(Absorbed_Tape)));
            fclose(fileID1);
            fclose(fileID2);
            fclose(fileID3);
            fclose(fileID4);
            fclose(fileID5);
            fclose(fileID6);
            fclose(fileID7);
            fclose(fileID8);
            fclose(fileID9);
            fclose(fileID11);
            fclose(fileID21);
            
        end
        
        if A
            fclose(fileID_selected_laser);
            
        end
            
            
        fclose(fileID_absorbed_objects);
        
        
        
    end

    function show_step_outputs(text_orientation,src,jobname,nip_point_M_all,...
            N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis_all,H_indentation,tv_def,...
            rb1,rb2,rb3,rb4,rb6,rb7,rb8,rb10,sld_transparency,Tape_Rack,Tape_Rack_steps)
        val = src.Value;
        str = src.String;
        %         str{val};
        %         disp(['Selection: ' str{val}]);
        
        % Fig_handles = guidata(gcbo);
        
        ID=str2num( get(text_box_ID,'String'));
        Power_Actual=Laser_Power_generator (Laser_head_Ax,Laser_head_Ay,Laser_head_nx,Laser_head_ny,ID);

        Power_Actual_all=zeros(Laser_head_ny,Laser_head_nx,Laser_modul_No);
for ii=1:Laser_modul_No
   
    Power_Actual_all(:,:,ii)=Power_Actual;
end


        
        ss=val;
        tv_def=tv_def_all(ss,:,:);
        
        % jobname='Example0';
        %     fileID1 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Cylinder_ints%d.txt'),jobname,ss),'r');   % file includes the xyz + ID
        %     fileID2 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Tape_ints%d.txt'),jobname,ss),'r');     % file includes the xyz + ID
        %     fileID3 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Roller_ints%d.txt'),jobname,ss),'r');
        
        fileID1 = fopen(strcat(UOT_pathfile, sprintf('Cylinder_ints%d.txt',ss)),'r');
        fileID2 = fopen(strcat(UOT_pathfile, sprintf('Tape_ints%d.txt',ss)),'r');
        fileID3 = fopen(strcat(UOT_pathfile, sprintf('Roller_ints%d.txt',ss)),'r');
        
          %TapeRack
fileID21 = fopen(strcat(UOT_pathfile, sprintf('TapeRack_ints%d.txt',ss)),'r');
        
        %%  instead of plotting in each loop, plot all afterwards
        %     fileID4 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\int_xyz%d.txt'),jobname,ss),'r');
        %     fileID5 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Laser_Rays%d.txt'),jobname,ss),'r');
        %     fileID6 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Mandrel%d.txt'),jobname,ss),'r');
        %     fileID7 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Mandrel%d.txt'),jobname,ss),'r');
        %     fileID8 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Normal_vectors_Tape%d.txt'),jobname,ss),'r');
        %     fileID9 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Reflection_vector_Tape%d.txt'),jobname,ss),'r');
        %
        %     fileID11 = fopen(sprintf(string('.\\Analysis_UOT\\%s\\Ref_int_xyz%d.txt'),jobname,ss),'r');
        fileID4 = fopen(strcat(UOT_pathfile, sprintf('int_xyz%d.txt',ss)),'r');
        fileID5 = fopen(strcat(UOT_pathfile, sprintf('Laser_Rays%d.txt',ss)),'r');
        fileID6 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Mandrel%d.txt',ss)),'r');
        fileID7 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Mandrel%d.txt',ss)),'r');
        fileID8 = fopen(strcat(UOT_pathfile, sprintf('Normal_vectors_Tape%d.txt',ss)),'r');
        fileID9 = fopen(strcat(UOT_pathfile, sprintf('Reflection_vector_Tape%d.txt',ss)),'r');
        fileID11 = fopen(strcat(UOT_pathfile, sprintf('Ref_int_xyz%d.txt',ss)),'r');
        
        
        
        Mandrel_E= textscan(fileID1,' %*f %*f %*f %f %f ','Delimiter',',','HeaderLines',1) ;
        Mandrel_E=cell2mat(Mandrel_E);
        
        Mandrel_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Mandrel_E(:,2)).*Mandrel_E(:,1);
        
        
        Tape_E = textscan(fileID2,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        Tape_E=cell2mat(Tape_E);
        
        Tape_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Tape_E(:,2)).*Tape_E(:,1);
        
        
        %         Rolelr_energy = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        %         Rolelr_energy=cell2mat(Rolelr_energy);
        
        Roller_E = textscan(fileID3,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        Roller_E=cell2mat(Roller_E);
        
        Roller_energy=(Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(Roller_E(:,2)).*Roller_E(:,1);
        

          % Tape Rack
                  TRack_E = textscan(fileID21,' %*f %*f %*f %f %f','Delimiter',',','HeaderLines',1) ;
        TRack_E=cell2mat(TRack_E);
        
          TRack_energy=(Total_energy/(Laser_head_nx*Laser_head_ny*Laser_modul_No))*Power_Actual_all(TRack_E(:,2)).*TRack_E(:,1);
          
        %%
        P_sub_in=Points_in_domain_all_sub{ss};
        P_Tape_in=Points_in_domain_all_Tape{ss};
        
        
        Absorbed_sub= (Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(P_sub_in(:,5)).*P_sub_in(:,4);
        Absorbed_Tape= (Total_energy/(Laser_head_nx*Laser_head_ny))*Power_Actual_all(P_Tape_in(:,5)).* P_Tape_in(:,4);
        
        
        
        
        set (text_orientation ,'String',sprintf('Total Absorbed Energy: \n Tape=%6.2f \n Mandrel=%6.2f \n Roller=%6.2f \n Substrate=%6.2f\n Tape=%6.2f \n TRack=%6.2f', ...
            sum(Tape_energy),sum(Mandrel_energy),sum(Roller_energy),...
            sum(Absorbed_sub), sum(Absorbed_Tape),sum(TRack_energy)));
        
        
          % update location and Vertices of Tape Rack
        if ~isempty(Tape_Rack_steps)
        Tape_Rack.Vertices=Tape_Rack_steps.Vertices{ss};
        end
        
        
        % h=figure(1);
        
        
        
        % delete([h_R])
        delete([Fig_handles.Sub,Fig_handles.h_R',Fig_handles.h_T,Fig_handles.Graphics_laser,Fig_handles.Graphics_reflection',Fig_handles.reflection_int,Fig_handles.Laser_int_points,Fig_handles.Normal_lines',...
            Fig_handles.P_sub_in]);
        % plot Mandrel
        % First is cylinder , second is bottom dome part, third is upper Dome part
        
        
        
        Roller_Pos_TV=nip_point_M_all(ss,:,:);
        
        
        tv(1)=Roller_Pos_TV(1)+(-tv_def(1));
        tv(2)=Roller_Pos_TV(2)-tv_def(2);
        tv(3)=Roller_Pos_TV(3)-tv_def(3);
        
        Rot_Roller_axis= reshape(Rot_Roller_axis_all(ss,:,:),[3 3]);
        % output is real center of roller without deformation, end of arrow of
        % the local axis
        [n1,n2,n3,x0,y0,z0,Fig_handles.h_R,Fig_handles.h_T,Center_roller]=Tape_plot_3D_Tshow(N_tape,W_tape,R_tape,L_flat,tv,th_y,thick_T,deg_tape,W_R,theta_ind,Rot_Roller_axis,H_indentation,tv_def);
        
        % assignin('base','h_R',Fig_handles.h_R);
        % assignin('base','h_T',Fig_handles.h_T);
        

        set(Fig_handles.h_R,...
    'FaceLighting','phong',...
'AmbientStrength',.4,'DiffuseStrength',.8,...
'SpecularStrength',.2,'SpecularExponent',25, ...
'SpecularColorReflectance',0.5,'BackFaceLighting','reverselit');
             set( Fig_handles.h_T(1:2),...
    'FaceLighting','phong',...
'AmbientStrength',.4,'DiffuseStrength',.8,...
'SpecularStrength',.2,'SpecularExponent',25, ...
'SpecularColorReflectance',0.5,'BackFaceLighting','reverselit');


             material(Fig_handles.h_R,'shiny');
material( Fig_handles.h_T(1:2),'shiny');
        
        mat_size=size(CV_mesh);
        
        X=reshape(CV_mesh(ss,1,:,:),mat_size(3:4));
        Y=reshape(CV_mesh(ss,2,:,:),mat_size(3:4));
        Z=reshape(CV_mesh(ss,3,:,:),mat_size(3:4));
        %         zeroval=zeros(size(Z))
        Fig_handles.Sub=mesh(X,Y,Z,'Facecolor','c','FaceAlpha',0.4);
        
        
        
        XYZ_int = textscan(fileID4,'%f %f %f','Delimiter',',','HeaderLines',0) ;
        XYZ_int=cell2mat(XYZ_int);
        Fig_handles.Laser_int_points=plot3( XYZ_int(:,1), XYZ_int(:,2), XYZ_int(:,3) ,'k.');
        
        XYZ_Rays = textscan(fileID5,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        XYZ_Rays=cell2mat(XYZ_Rays);
        Fig_handles.Graphics_laser= plot3( XYZ_Rays(:,1), XYZ_Rays(:,2), XYZ_Rays(:,3) ,'g:','Linewidth',0.5 );
        
        Normal_vectors_Mandrel = textscan(fileID6,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Normal_vectors_Mandrel=cell2mat(Normal_vectors_Mandrel);
        Fig_handles.Normal_lines1=plot3( Normal_vectors_Mandrel(:,1), Normal_vectors_Mandrel(:,2), Normal_vectors_Mandrel(:,3) ,'y:','Linewidth',0.25);
        
        
        Reflection_vector_Mandrel = textscan(fileID7,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Reflection_vector_Mandrel=cell2mat(Reflection_vector_Mandrel);
        Fig_handles.Graphics_reflection1=plot3( Reflection_vector_Mandrel(:,1), Reflection_vector_Mandrel(:,2), Reflection_vector_Mandrel(:,3) ,'r--','Linewidth',0.5);
        
        
        Normal_vectors_Tape = textscan(fileID8,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Normal_vectors_Tape=cell2mat(Normal_vectors_Tape);
        Fig_handles.Normal_lines2= plot3( Normal_vectors_Tape(:,1), Normal_vectors_Tape(:,2), Normal_vectors_Tape(:,3) ,'y:');
        
        
        Reflection_vector_Tape = textscan(fileID9,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Reflection_vector_Tape=cell2mat(Reflection_vector_Tape);
        Fig_handles.Graphics_reflection2=plot3( Reflection_vector_Tape(:,1), Reflection_vector_Tape(:,2), Reflection_vector_Tape(:,3) ,'b--');
        
        Ref_int_xyz = textscan(fileID11,' %f %f %f  ','Delimiter',',','HeaderLines',0) ;
        Ref_int_xyz=cell2mat(Ref_int_xyz);
        Fig_handles.reflection_int=plot3( Ref_int_xyz(:,1), Ref_int_xyz(:,2), Ref_int_xyz(:,3) ,'y.');
        
        
        
        
        Fig_handles.P_sub_in=plot3( P_sub_in(:,1), P_sub_in(:,2), P_sub_in(:,3) ,'cs');
        
        
        %% Changing the graphical object, Visibility or deleting
        
        %     Graphics_laser=findobj(h,'color','g','LineStyle',':');
        assignin('base','Graphics_laser',Fig_handles.Graphics_laser);
        % set(Graphics_laser,'Visible','off');
        %     Graphics_reflection1=findobj(h,'color','r');
        %     Graphics_reflection2=findobj(h,'color','b');
        Fig_handles.Graphics_reflection=[Fig_handles.Graphics_reflection1;Fig_handles.Graphics_reflection2];
        assignin('base','Graphics_reflection',Fig_handles.Graphics_reflection);
        
        %     reflection_int=findobj(h,'color','y','Marker','.'); % intersection from the reflection
        assignin('base','reflection_int',Fig_handles.reflection_int);
        
        assignin('base','P_sub_in',Fig_handles.P_sub_in);
        
        
        %     Laser_int_points=findobj(h,'color','k','Marker','.'); % intersection from the reflection
        assignin('base','Laser_int_points',Fig_handles.Laser_int_points);
        assignin('base','h_R',Fig_handles.h_R);
        
        
        %     Normal_lines=findobj(h,'color','y','LineStyle',':'); % intersection from the reflection
        Fig_handles.Normal_lines=[Fig_handles.Normal_lines1;Fig_handles.Normal_lines2];
        assignin('base','Normal_lines',Fig_handles.Normal_lines);
        
        % Thermal_points=findobj(h,'color','m','Marker','d'); % intersection from the reflection
        % assignin('base','Thermal_points',Thermal_points);
        
        Fig_handles.Tape_sub=[Fig_handles.Sub,Fig_handles.h_T];
        assignin('base','Tape_sub',Fig_handles.Tape_sub);
        
        
        
        
        % set(Graphics_reflection,'Visible','off');
        
        
        %        rb5 = uicontrol('Style', 'radiobutton','Position',[520 50 90 40],'Value',true,...
        %         'String','Tape-Substrate','Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]),'Backgroundcolor','w');
        
        %     btn5 = uicontrol('Style', 'pushbutton', 'String', 'Clear Tape-Sub',...
        %         'Position', [520 20 90 20],...
        %         'Callback', 'delete(Tape_sub)','Backgroundcolor','w');
        
        
        
        set(rb1,'Callback', @(rb1,event) Visible_plotButtonPushed(rb1,Fig_handles.Graphics_laser));
        set(rb2,'Callback', @(rb2,event) Visible_plotButtonPushed(rb2,Fig_handles.Graphics_reflection));
        % mandrel is not updating
        %set(rb3,'Callback', @(rb3,event) Visible_plotButtonPushed(rb1,Graphics_laser));
        set(rb4,'Callback', @(rb4,event) Visible_plotButtonPushed(rb4,Fig_handles.h_R));
        
        set(rb5,'Callback', @(rb5,event) Visible_plotButtonPushed(rb5,[Fig_handles.Sub,Fig_handles.h_T]));
        
        set(rb6,'Callback', @(rb6,event) Visible_plotButtonPushed(rb6,Fig_handles.reflection_int));
        set(rb7,'Callback', @(rb7,event) Visible_plotButtonPushed(rb7,Fig_handles.Laser_int_points));
        set(rb8,'Callback', @(rb8,event) Visible_plotButtonPushed(rb8,Fig_handles.Normal_lines));
        
        set(rb10,'Callback', @(rb10,event) Visible_plotButtonPushed(rb10,Fig_handles.P_sub_in));
        
        set(rb9,'Callback', @(rb9,event) Visible_plotButtonPushed(rb9,Tape_Rack));
        
        
        fclose(fileID1);
        fclose(fileID2);
        fclose(fileID3);
        fclose(fileID4);
        fclose(fileID5);
        fclose(fileID6);
        fclose(fileID7);
        fclose(fileID8);
        fclose(fileID9);
        fclose(fileID11);
        
          fclose(fileID21);
        
        
        
        % 4 transparency of Setup.fig
        data=get(gca,'ch');
        
        for i=1:7
            data(i).Color(4) = get(sld_transparency,'value');
            
        end
        
        
        guidata(gcbo,Fig_handles) ;
        
    end
end



% %%
%
% figure(33);
% subplot(2,1,1);
%   Y_unfolded= linspace(0,w,node_num_Z);
%
%     X_unfolded=linspace(0,L_prim-Nip_Mov,No_dev_L_Cut);
%
%
%
%     delta_x=X_unfolded(2)-X_unfolded(1);
%     delta_y=Y_unfolded(2)-Y_unfolded(1);
%     ynode=node_num_Z;
%     xnode=No_dev_L_Cut;
%
%
%     X=zeros(xnode,ynode);
%     Y=zeros(xnode,ynode);
%     Z=zeros(xnode,ynode);
%     Intensity=zeros(xnode,ynode);
%
%     for kk=1:xnode
%         index=((1:ynode)*xnode)-kk+1;
%         X(kk,1:end)=points(index,1)';
%         Y(kk,1:end)=points(index,2)';
%         Z(kk,1:end)=points(index,3)';
%         Intensity(kk,1:end)=E_points (index)'; % RHS in the thermal model
%     end
%     surf(X,Y,Z,Intensity);
%      axis equal;
%
%     % Tape
%        L_prim_T=R_tape*(deg-theta_ind)*(pi/180) +L_flat;
%
%
%
%
%
%
%
%
%
%      xnode_T=node_num_Z_T;
%     ynode_T=length  (th_v)+length(dl);
%
%     X_unfolded_T= linspace(0,W_tape,xnode_T);
%     Y_unfolded_T=linspace(0,L_prim_T,ynode_T);
%     Lx_T=W_tape;
%     Ly_T=L_prim_T;
%
%     delta_x_T=X_unfolded_T(2)-X_unfolded_T(1);
%     delta_y_T=Y_unfolded_T(2)-Y_unfolded_T(1);
%
%     % make number of xnode and ynode for Tape Thermal analysis
%
%
%
%
%     X_T=zeros(xnode_T,ynode_T);
%     Y_T=zeros(xnode_T,ynode_T);
%     Z_T=zeros(xnode_T,ynode_T);
%     Intensity_T=zeros(xnode_T,ynode_T);
%
%     for kk=1:xnode_T
%         index=((1:ynode_T)*xnode_T)-kk+1;
%         X_T(kk,1:end)=points_T_G(index,1)';
%         Y_T(kk,1:end)=points_T_G(index,2)';
%         Z_T(kk,1:end)=points_T_G(index,3)';
%         Intensity_T(kk,1:end)=E_points_T(index)';
%     end
%     subplot(2,1,2);
%     surf(X_T,Y_T,Z_T,Intensity_T);
%     axis equal