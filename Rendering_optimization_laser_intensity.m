


function   Rendering_optimization_laser_intensity(handles,...
              Laser_head);
Ax=Laser_head(1);
Ay=Laser_head(2);
nx=Laser_head(3);
ny=Laser_head(4);


step_SF(1)=str2num(handles.step_UOT_optimization{1});
step_SF(2)=str2num(handles.step_UOT_optimization{2});



Video = VideoWriter(strcat(handles.UOT_pathfile, sprintf('Laser_Intensity')),'MPEG-4');
% video.FrameRate = 10;
Video.FrameRate=2; Video.Quality = 99;


open(Video);

% UOT_pathfile='C:\Users\ZaamiA\Desktop\OTOM-Dev-17Juli-2019\Analysis_UOT\Dome89-deg-4Opt\';
  fileID30 = fopen(strcat(handles.UOT_pathfile, sprintf('Selected_Laser_ID.txt')),'r');
  
  ID= textscan(fileID30,' %f %f %f %f %f %f %f %f','Delimiter',',','HeaderLines',0) ;
    ID=cell2mat(ID);
[m,n]=size(ID);


% step_SF=[60 110];
figure;
xp_laser=linspace(-Ax,Ax,nx);
yp_laser=linspace(-Ay,Ay,ny);

[X,Y]=meshgrid(xp_laser,yp_laser);
for jj=step_SF(1):step_SF(2)




% for ii=1:1
Power_Actual=Laser_Power_generator (Ax,Ay,nx,ny,ID(jj-step_SF(1)+1,1:end-1));

 surf(X,Y,Power_Actual,'linestyle','none');

  view([0 -90]);
  colorbar;
   title(sprintf('Step=%d, Total laser power= %f [W]',jj,ID(jj-step_SF(1)+1,end)));
axis equal
grid off;
set(gcf,'color','w');
   javaFrame    = get(gcf,'JavaFrame');
    iconFilePath = 'OTOM-icon.png';
    javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));
    
 pause(0.01);
% end
    Mov(jj-step_SF(1)+1)=getframe(gcf);
%            pause(0.01);
        

%         writeVideo(Video,Mov);


end

 writeVideo(Video,Mov);

fclose(fileID30);
  close(Video);
