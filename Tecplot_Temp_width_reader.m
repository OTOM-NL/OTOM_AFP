


fileID1 = fopen('Temperature_Vel.plt');

fileID1 = fopen('Temperature_Vel_Tape.plt');
% node = textscan(fileID, '%*d %f %f  ','delimiter', ',','commentStyle', 'NODE');


Temp_cell = textscan(fileID1,' %f %f %f  ','Delimiter',',','HeaderLines',3) ;
Temp=cell2mat(Temp_cell);

Temp_width=zeros(10,1);
counter=0;

for ii=1:length(Temp)
   if(Temp(ii,1)==0)
       counter=counter+1;
    Temp_width(counter)=Temp(ii,3);
   end
    
end

x=linspace(0,1,length(Temp_width));


plot(x,Temp_width);
hold on;
xlabel('Normalized tape width');
ylabel('Temperature {\circ}C');
title('Temperature along Nip-point Line');

