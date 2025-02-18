
function Temp_Red_Box=Measuring_Box_Cal (T,Measure_Box,delta_x,delta_y,index_middle_long)



c=Measure_Box(1);
Blx=Measure_Box(2);
Bly=Measure_Box(3);


LLX= floor(((c-Blx )/delta_x ));
ULX=ceil((c+Blx )/delta_x );
LLY=floor( -(Bly /delta_y));
ULY=ceil( +(Bly /delta_y));

Mesure_Reader_mid_X=index_middle_long(end-ULX:end-LLX);

x_dir_len=length (Mesure_Reader_mid_X);
Offset_Y_Reader=(LLY:ULY);
Offset_Y_Reader_len=length(Offset_Y_Reader);
Measuring_index=zeros(Offset_Y_Reader_len,x_dir_len);

for ii=1:Offset_Y_Reader_len
Measuring_index(ii,1:x_dir_len)=Mesure_Reader_mid_X+Offset_Y_Reader(ii);
end


Temp_Red_Box=(mean(T(Measuring_index(1:end))));

% disp(Temp_Red_Box);