
%% Input dialogue for each section input


prompt = {'\theta_y:','Tape length','3 \circ C','4','5','6','7','8','9','10','11'};
dlg_title = 'Geometrical parameters';
num_lines = 1;
defaultans = {'20','56','3','4','5','6','7','8','9','10','11,30, 25'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);

[val status] = cell2mat(answer);  

%%
prompt = {'Velocity:','Laser intensity','3','4','5','6','7','8','9','10','11'};
dlg_title = 'Process parameters';
num_lines = 1;
defaultans = {'20','hsv','3','4','5','6','7','8','9','10','11'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);

%%

prompt = {'Laser source position Laser source position Laser source position','Laser direction','roller radius','roller width','tape width','tape long','mandrel radius','substarte long','9','10','11'};
dlg_title = 'Geometrical parameters';
num_lines = 1;
defaultans = {'10, 20, 40','0.02 0.01 0.07','3','4','5','6','7','8','9','10','11'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,options);
