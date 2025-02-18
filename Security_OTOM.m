
% codes for get mac address and check with existed data in the websie

function [Out,Access_level]=Security_OTOM(account)
% read name of the PC

Out=0; % default answer
Access_level=zeros(1,15);

[~,PC_name]=system('hostname');
PC_name=deblank(PC_name);
% cpu serial number
[~,CPU_name]=system('wmic bios get serialnumber');

% mac address
[~,status]=system('getmac');
% read whole sentences , 72 is the length of whole sentence
mac_all=strsplit(status,'\n');
mac_all=mac_all(4:end);

% mac_all=status(160:end);  > old , does not work based on all windows
% languages

% status(160:160+72-1)
% status(160+72-1+8:)

% current_mac_temp=strsplit(mac_all,'\n');

current_mac_temp=mac_all;

% deblank(current_mac(1))
current_mac=deblank(current_mac_temp(1:2));
% status(160:160+71)

% check how many Mac adress the PC has, if only 1 existed, use CPU serial number to ckeck

% remove space among characters
  current_mac{1}= current_mac{1}(~isspace(current_mac{1}));
  current_mac{2}= current_mac{2}(~isspace(current_mac{2}));
  CPU_name= CPU_name(~isspace(CPU_name));
  
if length(current_mac_temp) ==2
  
%     temp=deblank(CPU_name);
%     temp2=strsplit(temp,'\n');
  
      current_mac{2}=  CPU_name;
end




% In account section of OTOM

web_address=strcat('https://www.otomcomposite.com/wp-content/uploads/',account,'.xlsx');
% read from data file from website
% data=webread('https://www.otomcomposite.com/wp-content/uploads/Mac_Address.xlsx');
try
    data=webread(web_address);
    %convert to array
    cell_data=data{:,:};
    
    
    
    % temp = regexp(fileread('C:\Users\zaamiA\Desktop\Mac_Address.txt'), '\r?\n', 'split');
    
    % find a specified PC
    Index = find(contains(cell_data,PC_name));
    
    if Index
        
        % remove the sapce between characters
         cell_data{Index+1}= cell_data{Index+1}(~isspace(cell_data{Index+1}));
        cell_data{Index+2}= cell_data{Index+2}(~isspace(cell_data{Index+2}));
        
        Logic1=length(cell_data{Index+1})== length(current_mac{1}) ;
        Logic2=length(cell_data{Index+2})== length(current_mac{2}) ;
        
        if  Logic1 && Logic2
            
            if [cell_data{Index+1}== current_mac{1} , cell_data{Index+2}== current_mac{2}]
                h = msgbox('Licence is verified!','Success');
                Access_level=str2num(cell_data{Index,2});
                
                Out=1;
            else
               h = errordlg('Licence is Not verified!');
                Out=0;
            end
            
        else
           h = errordlg('Licence is Not verified!');
            Out=0;
        end
    else
        h =errordlg('Licence is Not verified!');
        Out=0;
    end
catch
    h =errordlg('Your account is not registered!','File Error');
    Out=0;
end

javaFrame    = get(h,'JavaFrame');
iconFilePath = 'OTOM-icon.png'; 
javaFrame.setFigureIcon(javax.swing.ImageIcon(iconFilePath));

%%
% code for linux
% system('/sbin/ifconfig eth0'); linux


%% send E-mail for pc-registration
% setpref('Internet','SMTP_Server','https://xs.utwente.nl');
%             setpref('Internet','E_mail','a.zaami@utwente.com');

% should be placed in Share section of OTOM
% message=sprintf('%s <br> \r\n %s  <br> \r\n %s  <br> \r\n %s',PC_name , mac_all(1:end/2),mac_all(1+(end/2):end), CPU_name );
%
% sendolmail('info@otomcomposite.com','PC-registration in OTOMcomposite',message);




