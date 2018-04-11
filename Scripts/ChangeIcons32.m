% Script to restore the desired order in overlay icons
% Copies the relevant registry, modifies it and replaces. A copy of both the 
% original and the altered is kept in case something went wrong and it is
% desired to roll back. Two batch files are created in line and executed
% through matlab. At the end explorer.exe is restarted. This script should
% be in a path with no spaces since they cannot be handled by DOS.
% Viva la revolution!

% VPE 3/5/2017
clc, clear all  %#ok<CLALL>

%% Choose the overlays to be kept on top
ActiveIcons={  'DropboxExt01'
               'DropboxExt02'
               'DropboxExt05'
               'DropboxExt07'               
               'Tortoise1Normal'
               'Tortoise2Modified'
               'Tortoise3Conflict'
               'Tortoise7Added'
               'Tortoise8Ignored'
               'Tortoise9Unversioned'
               'Tortoise4Locked'
               'Tortoise5ReadOnly'
               'Tortoise6Deleted'
               };
         
%% Create and run batch
path=pwd;
BATstr=['REGEDIT /e ' path '\OrigReg32.reg "HKEY_LOCAL_MACHINE\SOFTWARE\Wow6432Node\Microsoft\Windows\CurrentVersion\explorer\ShellIconOverlayIdentifiers"'];
batfile=[path '\Overlay32.bat'];
fid = fopen(batfile, 'w');
fprintf(fid,'%s\n', BATstr);
fprintf(fid,'%s\n', 'exit');
fclose(fid);
system( 'Overlay32.bat');

%% read in original registry file
warning off MATLAB:iofun:UnsupportedEncoding;
fid = fopen('OrigReg32.reg', 'r', 'l', 'UTF-16');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
%     tline = fgetl(fid);
    tline = fgetl(fid );
    
    if ischar(tline)
        A{i} = tline;
    end
end
fclose(fid);

%% Modify and rewrite registry file
fid = fopen('Modified32.reg', 'w' );
for i = 1:numel(A)
        if i<2 % dont modify header
            newstr = A{i};
        else
            curstr = A{i};
            if isempty(strfind(curstr,'\'))==0
               IndSlash   = strfind(curstr,'\');
               KeyNameInd = IndSlash(end);
               cnt        = 1;
               found      = 0;
               while found==0
                    if strcmp(curstr(KeyNameInd+cnt),' ') == 1
                        cnt = cnt+1;                   
                    else
                        StrLim = KeyNameInd+cnt;
                        found = 1;
                    end
               end
                newstr  =  [regexprep(curstr(1:(StrLim-1)),' +','') curstr((StrLim):end)];
            else
                newstr  = regexprep(curstr,' +','');  %remove all spaces 
            end
            % Check if the targeted name is included in current row and add a space to give priority            
            for j = 1:length(ActiveIcons)                 
               if isempty(strfind (newstr,ActiveIcons{j}))
                   CheckVec(j) = 0;
               else
                   CheckVec(j) = 1;                    %#ok<*SAGROW>
               end
            end
            if any(CheckVec)
                spaceind = strfind (newstr,ActiveIcons(CheckVec==1)) ;
                newstr(spaceind+1:end+1) = newstr(spaceind:end);
                newstr(spaceind) = ' '   ;            
            end
            
            if strcmp(newstr,curstr)==0
                disp (['modified key ' curstr ' to ' newstr])                   
            end
        end       
        fprintf(fid,'%s\n', newstr);  
        clear CheckVec
end
fclose(fid);

%% Import modify registry replacing the old one

Regstr1 = 'Windows Registry Editor Version 5.00';
Regstr2 = '';
Regstr3 = ['[-HKEY_LOCAL_MACHINE\SOFTWARE\Wow6432Node\Microsoft\Windows\CurrentVersion\explorer\ShellIconOverlayIdentifiers]']; %#ok<*NBRAK>
Regfile1=['RemoveReg32.reg'];
fid = fopen(Regfile1, 'w');
fprintf(fid,'%s\n', Regstr1);
fprintf(fid,'%s\n', Regstr2);
fprintf(fid,'%s\n', Regstr3);
fclose(fid);
system ('RemoveReg32.reg');
system ('Modified32.reg');

%% restart explorer.exe 
 system ('Taskkill /IM explorer.exe /F');
 system ('start explorer.exe');

