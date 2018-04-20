function [  ] = ct_gen_listname(  )

% generate evlist and listname for crazytremor
parentdir = pwd;
eventdir = strcat(parentdir,'/SAC');
path(path,genpath(fullfile(parentdir,'utilities')));

% generate evlist and listname for crazytremor
eventid = '*';
evlist = 'evlist.txt';
listname = 'list';

% note: component id should be unique. It is the only way to identify
% components for files. 
n = 0;
n = n + 1;
cmpname{n} = 'Z';
cmpid{n} = {'*.BHZ.*SAC','*.Z.*SAC','*.EHZ.*SAC','*.HHZ.*SAC','*.SHZ.*SAC','*.U.*SAC'};
n = n + 1;
cmpname{n} = 'N';
cmpid{n} = {'*.BHN.*SAC','*.N.*SAC','*.EHN.*SAC','*.HHN.*SAC','*.SHN.*SAC','*.EH1.SAC'};
n = n + 1;
cmpname{n} = 'E';
cmpid{n} = {'*.BHE.*SAC','*.E.*SAC','*.EHE.*SAC','*.HHE.*SAC','*.SHE.*SAC','*.EH2.SAC'};
n = n + 1;
cmpname{n} = 'R';
cmpid{n} = {'*.BHR.*SAC','*.R.*SAC','*.EHR.*SAC','*.HHR.*SAC','*.SHR.*SAC'};
n = n + 1;
cmpname{n} = 'T';
cmpid{n} = {'*.BHT.*SAC','*.T.*SAC','*.EHT.*SAC','*.HHT.*SAC','*.SHT.*SAC'};

events = dir(fullfile(eventdir,eventid));
fid = fopen(evlist,'w');
for i = 1:length(events)
    if strcmpi(events(i).name,'.') || strcmpi(events(i).name,'..')
        continue;
    end
    
    fprintf(fid,'%s\n',fullfile(eventdir,events(i).name));
    
    tr = [];
    num = 0;
    for j = 1:length(cmpname)
        for jj = 1:length(cmpid{j})
            files = dir(fullfile(eventdir,events(i).name,cmpid{j}{jj}));
            
            for k = 1:length(files)
                
                hd = irdsachd(fullfile(eventdir,events(i).name,files(k).name));
                
                if isempty(hd)
                    continue;
                end
                
                num = num + 1;
                tr(num).netwk = deblank(reshape(hd.knetwk,1,[]));
                tr(num).stnm = deblank(reshape(hd.kstnm,1,[]));
                tr(num).nstnm = [tr(num).netwk,'.',tr(num).stnm];
                tr(num).stla = hd.stla;
                tr(num).stlo = hd.stlo;
                tr(num).stel = hd.stel;
                tr(num).cmp = cmpname{j};
                tr(num).filename = files(k).name;
            end
        end
    end
    
    % find unique station names
    [nstnm,ia,ic] = unique({tr.nstnm});
    fid1 = fopen(fullfile(eventdir,events(i).name,listname),'w');
    fprintf(fid1,'#Netwk Stnm Stla Stlo Quality Group * Cmp1 File1 Cmp2 File2 ... *\n');
    for j = 1:length(nstnm)
        fprintf(fid1,'%s %s %f %f 0 U * ',tr(ia(j)).netwk,tr(ia(j)).stnm,tr(ia(j)).stla,tr(ia(j)).stlo);
        for k = 1:length(cmpname)
            ind = find(ic==j);
            idx = find(strcmpi(cmpname{k},{tr(ind).cmp}));
            if isempty(idx)
                fprintf(fid1,' %s EMPTY',cmpname{k});
            else
                idx = idx(1); % if there are multiple traces for each component of each station, choose the first one.
                fprintf(fid1,' %s %s',cmpname{k},tr(ind(idx)).filename);
            end
        end
        fprintf(fid1,' *\n');
    end
    fclose(fid1);
end
                           
fclose(fid);

end

