function [ varargout ] = ct_find_event_list( varargin )
%FIND_EVENT_LIST_PHASE Summary of this function goes here
%   Detailed explanation goes here

% find the event location
% example:
% find the index of nearest existed event with non-empty file list
% ind = ct_find_event_list( ievent, events,listname, listout, 'backward')
% ievent: the index of the reference event;
% events: a cell contains the path of all events
% listname: file list
% listout: out list file 
% (if listout not exist, read listname, else listout)
% backward: backward search; three options: 'forward','backward','exact'; default is 'exact'



varargout = {[],[]};
if nargin < 3
    fprintf('Number of input arguments less than 3!\n')
    return;
elseif nargin == 3
    ievent = varargin{1};
    events = varargin{2};
    listname = varargin{3};
    listout = [];
    search_type = 'exact';
elseif nargin == 4
    ievent = varargin{1};
    events = varargin{2};
    listname = varargin{3};
    listout = varargin{4};
    if strcmpi(listout,'forward') || strcmpi(listout,'backward') || strcmpi(listout,'exact')
        search_type = listout;
        listout = [];
    else
        search_type = 'exact';
    end
else
    ievent = varargin{1};
    events = varargin{2};
    listname = varargin{3};
    listout = varargin{4};
    search_type = varargin{5};
end

if ~isnumeric(ievent)
    fprintf('First parameter should be a number\n');
    return;
end

nevent = length(events);
% if ievent out of bound
if  ievent > nevent || ievent < 1
    fprintf('Event index exceeds the total event number!\n');
    return;
end

ind = ievent;
if strcmpi(search_type,'forward') % forward search
%     if ievent == nevent
%         fprintf('This is the last event!\n');
%     end

    while ind <= nevent
        if exist(char(events(ind)),'dir')
            if exist(fullfile(char(events(ind)),listname),'file')
                if ~exist(fullfile(char(events(ind)),listout),'file')
                    fid = fopen(fullfile(char(events(ind)),listname),'r');
                    C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                    Ccell = C{1};
                    fclose(fid);
                    if ~isempty(Ccell)
                        varargout{1} = ind;
                        varargout{2} = 1;
                        return;
                    end
                else
                    fid = fopen(fullfile(char(events(ind)),listout),'r');
                    C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                    Ccell = C{1};
                    fclose(fid);
                    if ~isempty(Ccell)
                        varargout{1} = ind;
                        varargout{2} = 2;
                        return;
                    end
                end
            end
        end
%         fprintf('No data found in %s. Go to next event!\n',char(events(ind)));
        ind = ind + 1;
    end
%     if ind > nevent
%         fprintf('No next event found!\n');
%     end
elseif strcmpi(search_type,'backward') % backward search
%     if ievent == 1
%         fprintf('This is the first event!\n');
%     end
    while ind >= 1
        if exist(char(events(ind)),'dir')
            if exist(fullfile(char(events(ind)),listname),'file')
                if ~exist(fullfile(char(events(ind)),listout),'file')
                    fid = fopen(fullfile(char(events(ind)),listname),'r');
                    C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                    Ccell = C{1};
                    fclose(fid);
                    if ~isempty(Ccell)
                        varargout{1} = ind;
                        varargout{2} = 1;
                        return;
                    end
                else
                    fid = fopen(fullfile(char(events(ind)),listout),'r');
                    C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                    Ccell = C{1};
                    fclose(fid);
                    if ~isempty(Ccell)
                        varargout{1} = ind;
                        varargout{2} = 2;
                        return;
                    end
                end
            end
        end
%         fprintf('No data found in %s. Go to previous event!\n',char(events(ind)));
        ind = ind - 1;
    end
%     if ind > nevent
%         fprintf('No previous event found!\n');
%     end
elseif strcmpi(search_type,'exact') % exact the specified event
    if exist(char(events(ind)),'dir')
        if exist(fullfile(char(events(ind)),listname),'file')
            if ~exist(fullfile(char(events(ind)),listout),'file')
                fid = fopen(fullfile(char(events(ind)),listname),'r');
                C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                Ccell = C{1};
                fclose(fid);
                if ~isempty(Ccell)
                    varargout{1} = ind;
                    varargout{2} = 1;
                    return;
                end
            else
                fid = fopen(fullfile(char(events(ind)),listout),'r');
                C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
                Ccell = C{1};
                fclose(fid);
                if ~isempty(Ccell)
                    varargout{1} = ind;
                    varargout{2} = 2;
                    return;
                end
            end
        end
    end
else
    fprintf('Search type must be: forward, backward, or exact\n');
end



end

