function [  ] = y_CrazyTremor_v1(  )
%
% CrazyTremor for triggered tremor analysis
% Please see online Manual for detailed description of the software:
% http://www.kevinchao.com/crazytremor
%
% version 1.1 - Nov. 13, 2017
% version 1.0 - Sep. 24, 2015
%
% Written by Chunquan Yu and Kevin Chao
% Contact: yucq@gps.caltech.edu or kchao@northwestern.edu
%
% Reference: 
% Chao, K. and C. Yu, A MATLAB GUI-based Package for Examining Triggered Tremor: 
% A Case Study in New Zealand (2018), Seismol. Res. Lett., under review.
%

%% add path
parentdir = pwd;
path(path,genpath(fullfile(parentdir,'utilities')));


%% set parameters
para.ftsize = 12;
para.axesftsize = 12;
para.panelftsize = 14;
para.ftname = 'Arial';
set(0,'defaultaxesfontsize',para.axesftsize);
set(0,'defaulttextfontsize',para.axesftsize);
set(0,'DefaultUicontrolFontsize',para.ftsize);
set(0,'DefaultUicontrolFontname',para.ftname);
set(0,'DefaultUicontrolUnits','normalized');

% Uicontrol background color
para.bcolor = 0.9*[1 1 1]; % edit box background color

% Trace window
para.linewidth_show = 1; % line width of regular traces
para.linewidth_selected = 2; % line width of selected traces
% Arrival time - plot parameters
para.tmark_linewidth = 1;
para.tmark_color = 'k';
para.tmark_fontsize = 12;

% Map window
% station plot parameters
para.markersize_st = 8;
para.markercolor_st = [0.7 0.7 0.7];
para.markerlinewidth_st = 2;
% statons in box
para.markersize_box = 10;
para.markerlinewidth_box = 2;
% stations selected
para.markersize_sel = 16;
para.markerlinewidth_sel = 3;
% Location - plot parameters
para.markersize_loc = 8;
para.marker_loc = 'p';
para.markercolor_loc = 'k';
para.markercolor_noloc = [0.7 0.7 0.7];
% reference point (for sorting)
para.markercolor_refpoint = 'b';
para.markersize_refpoint = 8;
para.marker_refpoint = '+';
para.markerlinewidth_refpoint = 2;

% set defaults
para.evlist = fullfile(pwd,'evlist.txt');
[evpathname, name1, name2] = fileparts(para.evlist);
para.evpathname = evpathname;
para.evlistname = [name1,name2];

para.listname = 'list';
para.outlistid = '_o';
para.outlist = [para.listname,para.outlistid];

% output location file
para.outlocfile = 'loc.txt';

para.ievent = 1;
para.iframe = 1;
para.nframes = 1;

% number of traces per frame
para.n_per_frame = 10;

% sort type
para.sortlist = {'stla','stlo','point'};

% filter list
para.filterlist = {'bp','lp','hp'}; % band pass, low pass, high pass
% filtering parameters for traces in panel 1 - usually higher frequency
para.fl1 = 2;
para.fh1 = 8;
para.order1 = 2;
% filtering parameters for traces in panel 2 - usually lower frequency
para.fl2 = 0;
para.fh2 = 0;
para.order2 = 2;

% delta for plotting
para.delta = 0.1;

% time window
para.timewin = [0 0];

% plot scale
para.scale = 1;

% is plot envelope
para.isenv = 0;

% is log plot for amplitude
% this is only valid for envelope plot
para.issemilogy = 0;

% is smooth traces
% this is only valid for envelope plot
para.issmooth = 0;
para.smoothtime = 10; % smooth time length in seconds

% is spectrogram for subplot2
para.isspectrogram2 = 0;
% spectrogram plot parameters
para.spec_twinl = 10;
para.spec_tstep = 5;

% all components
para.cmps = {'Z','N','E','R','T'};
para.cmps_colors = {[0 0 0],'m','c','r','g'};
% para.cmps_colors_sel = {'k','k','k','k','k'};
% para.cmps_colors = {[1 0 0],[0 1 0],[0 0 1],[0.5 0 0.9],[0.9 0.5 0]};
para.cmps_colors_sel = {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]};

% quality and group list
para.quality_list = {'0','1','2','3','4','9'};
para.group_list = {'U','A','B','C','D','E'};
para.marker_quality = {'o','^','s','d','+','x'};
para.quality_liststr = {'0 - o','1 - ^','2 - s','3 - d','4 - +','9 - x'};
para.markercolor_group = {[0 0 0],[1 0 0],[0 1 0],[0 0 1],[0.5 0 0.9],[0.9 0.5 0]};

% base map
% tmp = load('coast');
tmp = load('coastlines');
para.long = tmp.coastlon;
para.lat = tmp.coastlat;

% for location
para.vmodelfile = fullfile(parentdir,'utilities','vmodels','vmodel1.txt');
% search region beyond the current station lat/lon limits
para.ostla = 0.5;
para.ostlo = 0.5;
% grid interval
para.ddelta = 0.01;
% default search depth
para.xdep = 35;

tr = [];

%
pre = '<HTML><FONT color="';
post = '</FONT></HTML>';
para.group_liststr = cell(length( para.markercolor_group ),1);
for tmp0 = 1:length(para.markercolor_group)
   str = [pre reshape( dec2hex(round(255*para.markercolor_group{tmp0}), 2 )',1, 6) '">' para.group_list{tmp0} post];
   para.group_liststr{tmp0} = str;
end

%%
handle.f1               =	figure('Toolbar','figure','Units','normalized','keypressfcn', @CT_short_cut,'Position',[.1 .1 .8 .8]);
                            set(handle.f1,'name','CrazyTremor','NumberTitle','off');
handle.hax1             =	axes('pos', [.15 .3 .5 .6]);
handle.hax2             =	axes('pos', [.15 .075 .5 .15]);
handle.hax3             =	axes('pos', [.78 .5 .2 .4]);

% hot keys
handle.h_hot            =	uicontrol('String','<html><b>Crazy Tremor</b></html>','keypressfcn',@CT_short_cut,'callback',@CT_copyright,'Position',[.02 .91 .1 .08]);

%% control panel for subplot 1
handle.control1         =   uipanel('parent', handle.f1,'title', 'Plot', 'Position', [0.02 .3 .1 .6],'fontsize',para.panelftsize);

% plot data
handle.h_iniplot           =   uicontrol(handle.control1,'String','Ini Plot(i)','callback',@CT_callback_iniplot,'Position',[.1 .90 .8 .1]);
% components to show
                            uicontrol(handle.control1,'Style','text','String','Z','Position',[.0 .84 .2 .05]);
                            uicontrol(handle.control1,'Style','text','String','N','Position',[.2 .84 .2 .05]);
                            uicontrol(handle.control1,'Style','text','String','E','Position',[.4 .84 .2 .05]);
                            uicontrol(handle.control1,'Style','text','String','R','Position',[.6 .84 .2 .05]);
                            uicontrol(handle.control1,'Style','text','String','T','Position',[.8 .84 .2 .05]);
handle.h_cmp_Z          =   uicontrol(handle.control1,'Style','checkbox','Value',1,'callback',@CT_callback_ZNERT,'Position',[.0 .81 .2 .03],'BackgroundColor',para.cmps_colors{1});
handle.h_cmp_N          =   uicontrol(handle.control1,'Style','checkbox','Value',1,'callback',@CT_callback_ZNERT,'Position',[.2 .81 .2 .03],'BackgroundColor',para.cmps_colors{2});
handle.h_cmp_E          =   uicontrol(handle.control1,'Style','checkbox','Value',1,'callback',@CT_callback_ZNERT,'Position',[.4 .81 .2 .03],'BackgroundColor',para.cmps_colors{3});
handle.h_cmp_R          =   uicontrol(handle.control1,'Style','checkbox','Value',0,'callback',@CT_callback_ZNERT,'Position',[.6 .81 .2 .03],'BackgroundColor',para.cmps_colors{4});
handle.h_cmp_T          =   uicontrol(handle.control1,'Style','checkbox','Value',0,'callback',@CT_callback_ZNERT,'Position',[.8 .81 .2 .03],'BackgroundColor',para.cmps_colors{5});
% number of trace
                            uicontrol(handle.control1,'Style','text','String','Ntrace','Position',[.0 .75 .5 .04]);
handle.h_ntrace_num     =   uicontrol(handle.control1,'Style','edit','String',num2str(para.n_per_frame),'callback',@CT_callback_tracenumber,'Position',[.5 .75 .5 .04],'BackgroundColor',para.bcolor);
% sort traces
                            uicontrol(handle.control1,'Style','text','String','Sort','Position',[.0 .7 .5 .04]);
handle.h_sort_list      =   uicontrol(handle.control1,'Style','popup','String',para.sortlist,'callback',@CT_callback_sort,'Value',1,'Position',[.5 .7 .5 .04]);

% filter 1
                            uicontrol(handle.control1,'Style','text','String','Filter','Position',[.0 .65 .5 .04]);
                            uicontrol(handle.control1,'Style','text','String','Order','Position',[.5 .65 .5 .04]);
handle.h_filter_list1   =   uicontrol(handle.control1,'Style','popup','String',para.filterlist,'Value',1,'callback',@CT_callback_filter1,'Value',1,'Position',[.0 .60 .5 .04]);                           
handle.h_filter_order1  =   uicontrol(handle.control1,'Style','edit','String',num2str(para.order1),'callback',@CT_callback_filter1,'Position',[.5 .60 .5 .04],'BackgroundColor',para.bcolor);
handle.h_filter_fl1     =   uicontrol(handle.control1,'Style','edit','String',num2str(para.fl1),'callback',@CT_callback_filter1,'Position',[.1 .55 .3 .04],'BackgroundColor',para.bcolor);
handle.h_filter_fh1     =   uicontrol(handle.control1,'Style','edit','String',num2str(para.fh1),'callback',@CT_callback_filter1,'Position',[.5 .55 .3 .04],'BackgroundColor',para.bcolor);
                            uicontrol(handle.control1,'Style','text','String','-','Position',[.4 .55 .1 .04]);
                            uicontrol(handle.control1,'Style','text','String','Hz','Position',[.8 .55 .2 .04]);

% time window
                            uicontrol(handle.control1,'Style','text','String','Delta','Position',[.0 .5 .5 .04]);
handle.h_delta          =   uicontrol(handle.control1,'Style','edit','String',num2str(para.delta),'callback',@CT_callback_replot,'Position',[.5 .5 .5 .04],'BackgroundColor',para.bcolor);
                            uicontrol(handle.control1,'Style','text','String','Time win.','Position',[.0 .45 1 .04]);
handle.h_timewin_L      =   uicontrol(handle.control1,'Style','edit','String',num2str(para.timewin(1)),'callback',@CT_callback_replot,'Position',[.0 .4 .5 .05],'BackgroundColor',para.bcolor);
handle.h_timewin_R      =   uicontrol(handle.control1,'Style','edit','String',num2str(para.timewin(2)),'callback',@CT_callback_replot,'Position',[.5 .4 .5 .05],'BackgroundColor',para.bcolor);
% amplitude up and down
                            uicontrol(handle.control1,'Style','text','String','Amplitude','Position',[.0 .35 1 .04]);
                            uicontrol(handle.control1,'String','+(=)','callback',@CT_callback_ampup,'Position',[.0 .3 .5 .05]);
                            uicontrol(handle.control1,'String','-(-)', 'callback',@CT_callback_ampdown,'Position',[.5 .3 .5 .05]);
% envelope
handle.h_envelope       =   uicontrol(handle.control1,'Style','togglebutton','String','Envelope (e)','Value',para.isenv,'callback',@CT_callback_envelope,'Position',[0.0 0.225 1 0.05]);
% semilogy
handle.h_semilogy       =   uicontrol(handle.control1,'Style','togglebutton','String','Log y','Value',para.issemilogy,'callback',@CT_callback_envelope,'Position',[0.0 0.175 0.5 0.05]);
% smooth
handle.h_smooth         =   uicontrol(handle.control1,'Style','togglebutton','String','Smooth','Value',para.issmooth,'callback',@CT_callback_envelope,'Position',[0.5 0.175 0.5 0.05]);

handle.h_timewin_zoom   =   uicontrol(handle.control1,'Style','pushbutton','String','<html>Zoom(z)<br/>(exit q)','callback',@CT_callback_windowzoom,'Position',[.0 .06 .5 .1],'max',2);
handle.h_timewin_full   =   uicontrol(handle.control1,'Style','pushbutton','String','Full(f)','callback',@CT_callback_windowfull,'Position',[.5 .06 .5 .1]);
% input parameters
                            uicontrol(handle.control1,'Style','pushbutton','String','Parameters','callback',@CT_callback_inputpara,'Position',[.0 .0 1 .05]);


%% control panel for subplot 2
handle.control2         =   uipanel('parent', handle.f1,'title', 'Plot2', 'Position', [0.02 .05 .1 .2],'fontsize',para.panelftsize);

% filter 2
                            uicontrol(handle.control2,'Style','text','String','Filter','Position',[.0 .8 .5 .15]);
                            uicontrol(handle.control2,'Style','text','String','Order','Position',[.5 .8 .5 .15]);
handle.h_filter_list2   =   uicontrol(handle.control2,'Style','popup','String',para.filterlist,'Value',1,'callback',@CT_callback_plottrace2,'Value',1,'Position',[.0 .65 .5 .15]);                           
handle.h_filter_order2  =   uicontrol(handle.control2,'Style','edit','String',num2str(para.order2),'callback',@CT_callback_plottrace2,'Position',[.5 .65 .5 .15],'BackgroundColor',para.bcolor);
handle.h_filter_fl2     =   uicontrol(handle.control2,'Style','edit','String',num2str(para.fl2),'callback',@CT_callback_plottrace2,'Position',[.1 .5 .3 .15],'BackgroundColor',para.bcolor);
handle.h_filter_fh2     =   uicontrol(handle.control2,'Style','edit','String',num2str(para.fh2),'callback',@CT_callback_plottrace2,'Position',[.5 .5 .3 .15],'BackgroundColor',para.bcolor);
                            uicontrol(handle.control2,'Style','text','String','-','Position',[.4 .5 .1 .15]);
                            uicontrol(handle.control2,'Style','text','String','Hz','Position',[.8 .5 .2 .15]);

% spectrogram
handle.h_spectrogram2   =   uicontrol(handle.control2,'Style','togglebutton','String','Spectrogram','Value',para.isspectrogram2,'callback',@CT_callback_plottrace2,'Position',[.0 .35 1 .15]);
                            uicontrol(handle.control2,'Style','text','String','Twinl','Position',[.0 .2 .5 .15]);
handle.h_spec_twinl     =   uicontrol(handle.control2,'Style','edit','String',num2str(para.spec_twinl),'Position',[.5 .2 .5 .15],'BackgroundColor',para.bcolor);
                            uicontrol(handle.control2,'Style','text','String','Tstep','Position',[.0 .05 .5 .15]);
handle.h_spec_tstep     =   uicontrol(handle.control2,'Style','edit','String',num2str(para.spec_tstep),'Position',[.5 .05 .5 .15],'BackgroundColor',para.bcolor);

%% control panel for subplot 3
handle.control3         =   uipanel('parent', handle.f1,'title', 'Map', 'Position', [0.78 .3 .2 .15],'fontsize',para.panelftsize);
% load in map
handle.h_loadmap        =   uicontrol(handle.control3,'String','Reset map(m)','callback',@CT_callback_loadmap,'Position',[.525 .525 .45 .45]);

handle.h_mapsel_all     =   uicontrol(handle.control3,'Style','pushbutton','String','All','callback',@CT_callback_mapselect,'Position',[.025 .675 .45 .3]);
handle.h_mapsel_box     =   uicontrol(handle.control3,'Style','pushbutton','String','Box','callback',@CT_callback_mapselect,'Position',[.025 .35 .45 .3]);
handle.h_mapsel_circle  =   uicontrol(handle.control3,'Style','pushbutton','String','Circle','callback',@CT_callback_mapselect,'Position',[.025 .025 .45 .3]);
% location event
handle.h_location       =   uicontrol(handle.control3,'Style','pushbutton','String','Locate','callback',@CT_callback_location,'Position',[.525 .025 .45 .45]); 

%% Input/Output panel
handle.io_panel         =   uipanel('parent', handle.f1,'title', 'I/O', 'pos', [0.7 .05 .28 .2],'fontsize',para.panelftsize);
% read in new evlist
handle.h_evlist         =   uicontrol(handle.io_panel,'Style','edit','String',para.evlistname,'callback',@CT_callback_load_evlist_2,'Position',[.25 .75 .5 .25],'BackgroundColor',para.bcolor);
                            uicontrol(handle.io_panel,'callback',@CT_callback_load_evlist,'String','Load evlist','Position',[.25 .5 .5 .25]);
% read in new listname
handle.h_listname       =   uicontrol(handle.io_panel,'Style','edit','String',para.listname,'callback',@CT_callback_load_listname_2,'Position',[.25 .25 .5 .25],'BackgroundColor',para.bcolor);
                            uicontrol(handle.io_panel,'Style','pushbutton','String','Load listname','callback',@CT_callback_load_listname,'Position',[.25 .0 .5 .25]);
% save traces
                            uicontrol(handle.io_panel,'Style','pushbutton','String','Save(s)','callback',@CT_callback_save,'Position',[.0 .0 .25 .5]);
% save figure
                            uicontrol(handle.io_panel,'Style','pushbutton','String','Save Fig','callback',@CT_callback_savefig,'Position',[.0 .5 .25 .5]);
% delete event
                            uicontrol(handle.io_panel,'Style','pushbutton','String','Delete(Ctrl+d)','callback',@CT_callback_del_event,'Position',[.75 .0 .25 .5]);
% reset event
                            uicontrol(handle.io_panel,'Style','pushbutton','String','Reset (r)','callback',@CT_callback_reset_event,'Position',[.75 .5 .25 .5]);

%% Event and pages
                            uicontrol(handle.f1,'Style','pushbutton','String','1st','callback',@CT_callback_firstevent,'Position', [.15 .95 .05 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','pre_ev(b)', 'callback',@CT_callback_preevent,'Position', [.20 .95 .1 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','next_ev(n)', 'callback',@CT_callback_nextevent,'Position', [.5 .95 .1 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','last','callback',@CT_callback_lastevent,'Position', [.6 .95 .05 .04]);
% pages
                            uicontrol(handle.f1,'Style','pushbutton','String','1st','callback',@CT_callback_firstpage,'Position', [.15 .90 .05 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','pre_ls(,)','callback',@CT_callback_prepage,'Position', [.20 .90 .1 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','next_ls(.)','callback',@CT_callback_nextpage,'Position', [.5 .90 .1 .04]);
                            uicontrol(handle.f1,'Style','pushbutton','String','last','callback',@CT_callback_lastpage,'Position', [.6 .90 .05 .04]);
% ievent
                            uicontrol(handle.f1,'Style','text','String','ievent','Position',[.65 .97 .05 .02]);
handle.h_ievent         =   uicontrol(handle.f1,'Style','edit','String','1','callback',@CT_callback_ievent,'Position',[.65 .95 .05 .02],'BackgroundColor',para.bcolor);
% iframe
                            uicontrol(handle.f1,'Style','text','String','iframe','Position',[.65 .92 .05 .02]);
handle.h_iframe         =   uicontrol(handle.f1,'Style','edit','String','1','callback',@CT_callback_iframe,'Position',[.65 .90 .05 .02],'BackgroundColor',para.bcolor);

%% list box and delete trace
handle.h_listbox        =   uicontrol(handle.f1,'Style','listbox','Value',0,'callback',@CT_callback_listbox,'keypressfcn', @CT_callback_listkey,'Position',[.65 .3 .05 .6],'max',1000);

                            uicontrol(handle.f1,'Style','text','String','Quality','Position',[.7 .88 .04 .02]);
handle.h_quality        =   uicontrol(handle.f1,'Style','popup','String',para.quality_list,'Value',1,'callback',@CT_callback_quality,'Position',[.7 .85 .04 .03]);
                            uicontrol(handle.f1,'Style','text','String','Group','Position',[.7 .83 .04 .02]);
handle.h_group          =   uicontrol(handle.f1,'Style','popup','String',para.group_list,'Value',1,'callback',@CT_callback_group,'Position',[.7 .80 .04 .03]);
                            uicontrol(handle.f1,'Style','text','String','Q show','Position',[.7 .77 .04 .02]);
handle.h_qshow          =   uicontrol(handle.f1,'Style','listbox','String',para.quality_liststr,'Value',1:length(para.quality_liststr),'Position',[.7 .66 .04 .11],'max',100);
                            uicontrol(handle.f1,'Style','text','String','G show','Position',[.7 .63 .04 .02]);
handle.h_gshow          =   uicontrol(handle.f1,'Style','listbox','String',para.group_liststr,'Value',1:length(para.group_list),'Position',[.7 .52 .04 .11],'max',100);
%      e.h_mquality       =   uicontrol(handle.f1,'Style','popup','String',para.quality_list,'Value',length(para.quality_list)-1,'callback',@CT_callback_mquality,'Position',[.7 .7 .04 .03]);
                            uicontrol(handle.f1,'Style','pushbutton','String','Pick(p)','callback',@CT_callback_pickt,'Position',[.7 .46 .04 .05]);
handle.h_tlistbox       =   uicontrol(handle.f1,'Style','listbox','callback',@CT_callback_tlistbox, 'Position',[.7 .36 .04 .10],'max',100);
                            uicontrol(handle.f1,'Style','pushbutton','String','AddT','callback',@CT_callback_addt,'Position',[.7 .33 .04 .03]);
                            uicontrol(handle.f1,'Style','pushbutton','String','DelT','callback',@CT_callback_delt,'Position',[.7 .3 .04 .03]);

                            
%% Basic functions
    function CT_plotmap(h,dummy)
            
        if isfield(handle,'h_map_base');
            delete(handle.h_map_base);
            handle = rmfield(handle,'h_map_base');
        end  
        if isfield(handle,'h_map_st');
            delete(handle.h_map_st);
            handle = rmfield(handle,'h_map_st');
        end   
        if isfield(handle,'h_map_sel');
            delete(handle.h_map_sel);
            handle = rmfield(handle,'h_map_sel');
        end        
        if isfield(handle,'h_plot_loc')
            delete(handle.h_plot_loc);
            delete(handle.h_text_loc);
            handle = rmfield(handle,'h_plot_loc');
            handle = rmfield(handle,'h_text_loc');
        end    
        if isfield(handle,'h_plot_ref')
            delete(handle.h_plot_ref);
            handle = rmfield(handle,'h_plot_ref');
        end
        cla(handle.hax3,'reset');
        hold(handle.hax3, 'on');             
        
        handle.h_map_base = plot(handle.hax3,para.long,para.lat,'k-');
        for j = 1:length(tr)
            handle.h_map_st(j) = plot(handle.hax3,tr(j).stlo,tr(j).stla,'marker',para.marker_quality{tr(j).quality_ind},'markersize',para.markersize_st,'color',para.markercolor_st,'linewidth',para.markerlinewidth_st);
        end
                
        % change plot properties of stations in selected region
        if isfield(para,'idx')
            for j = 1:length(para.idx)
                set(handle.h_map_st(para.idx(j)),'color',para.markercolor_group{tr(para.idx(j)).group_ind},'marker',para.marker_quality{tr(para.idx(j)).quality_ind},'markersize',para.markersize_box,'linewidth',para.markerlinewidth_box);
            end
        end
        
        if ~isfield(para,'x0') || isempty(para.x0)
        else
            handle.h_plot_ref = plot(handle.hax3,para.x0,para.y0,'color',para.markercolor_refpoint,'marker',para.marker_refpoint,'markersize',para.markersize_refpoint,'linewidth',para.markerlinewidth_refpoint);
        end
        
        % axis lim
        if isfield(para,'slat')
            handle.h_map_sel = plot(handle.hax3,para.slon,para.slat,'color','k','linewidth',para.markerlinewidth_box);
        end
%         xmin = min([para.slon]);
%         xmax = max([para.slon]);
%         ymin = min([para.slat]);
%         ymax = max([para.slat]);
        xmin = min([tr.stlo]);
        xmax = max([tr.stlo]);
        ymin = min([tr.stla]);
        ymax = max([tr.stla]);
        
        dx = (xmax-xmin)/10;
        dy = (ymax-ymin)/10;
        if dx == 0
            dx = 1;
        end
        if dy == 0
            dy = 1;
        end
        
        xlim(handle.hax3,[xmin-dx xmax+dx]);
        ylim(handle.hax3,[ymin-dy ymax+dy]);
        set(handle.hax3,'clipping','on');
        
        [aa,bb,cc] = fileparts(para.events{para.ievent});
        eventname = regexprep([bb,cc],'_','\\_');
        title(handle.hax1,eventname);
        
        box(handle.hax3, 'on'); 
        box(handle.hax1, 'on'); 
        box(handle.hax2, 'on'); 
        uicontrol(handle.h_hot);
    end
    
    function CT_loaddata(h,dummy)
        
        mintime = [];
        maxtime = [];
        num = 0;        
        for i = 1:length(para.idx)
            if ~isfield(tr(para.idx(i)),'cmps_isexist') || isempty(tr(para.idx(i)).cmps_isexist) % first time load in
                for j = 1:length(para.cmps)
                    tr(para.idx(i)).cmps_isexist(j) = 0;
                    tr(para.idx(i)).cmps_isshow(j) = 0;
                    tr(para.idx(i)).data{j} = [];
                    tr(para.idx(i)).time{j} = [];
                    tr(para.idx(i)).delta{j} = [];
                    tr(para.idx(i)).A0{j} = [];
                                      
                    ind = find(strcmpi(para.cmps{j},tr(para.idx(i)).cmps)); % find the component from data
                    if isempty(ind) % component does not exist
                        continue;
                    end
                    ind = ind(1); % choose the first one if there are multiple matches
                    if strcmpi(tr(para.idx(i)).filename{ind},'EMPTY') || strcmpi(tr(para.idx(i)).filename{ind},'NAN') % empty component
                        continue;
                    end
                    tr(para.idx(i)).cmps_isexist(j) = 1;
                    if para.cmps_isshow(j) == 0 % do not show component
                        continue;
                    end  
                    
                    % load in data
                    [hd,data] = irdsac(fullfile(para.events{para.ievent},tr(para.idx(i)).filename{ind}));
                    
                    % bad trace
                    if isempty(hd) || isnan(hd.npts) || isnan(hd.b) || isnan(hd.delta) || isnan(hd.o) || max(isnan(data))==1 || (max(data)-min(data)) == 0
                        continue
                    end
                    
                    % detrend 
                    data = detrend(data);
                    
                    % normalization data
                    A0 = max(data)-min(data);
                    tr(para.idx(i)).cmps_isshow(j) = 1;
                    tr(para.idx(i)).A0{j} = A0;
                    tr(para.idx(i)).data{j} = data/A0;
                    tr(para.idx(i)).time{j} = reshape(linspace(hd.b,hd.e,hd.npts)-hd.o,[],1); % Reference time is o
                    tr(para.idx(i)).delta{j} = hd.delta;
                    
                    num = num + 1;
                    mintime(num) = hd.b-hd.o;
                    maxtime(num) = hd.e-hd.o;
                end
            else
                for j = 1:length(para.cmps)                    
                    if para.cmps_isshow(j) == 0 || tr(para.idx(i)).cmps_isexist(j) == 0  %do not show or component does not exist
                        tr(para.idx(i)).cmps_isshow(j) = 0;
                        continue;
                    end
                    if tr(para.idx(i)).cmps_isshow(j) == 1 % already loaded in
                        num = num + 1;
                        mintime(num) = tr(para.idx(i)).time{j}(1);
                        maxtime(num) = tr(para.idx(i)).time{j}(end);
                        continue;
                    end
                    ind = find(strcmpi(para.cmps(j),tr(para.idx(i)).cmps)); % find the component from data
                    if isempty(ind) % component does not exist
                        continue;
                    end
                    ind = ind(1); % choose the first one if there are multiple matches
                    if strcmpi(tr(para.idx(i)).filename{ind},'EMPTY') || strcmpi(tr(para.idx(i)).filename{ind},'NAN') % empty component
                        continue;
                    end
                    
                    % load in data
                    [hd,data] = irdsac(fullfile(para.events{para.ievent},tr(para.idx(i)).filename{ind}));
                    
                    % bad trace
                    if isempty(hd) || isnan(hd.npts) || isnan(hd.b) || isnan(hd.delta) || isnan(hd.o) || max(isnan(data))==1 || (max(data)-min(data)) == 0
                        continue
                    end
                    
                    % normalization data
                    A0 = max(data)-min(data);
                    tr(para.idx(i)).cmps_isshow(j) = 1;
                    tr(para.idx(i)).A0{j} = A0;
                    tr(para.idx(i)).data{j} = data/A0;
                    tr(para.idx(i)).time{j} = reshape(linspace(hd.b,hd.e,hd.npts)-hd.o,[],1); % Reference time is o
                    tr(para.idx(i)).delta{j} = hd.delta;
                    
                    num = num + 1;
                    mintime(num) = hd.b-hd.o;
                    maxtime(num) = hd.e-hd.o;
                end
            end
            %             end
        end
        if ~isfield(para,'timewin_full') || isempty(para.timewin_full)
            para.timewin_full = [min(mintime) max(maxtime)];
        else
            para.timewin_full = [min([mintime para.timewin_full(1)]) max([maxtime para.timewin_full(2)])];
        end
        
        if para.timewin(1)==0 && para.timewin(2)==0
            para.timewin = para.timewin_full;
        end
        
    end

    function CT_filtering1(h,dummy)
        
        tmpa = 'NAN';
        tmpb = 'NAN';
        if strcmpi(para.filterlist{para.filterlist_ind1},'bp')
            if para.fl1 > 0 && para.fh1 > para.fl1 % band pass filtering
                tmpa = para.fl1;
                tmpb = para.fh1;
            elseif para.fh1 > 0 && para.fl1 == 0 % low pass filtering
                tmpa = para.fh1;
                tmpb = 'low';
            elseif para.fl1 > 0 && para.fh1 == 0 % high pass filtering
                tmpa = para.fl1;
                tmpb = 'high';
            end
        elseif strcmpi(para.filterlist{para.filterlist_ind1},'lp')
            if para.fl1 >0 || para.fh1 >0
                tmpa = max([para.fl1,para.fh1]);
                tmpb = 'low';
            end 
        elseif strcmpi(para.filterlist{para.filterlist_ind1},'hp')
            if para.fl1 >0 && para.fh1 >0
                tmpa = min([para.fl1,para.fh1]);
                tmpb = 'high';
            elseif para.fl1>0 || para.fh1>0
                tmpa = max([para.fl1,para.fh1]);
                tmpb = 'high';
            end 
        end
        
        for i = 1:length(para.idx)
            for j = 1:length(tr(para.idx(i)).data)
                if isempty(tr(para.idx(i)).delta{j})
                    tr(para.idx(i)).data_f1{j} = [];
                    tr(para.idx(i)).A0_f1{j} = [];
                    continue;
                end
                if (isnumeric(tmpa) && tmpa > 1/tr(para.idx(i)).delta{j}/2) || (isnumeric(tmpb) && tmpb > 1/tr(para.idx(i)).delta{j}/2)
                    data_tmp = tr(para.idx(i)).data{j}; % exceeds Nyquist frequency -- no filtering
                else
                    data_tmp = filtering(tr(para.idx(i)).data{j},tr(para.idx(i)).delta{j},tmpa,tmpb,para.order1);
                end
                
                % is envelope
                if para.isenv == 1
                    data_tmp = hilbert(data_tmp);
                    data_tmp = abs(data_tmp);
                    if para.issemilogy == 1
                        data_tmp = real(log10(data_tmp));
                        data_tmp = detrend(data_tmp);
                    end
                    if para.issmooth == 1
                        nt = 2*round(para.smoothtime/tr(para.idx(i)).delta{j}/2) + 1;
                        data_tmp = smooth(data_tmp,nt);
                    end
                end             
                                
                % normalization
                tr(para.idx(i)).A0_f1{j} = max(data_tmp) - min(data_tmp);
                tr(para.idx(i)).data_f1{j} = data_tmp/tr(para.idx(i)).A0_f1{j};
            end
        end
        
    end

    function CT_sort(h,dummy)
        if isempty(tr)
            return;
        end
        para.sortind = get(handle.h_sort_list,'Value');
        if strcmpi(para.sortlist{para.sortind},'stla')
            offsettmp = [tr(para.idx).stla];
            para.x0 = [];
            para.y0 = [];
        elseif strcmpi(para.sortlist{para.sortind},'stlo')
            offsettmp = [tr(para.idx).stlo];
            para.x0 = [];
            para.y0 = [];
        elseif strcmpi(para.sortlist{para.sortind},'point')
            fprintf('Pick a point on the map. Stations are sorted by the distance to this point\n');
%             if ~isfield(para,'x0') || isempty(para.x0)
                axes(handle.hax3);
                [x1,y1,but1] = myginput(1,'crosshair');
                para.x0 = x1;
                para.y0 = y1;
%             end
            offsettmp = distance([tr(para.idx).stla],[tr(para.idx).stlo],para.y0,para.x0);
        else
            offsettmp = 1:length(para.idx);
        end
        
        [B,IX] = sort(offsettmp,'ascend');
        if ~isfield(para,'idx')
            return;
        end
        para.idx = para.idx(IX);
        para.idx_one = para.idx(1);
        para.iframe = 1;
        
        CT_callback_para2GUI(h,dummy);  
    end                        

%% Callback functions

	%% load, select map
    function CT_callback_loadmap(h,dummy)
        
        uicontrol(handle.h_loadmap);
        CT_callback_reset(h,dummy);
        
        if isfield(para,'evlist')
            para.events = textread(para.evlist,'%s');
            para.nevent = length(para.events);
        end                 
        
        if isfield(para,'tname')
            para = rmfield(para,'tname');
        end
        
        [ind, flag_file] = ct_find_event_list(para.ievent, para.events, para.listname, para.outlist,'forward');
        if isempty(ind)
            para.ievent = para.nevent;
            CT_callback_lastevent (h, dummy);
        else
            para.ievent = ind;
        end
        
        if flag_file == 2
            fid = fopen(fullfile(para.events{para.ievent},para.outlist),'r');
        elseif flag_file == 1
            fid = fopen(fullfile(para.events{para.ievent},para.listname),'r');
        else
        end
        C = textscan(fid,'%s',-1,'delimiter','\n','commentstyle','#');
        Ccell = C{1};
        fclose(fid);
        num = 0;
        tr = [];
        for i = 1:length(Ccell)
            Cstr = deblank(Ccell{i});
            parts = textscan(Cstr,'%s','Whitespace',' \b\t');
            partstr = parts{1};
            num = num + 1;
            tr(num).netwk = partstr{1};
            tr(num).stnm = partstr{2};
            tr(num).stname = [partstr{1},'.',partstr{2}];
            tr(num).stla = str2num(partstr{3});
            tr(num).stlo = str2num(partstr{4});
            tr(num).quality = partstr{5};
            tr(num).group = partstr{6};            
            ind = find(strcmpi(para.quality_list,tr(num).quality));
            if isempty(ind)
                tr(num).quality_ind = length(para.quality_list);
                tr(num).quality = para.quality_list{tr(num).quality_ind};
            else
                tr(num).quality_ind = ind;
                tr(num).quality = para.quality_list{tr(num).quality_ind};
            end
            ind = find(strcmpi(para.group_list,tr(num).group));
            if isempty(ind)
                tr(num).group_ind = length(para.group_list);
                tr(num).group = para.group_list{tr(num).group_ind};
            else
                tr(num).group_ind = ind;
                tr(num).group = para.group_list{tr(num).group_ind};
            end
                
            % load in components and corresponding filenames
            j = 1;
            while 1
                if strcmpi(partstr{7+2*j-1},'*')
                    break;
                else
                    tr(num).cmps{j}=partstr{7+2*j-1};
                    tr(num).filename{j}=partstr{7+2*j};
                    j = j + 1;
                end
            end
            % read in marked arrival times
            tr(num).nt = length(partstr)-(7+2*j-1);
            for k = 1:length(partstr)-(7+2*j-1)
                ttmp = str2num(lower(partstr{7+2*j-1+k}));
                if ~isempty(ttmp)
                    tr(num).t(k) = ttmp;
                else
                    tr(num).t(k) = nan;
                end
            end            
        end   
        
        % resize tr.t
        if ~isempty(tr)
            nt_max = max([tr.nt]);
            for i = 1:length(tr)
                tr(i).t(tr(i).nt+1:nt_max) = nan;
                tr(i).nt = nt_max;
            end
            for i = 1:nt_max
                para.tname{i} = ['T',num2str(i)];
            end
        end
        
        CT_plotmap(h,dummy); 
    end

    function CT_callback_mapselect(h,dummy)

        if ~isfield(para,'mapresel') || para.mapresel == 0 || ~isfield(para,'slat')
            if isempty(tr)
                return;
            end
            
            h_mapsel_all = get(handle.h_mapsel_all,'Value');
            h_mapsel_circle = get(handle.h_mapsel_circle,'Value');
%             h_mapsel_box = get(handle.h_mapsel_box,'Value');
            
            if h_mapsel_all==1
                x1 = min([tr.stlo]);
                x2 = max([tr.stlo]);
                y1 = min([tr.stla]);
                y2 = max([tr.stla]);
                para.slat = [y1 y1 y2 y2 y1];
                para.slon = [x1 x2 x2 x1 x1];
            else
                fprintf('Please select a region on map\n');
                
                [x1,y1,but1] = myginput(1,'crosshair');
                htmp = plot(x1,y1,'+');
                [x2,y2,but2] = myginput(1,'crosshair');
                delete(htmp);
                if h_mapsel_circle==1
                    [para.slat,para.slon] = scircle2(y1,x1,y2,x2);
                else
                    para.slat = [y1 y1 y2 y2 y1];
                    para.slon = [x1 x2 x2 x1 x1];
                end
            end
        end
        
        slat = para.slat;
        slon = para.slon;
        CT_callback_reset(h,dummy);
        para.slat = slat;
        para.slon = slon;
        
        para.mapresel = 0;
        IN = inpolygon([tr.stlo],[tr.stla],para.slon,para.slat);
        para.idx = find(IN);
        
        CT_callback_GUI2para(h,dummy);
        % keep stations with selected quality and group
        para.idx = para.idx(ismember([tr(para.idx).quality_ind],para.qshow_ind));
        para.idx = para.idx(ismember([tr(para.idx).group_ind],para.gshow_ind));
        
        CT_plotmap(h,dummy);
        
        if isempty(para.idx)
            return;
        end
        para.nframes = ceil(length(para.idx)/para.n_per_frame);
        para.iframe = 1;
        para.offset = 1:length(para.idx);
        para.idx_one = para.idx(1);
        
        if isfield(para,'idx_frame')
            para = rmfield(para,'idx_frame');
        end
        
        uicontrol(handle.h_hot);
    end

    %% locate events
    function CT_callback_location(h,dummy)
        uicontrol(handle.h_hot);
        [vmodel(:,1) vmodel(:,2)] = textread(para.vmodelfile,'%f %f','commentstyle','shell');
        input.vmodel = vmodel;
        input.xdep = para.xdep;
        
        % for selected times
        tstr = get(handle.h_tlistbox,'String');
        if isempty(tstr)
            return;
        end
        tind = get(handle.h_tlistbox,'Value');
        num = 0;
        opt = [];
        fprintf('Now locating ...\n');  
%         handle.hmsg = msgbox('Now locating ...','Message','modal'); 
        for i = 1:length(tind)
            ix = tind(i);
            r = cellfun(@(x) x(ix),{tr(para.idx).t},'uni',false);
            ind1 = find(~isnan([r{:}]));
            
            % for all groups
            groups = unique({tr(para.idx(ind1)).group});            
            
            if isempty(groups)
                continue;
            end
            
            for j = 1:length(groups)
                ind2 = find(strcmpi({tr(para.idx(ind1)).group},groups(j)));
                input.stla = [tr(para.idx(ind1(ind2))).stla]';
                input.stlo = [tr(para.idx(ind1(ind2))).stlo]';
                input.tt = [r{ind1(ind2)}]';
%                 input.tt = tr(para.idx(ind1(ind2))).t(tind(i));
                input.xlat = min(input.stla)-para.ostla:para.ddelta:max(input.stla)+para.ostla;
                input.xlon = min(input.stlo)-para.ostlo:para.ddelta:max(input.stlo)+para.ostlo;
                
                num = num + 1;
                opt(num).group = groups{j};
                opt(num).tname = para.tname{tind(i)};
                opt(num).isout = 1;
                opt(num).stla = input.stla;
                opt(num).stlo = input.stlo;
                opt(num).netwk = {tr(para.idx(ind1(ind2))).netwk}';
                opt(num).stnm = {tr(para.idx(ind1(ind2))).stnm}';
                if length(input.stla)<=2
                    continue;
                end
                [ loc ] = CT_location_v2( input );
                opt(num).loc = loc;
                if opt(num).loc.lat==input.xlat(1) || opt(num).loc.lat==input.xlat(2) || opt(num).loc.lon==input.xlon(1) || opt(num).loc.lon==input.xlon(2)
                else
                    opt(num).isout = 0;
                end
            end
        end
        
        if isfield(handle,'h_plot_loc')
            delete(handle.h_plot_loc);
            delete(handle.h_text_loc);
            handle = rmfield(handle,'h_plot_loc');
            handle = rmfield(handle,'h_text_loc');
        end
        
        fprintf('Plot location results\n');       
        axes(handle.hax3);
        % plot on map
        for i = 1:length(opt)
            num = num + 1;
            if opt(i).isout == 1 % out of box
                fprintf('No result for time %s group %s\n',opt(i).tname,opt(i).group);
                if ~isfield(opt,'loc') || isempty(opt(i).loc)
                    handle.h_plot_loc(num,1) = plot(handle.hax3,mean(opt(i).stlo),mean(opt(i).stla),'color',para.markercolor_noloc,'Marker','s','Markersize',para.markersize_loc/2,'linewidth',para.markerlinewidth_st/2);
                    handle.h_text_loc(num) = text(mean(opt(i).stlo),mean(opt(i).stla),[' ',opt(i).tname,'.',opt(i).group,': No result'],'color',para.markercolor_noloc,'Fontsize',para.tmark_fontsize/2);
                else
                    handle.h_plot_loc(num,1) = plot(handle.hax3,opt(i).loc.lon,opt(i).loc.lat,'color',para.markercolor_noloc,'Marker','s','Markersize',para.markersize_loc/2,'linewidth',para.markerlinewidth_st/2);
                    handle.h_text_loc(num) = text(opt(i).loc.lon,opt(i).loc.lat,[' ',opt(i).tname,'.',opt(i).group,': No result'],'color',para.markercolor_noloc,'Fontsize',para.tmark_fontsize/2);
                end
            else
                handle.h_plot_loc(num,1) = plot(handle.hax3,opt(i).loc.lon,opt(i).loc.lat,'color',para.markercolor_loc,'Marker',para.marker_loc,'Markersize',para.markersize_loc,'linewidth',para.markerlinewidth_st);
                handle.h_plot_loc(num,2) = plot(handle.hax3,[opt(i).loc.lon opt(i).loc.lon],[opt(i).loc.lat-opt(i).loc.lat_err opt(i).loc.lat+opt(i).loc.lat_err],'-','color',para.markercolor_loc,'Markersize',para.markersize_loc,'linewidth',para.markerlinewidth_st);
                handle.h_plot_loc(num,3) = plot(handle.hax3,[opt(i).loc.lon-opt(i).loc.lon_err opt(i).loc.lon+opt(i).loc.lon_err],[opt(i).loc.lat opt(i).loc.lat],'-','color',para.markercolor_loc,'Markersize',para.markersize_loc,'linewidth',para.markerlinewidth_st);
                handle.h_text_loc(num) = text(opt(i).loc.lon,opt(i).loc.lat,[' ',opt(i).tname,'.',opt(i).group],'color',para.tmark_color,'Fontsize',para.tmark_fontsize);
            end
        end
        
        % output location file
        fid = fopen(fullfile(para.events{para.ievent},para.outlocfile),'w');
        fprintf(fid,'#tname group nstation evla evlo evdp(fixed) lat_err lon_err in/out\n');
        for i = 1:length(opt)
            if isfield(opt,'loc') && ~isempty(opt(i).loc) 
                if opt(i).isout==0
                    fprintf(fid,'%s %s %d %f %f %f %f %f in\n',opt(i).tname, opt(i).group, length(opt(i).stla), opt(i).loc.lat, opt(i).loc.lon, opt(i).loc.dep, opt(i).loc.lat_err, opt(i).loc.lon_err);
                else
                    fprintf(fid,'%s %s %d %f %f %f %f %f out\n',opt(i).tname, opt(i).group, length(opt(i).stla), opt(i).loc.lat, opt(i).loc.lon, opt(i).loc.dep, opt(i).loc.lat_err, opt(i).loc.lon_err);
                end
            end
        end
        fclose(fid);
        fprintf('Locations are saved to %s\n',fullfile(para.events{para.ievent},para.outlocfile));
%         delete(handle.hmsg);
    end

    %% plot data
    function CT_callback_iniplot(h,dummy)
        
        CT_callback_GUI2para(h,dummy);
        para.mapresel = 1;
        CT_callback_mapselect(h,dummy);
        if ~isfield(para,'idx') || isempty(para.idx)
            return;
        end
        CT_loaddata(h,dummy);
        CT_filtering1(h,dummy);
        CT_sort(h,dummy);
        set(handle.h_listbox,'Value',1); 
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);  
        CT_plotmap(h,dummy);
    end

    function CT_callback_ZNERT(h,dummy)
        if ~isfield(para,'idx') || isempty(para.idx)
            return;
        end
        CT_callback_GUI2para(h,dummy);
        CT_loaddata(h,dummy);
        CT_filtering1(h,dummy);  
        set(handle.h_listbox,'Value',1);      
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_tracenumber(h,dummy)
        para.n_per_frame = floor(str2num(get(handle.h_ntrace_num,'String')));
        if ~isfield(para,'idx')
            uicontrol(handle.h_hot);
            return;
        end
        para.nframes = ceil(length(para.idx)/para.n_per_frame);
        para.iframe = 1;
        set(handle.h_iframe,'String',num2str(para.iframe));
        set(handle.h_iframe,'Value',1);
        fprintf('Number of trace per frame set to %d\n',para.n_per_frame);  
        set(handle.h_listbox,'Value',1); 
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_sort(h,dummy)
        if ~isfield(para,'idx') || isempty(para.idx)
            return;
        end
        set(handle.h_listbox,'Value',1);
        CT_sort(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);  
        CT_plotmap(h,dummy);
    end

    function CT_callback_filter1(h,dummy)
        if ~isfield(para,'idx') || isempty(para.idx)
            return;
        end
        CT_callback_GUI2para(h,dummy);
        CT_filtering1(h,dummy);        
        CT_callback_plottrace1(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_replot(h,dummy)
        CT_callback_GUI2para(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_plottrace1(h,dummy)
        
        % no data
        if ~isfield(tr,'data') || diff(para.timewin) == 0 
            return;
        end
        if isfield(handle,'h_trace_plot1');
            for kk = 1:length(handle.h_trace_plot1)
                delete(handle.h_trace_plot1{kk});
            end
            handle = rmfield(handle,'h_trace_plot1');
        end
        if isfield(handle,'h_plot_tmark');
            delete(handle.h_plot_tmark);
            handle = rmfield(handle,'h_plot_tmark');
            delete(handle.h_text_tmark);
            handle = rmfield(handle,'h_text_tmark');
        end
        cla(handle.hax1,'reset');
        hold(handle.hax1,'on');
              
%         CT_callback_GUI2para(h,dummy);              
        ind_tmp = (para.iframe-1)*para.n_per_frame+1:min([para.iframe*para.n_per_frame length(para.idx)]);
        para.idx_frame = para.idx(ind_tmp);
        para.offset_frame = para.offset(ind_tmp);
        
        C = intersect(para.idx_one,para.idx_frame);
        if isempty(C)
            para.idx_one = para.idx_frame(1);
        end
                
%         fprintf('Plot traces... ');          
%         handle.hmsg = msgbox('Plotting ...','Message','modal'); 
        for i = 1:length(para.idx_frame)
            num = 0;
            for j = 1:para.ncmps_show
                ind = para.cmps_showind(j);
                if tr(para.idx_frame(i)).cmps_isshow(ind) == 0
                    continue;
                end
                time_tmp = para.timewin(1):para.delta:para.timewin(2);
%                 time_tmp = para.timewin(1):tr(para.idx_frame(i)).delta{ind}:para.timewin(2);
                data_tmp = interp1(tr(para.idx_frame(i)).time{ind},tr(para.idx_frame(i)).data_f1{ind},time_tmp,'linear',0);
                            
                % remove mean
                data_tmp = detrend(data_tmp,'constant');
                
                % normalize
                A0_tmp = max(data_tmp)-min(data_tmp);
                data_tmp = data_tmp/A0_tmp;
                
                num = num+1;
                handle.h_trace_plot1{i}(num) = plot(handle.hax1, time_tmp, -para.scale*data_tmp/para.ncmps_show + para.offset_frame(i) + j/(para.ncmps_show+1),'color',para.cmps_colors_show{j},'linewidth',para.linewidth_show);            
            end
        end
        
        set(handle.h_listbox, 'String', {tr(para.idx_frame).stname});
        if isfield(para,'tname')
            set(handle.h_tlistbox,'String', para.tname);
        end
        set(handle.hax1,'Ytick',min(para.offset_frame):1:max(para.offset_frame));
        box(handle.hax1,'on');  
        grid(handle.hax1,'on');               
        xlim(handle.hax1,para.timewin);
        ylim(handle.hax1,[min(para.offset_frame) max(para.offset_frame)+1]);
        ylabel(handle.hax1,'#');
        set(handle.hax1,'Ydir','reverse');   
        [aa,bb,cc] = fileparts(para.events{para.ievent});
        eventname = regexprep([bb,cc],'_','\\_');
        title(handle.hax1,eventname);
        
        CT_callback_plottmark(h,dummy);
    end
               
    function CT_callback_plottrace2(h,dummy)
        
        % no data
        if ~isfield(tr,'data') || diff(para.timewin) == 0 
            return;
        end
        
        if isfield(handle,'h_trace_plot2');
            delete(handle.h_trace_plot2);
            handle = rmfield(handle,'h_trace_plot2');
        end
        cla(handle.hax2,'reset');
        hold(handle.hax2,'on');
        axes(handle.hax2);
        
        % choose one station to show in panel2
        if ~isfield(para,'idx_one')
            para.idx_one = para.idx_frame(1);
        end
        
        para.fl2 = str2num(get(handle.h_filter_fl2,'String'));
        para.fh2 = str2num(get(handle.h_filter_fh2,'String'));
        para.order2 = str2num(get(handle.h_filter_order2,'String'));  
        if isempty(para.fl2),	para.fl2 = 0;	end
        if isempty(para.fh2),	para.fh2 = 0;	end
        if isempty(para.order2),	para.order2 = 2;	end
        
        para.filterlist_ind2 = get(handle.h_filter_list2,'Value');     
        
        tmpa = 'NAN';
        tmpb = 'NAN';
        if strcmpi(para.filterlist{para.filterlist_ind2},'bp')
            if para.fl2 > 0 && para.fh2 > para.fl2 % band pass filtering
                tmpa = para.fl2;
                tmpb = para.fh2;
            elseif para.fh2 > 0 && para.fl2 == 0 % low pass filtering
                tmpa = para.fh2;
                tmpb = 'low';
            elseif para.fl2 > 0 && para.fh2 == 0 % high pass filtering
                tmpa = para.fl2;
                tmpb = 'high';
            end
        elseif strcmpi(para.filterlist{para.filterlist_ind2},'lp')
            if para.fl2 >0 || para.fh2 >0
                tmpa = max([para.fl2,para.fh2]);
                tmpb = 'low';
            end 
        elseif strcmpi(para.filterlist{para.filterlist_ind2},'hp')
            if para.fl2 >0 && para.fh2 >0
                tmpa = min([para.fl2,para.fh2]);
                tmpb = 'high';
            elseif para.fl2>0 || para.fh2>0
                tmpa = max([para.fl2,para.fh2]);
                tmpb = 'high';
            end 
        end
        
        para.isspectrogram2 = get(handle.h_spectrogram2,'Value');
        para.spec_twinl = str2num(get(handle.h_spec_twinl,'String'));
        para.spec_tstep = str2num(get(handle.h_spec_tstep,'String'));
        window = round(para.spec_twinl/para.delta);
        noverlap = round((para.spec_twinl-para.spec_tstep)/para.delta);
        P_stack = [];
        time_tmp = para.timewin(1):para.delta:para.timewin(2);
        if para.isspectrogram2 == 1 && length(time_tmp) >= window
            for j = 1:para.ncmps_show
                ind = para.cmps_showind(j);
                if tr(para.idx_one).cmps_isshow(ind) == 0
                    continue;
                end
                
                % filtering
                data_tmp = filtering(tr(para.idx_one).data{ind},tr(para.idx_one).delta{ind},tmpa,tmpb,para.order2);
                % cut
                data_tmp = interp1(tr(para.idx_one).time{ind},data_tmp,time_tmp,'linear',0);
                
                [S,F,T,P] = spectrogram(data_tmp,window,noverlap,[],1/para.delta);
                if isempty(P_stack)
                    P_stack = P;
                else
                    P_stack = P_stack + P;
                end
                                
            end  
            handle.h_trace_plot2(1) = imagesc(T+time_tmp(1),F,10*log10(abs(P_stack)));    
        
            x1 = get(handle.hax2,'Position');
            handle.h_trace_plot2(2) = colorbar;
            ylabel(handle.h_trace_plot2(2),'dB');
            x = get(handle.h_trace_plot2(2),'Position');
            dcx = 0.005;
            x(3) = 2*dcx;
            x(1) = x1(1) + x1(3) + dcx/2;
            set(handle.h_trace_plot2(2),'Position',x);
            
            ylim(handle.hax2,[F(1) F(end)]);
            ylabel(handle.hax2,'Frequency (Hz)');
        else
            if para.isspectrogram2 == 1 && length(time_tmp) < window
                fprintf('Data length must be longer than window length\n');
                para.isspectrogram2 = 0;
                set(handle.h_spectrogram2,'Value',0);
            end
            num = 0;
            for j = 1:para.ncmps_show
                ind = para.cmps_showind(j);
                if tr(para.idx_one).cmps_isshow(ind) == 0
                    continue;
                end
                
                % filtering
                data_tmp = filtering(tr(para.idx_one).data{ind},tr(para.idx_one).delta{ind},tmpa,tmpb,para.order2);
                
                % cut
                data_tmp = interp1(tr(para.idx_one).time{ind},data_tmp,time_tmp,'linear',0);
                
                A0_tmp = max(data_tmp)-min(data_tmp);
                data_tmp = data_tmp/A0_tmp;
                
                num = num + 1;
                handle.h_trace_plot2(num) = plot(handle.hax2, time_tmp, -para.scale*data_tmp/para.ncmps_show + 1 + j/(para.ncmps_show+1),'color',para.cmps_colors_show{j},'linewidth',para.linewidth_show);
                
            end
            set(handle.hax2,'Ytick',[1 2]);
            ylim(handle.hax2,[1 2]);
            ylabel(handle.hax2,'#');
            set(handle.hax2,'Ydir','reverse');
        end            
        
        xlim(handle.hax2,para.timewin);
        xlabel(handle.hax2,'Time (s)');
        title(handle.hax2,[tr(para.idx_one).stname,'  (Press "Enter" to change)']);
        box(handle.hax2,'on'); 
        grid(handle.hax2,'on');  
    end

    function CT_callback_plottmark(h,dummy)
       
        if isfield(handle,'h_plot_tmark')
            delete(handle.h_plot_tmark);
            handle = rmfield(handle,'h_plot_tmark');
            delete(handle.h_text_tmark);
            handle = rmfield(handle,'h_text_tmark');
        end
        
        if ~isfield(para,'tname') || ~isfield(para,'idx_frame')
            return;
        end
        ind = get(handle.h_tlistbox,'Value');
        
        hold(handle.hax1,'on');
        axes(handle.hax1);
        num = 0;
        for i = 1:length(para.idx_frame)
            for j = 1:length(ind)
                if isnan(tr(para.idx_frame(i)).t(ind(j)))
                    continue;
                end
                if tr(para.idx_frame(i)).t(ind(j)) < para.timewin(1) || tr(para.idx_frame(i)).t(ind(j)) > para.timewin(2) % out of plot window
                    continue;
                end
                num = num + 1;
                handle.h_plot_tmark(num) = plot(handle.hax1,[tr(para.idx_frame(i)).t(ind(j)) tr(para.idx_frame(i)).t(ind(j))],para.offset_frame(i)+[0 1],'color',para.tmark_color,'linewidth',para.tmark_linewidth);
                handle.h_text_tmark(num) = text(tr(para.idx_frame(i)).t(ind(j)),para.offset_frame(i)+0.1,['T',num2str(ind(j))],'color',para.tmark_color,'Fontsize',para.tmark_fontsize);
            end
        end
                
    end

    function CT_callback_windowzoom(h,dummy)
        
        if ~isfield(handle,'h_trace_plot1')
            return;
        end
        
        % get time axis
        timewin_L = str2num(get(handle.h_timewin_L,'String'));
        timewin_R = str2num(get(handle.h_timewin_R,'String'));

        % change the figure focus        
        hObject = handle.h_timewin_zoom; 
        uicontrol(hObject); 
        set(hObject, 'enable', 'off');
        drawnow;
        set(hObject, 'enable', 'on');

        while 1
            [x1,y1,but1] = myginput(1,'crosshair');
            if but1 == 2 || strcmpi(char(but1),'q') || but1 == 27
                break;
            elseif but1 == 3 || strcmpi(char(but1),'o')
                x1 = timewin_L;
                x2 = timewin_R;
            else
                [x2,y2,but2] = myginput(1,'crosshair');
                if but2 == 2 || strcmpi(char(but2),'q') || but1 == 27
                    break;
                elseif but2 == 3 || strcmpi(char(but2),'o')
                    x1 = timewin_L;
                    x2 = timewin_R;                    
                end
            end
            
            para.timewin = [min([x1 x2]) max([x1 x2])];
            set(handle.h_timewin_L,'String',num2str(para.timewin(1)));
            set(handle.h_timewin_R,'String',num2str(para.timewin(2)));            
            CT_callback_plottrace1(h,dummy);
            CT_callback_plottrace2(h,dummy);            
        end                
        
    end

    function CT_callback_windowfull(h,dummy)
        if ~isfield(para,'timewin_full')
            return;
        end
        para.timewin = para.timewin_full;
        set(handle.h_timewin_L,'String',num2str(para.timewin(1)));
        set(handle.h_timewin_R,'String',num2str(para.timewin(2)));        
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy); 
        uicontrol(handle.h_hot);
    end
    
    function CT_callback_ampup(h,dummy)
        para.scale = para.scale*1.25;     
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);  
        uicontrol(handle.h_hot);      
    end

    function CT_callback_ampdown(h,dummy)
        para.scale = para.scale*0.8;     
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);  
        uicontrol(handle.h_hot);      
    end

    function CT_callback_envelope(h,dummy)
        uicontrol(handle.h_hot)
        if ~isfield(tr,'data')
            return;
        end
        CT_callback_GUI2para(h,dummy);
        CT_filtering1(h,dummy);
        CT_callback_plottrace1(h,dummy);  
    end

    function CT_callback_inputpara(h,dummy)
        dlg_title = 'Plotting parameters';
        prompt = {'Smooth timelength'};
        def = {num2str(para.smoothtime)};
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        para.smoothtime = str2num(answer{1});
        CT_callback_plottrace1(h,dummy);
    end

%% event and frame navigation
    function CT_callback_firstevent(h,dummy)   
        uicontrol(handle.h_hot);
        if ~isfield(para,'nevent')
            return;
        end
        para.ievent = 1;
        [ind, flag_file] = ct_find_event_list(para.ievent, para.events, para.listname, para.outlist,'forward');
        if isempty(ind)
            return;
        else
            para.ievent = ind;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_loadmap(h,dummy);        
    end

    function CT_callback_preevent(h,dummy)   
        uicontrol(handle.h_hot);
        if ~isfield(para,'nevent')
            return;
        end
        if para.ievent ~= 1
            [ind, flag_file] = ct_find_event_list(para.ievent-1, para.events, para.listname, para.outlist,'backward');
        else
            [ind, flag_file] = ct_find_event_list(para.ievent, para.events, para.listname, para.outlist,'backward');
        end

        if isempty(ind)
            fprintf('This is the first event!\n');
            CT_callback_firstevent(h,dummy);
        else
            para.ievent = ind;
            CT_callback_para2GUI(h,dummy);
            CT_callback_loadmap(h,dummy);
        end
    end

    function CT_callback_nextevent(h,dummy)   
        uicontrol(handle.h_hot); 
        if ~isfield(para,'nevent')
            return;
        end
        if para.ievent ~= para.nevent
            [ind, flag_file] = ct_find_event_list(para.ievent+1, para.events, para.listname, para.outlist,'forward');
        else
            [ind, flag_file] = ct_find_event_list(para.ievent, para.events, para.listname, para.outlist,'forward');
        end

        if isempty(ind)
            fprintf('This is the last event!\n');
            CT_callback_lastevent(h,dummy);
        else
            para.ievent = ind;
            CT_callback_para2GUI(h,dummy);
            CT_callback_loadmap(h,dummy);
        end
    end

    function CT_callback_lastevent(h,dummy) 
        uicontrol(handle.h_hot);
        if ~isfield(para,'nevent')
            return;
        end
        para.ievent = para.nevent;
        [ind, flag_file] = ct_find_event_list(para.ievent, para.events, para.listname, para.outlist,'backward');
        if isempty(ind)
            return;
        else
            para.ievent = ind;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_loadmap(h,dummy);
    end

    function CT_callback_firstpage(h,dummy)
        uicontrol(handle.h_hot);
        set(handle.h_listbox,'Value',1);
        para.iframe = 1;
        if diff(para.timewin)==0
            para.iframe = 1;
            return;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
    end

    function CT_callback_prepage(h,dummy)
        uicontrol(handle.h_hot);
        set(handle.h_listbox,'Value',1);
        para.iframe = para.iframe - 1;
        if para.iframe < 1 ,
            para.iframe = 1;
        end
        if para.iframe > para.nframes,
            para.iframe = para.nframes;
        end
        if diff(para.timewin)==0
            para.iframe = 1;
            return;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
    end

    function CT_callback_nextpage(h,dummy)
        uicontrol(handle.h_hot);
        set(handle.h_listbox,'Value',1);
        para.iframe = para.iframe + 1;
        if para.iframe < 1 ,
            para.iframe = 1;
        end
        if para.iframe > para.nframes,
            para.iframe = para.nframes;
        end
        if diff(para.timewin)==0
            para.iframe = 1;
            return;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);        
    end

    function CT_callback_lastpage(h,dummy)
        uicontrol(handle.h_hot);
        set(handle.h_listbox,'Value',1);
        para.iframe = para.nframes;
        if diff(para.timewin)==0
            para.iframe = 1;
            return;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
    end

    function CT_callback_ievent(h,dummy)
        uicontrol(handle.h_hot);
        if ~isfield(para,'nevent')
            set(handle.h_ievent,'String','1');
            return;
        end
        para.ievent = floor(str2num(get(handle.h_ievent,'String')));
        if para.ievent <= 0            
            para.ievent = 1;
        elseif para.ievent > para.nevent
            para.ievent = para.nevent;
        end
        CT_callback_loadmap(h,dummy);
    end

    function CT_callback_iframe (h, dummy)
        uicontrol(handle.h_hot);
        set(handle.h_listbox,'Value',1);
        if ~isfield(para,'nevent')
            set(handle.h_iframe,'String','1');
            return;
        end
        para.iframe = floor(str2num(get(handle.h_iframe,'String')));
        if para.iframe <= 0            
            para.iframe = 1;
        elseif para.iframe > para.nframes
            para.iframe = para.nframes;
        end
        if diff(para.timewin)==0
            para.iframe = 1;
            return;
        end
        CT_callback_para2GUI(h,dummy);
        CT_callback_plottrace1(h,dummy);
        CT_callback_plottrace2(h,dummy);
    end


%% I/O
    function CT_callback_load_evlist(h,dummy)
                
        [templist, temppath] = uigetfile({'*.txt','Txt Files (*.txt)';'*.m;*.fig;*.mat;*.mdl','Matlab Files (*.m,*.fig,*.mat,*.mdl)';'*.*', 'All Files (*.*)'}, 'Load event list',para.evlistname);
        if templist == 0
        else
            para.evpathname = temppath;
            para.evlistname = templist;
            set(handle.h_evlist,'String', para.evlistname);
            para.evlist = fullfile(para.evpathname,para.evlistname);
            para.ievent = 1;            
            CT_callback_reset(h,dummy);
            CT_callback_loadmap(h,dummy);
        end        
        uicontrol(handle.h_hot);
    end

    function CT_callback_load_evlist_2(h,dummy)
        para.evlistname = get(handle.h_evlist,'String');
        set(handle.h_evlist,'String', para.evlistname);
        para.evlist = fullfile(para.evpathname,para.evlistname);
        para.ievent = 1;        
        CT_callback_reset(h,dummy);
        CT_callback_loadmap(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_load_listname(h,dummy)
                
        if ~isfield(para,'nevent')
            return;
        end
        [templist, temppath] = uigetfile({'*list*', 'List Files (*list*)';'*.txt','Txt Files (*.txt)';'*.m;*.fig;*.mat;*.mdl','Matlab Files (*.m,*.fig,*.mat,*.mdl)';'*.*', 'All Files (*.*)'}, 'Load event list',para.events{para.ievent});
        
        if templist ~= 0
            if strcmp(temppath(end),'/') || strcmp(temppath(end),'\')
                temppath = temppath(1:end-1);
            end
            [temp1,temp2,temp3] = fileparts(temppath);
            evname = [temp2,temp3];
            
            ind = strfind(para.events,evname);
            tempev = find(~cellfun(@isempty,ind));
            if ~isempty(tempev)
                if strcmpi(templist(end-1:end),para.outlistid);
                    para.listname = templist(1:end-2);
                    para.outlistname = templist;
                else
                    para.listname = templist;
                    para.outlistname = [templist,para.outlistid];
                end                    
                set(handle.h_listname,'String', para.listname);
                para.ievent = tempev(1);            % change event number
                
                CT_callback_reset(h,dummy);
                CT_callback_loadmap(h,dummy);
            end
        end
        uicontrol(handle.h_hot)
    end

    function CT_callback_load_listname_2 (h, dummy)
                
        templist = get(handle.h_listname,'String');
        if strcmpi(templist(end-1:end),para.outlistid);
            para.listname = templist(1:end-2);
            para.outlistname = templist;
        else
            para.listname = templist;
            para.outlistname = [templist,para.outlistid];
        end
        
        CT_callback_reset(h,dummy);
        CT_callback_loadmap(h,dummy);        
        uicontrol(handle.h_hot);
    end

    function CT_callback_reset_event(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        if exist(fullfile(para.events{para.ievent},para.outlist),'file')
            delete(fullfile(para.events{para.ievent},para.outlist));
        end
        fprintf('Reset Event: %s\n',para.events{para.ievent});
        CT_callback_loadmap(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_del_event(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        fid = fopen(fullfile(para.events{para.ievent},para.outlist),'w');
        fclose(fid);
        CT_callback_nextevent(h,dummy);
        uicontrol(handle.h_hot);
    end

    function CT_callback_save(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        if ~isfield(tr,'t')
            ind = [];
        else
            r = isnan(reshape([tr.t],[],length(tr))');
            ind = find(~min(r,[],1));
        end
        
        fid = fopen(fullfile(para.events{para.ievent},para.outlist),'w');
        fprintf(fid,'#Netwk Stnm Stla Stlo Quality Group * Cmp1 File1 Cmp2 File2 ... *\n');
        for j = 1:length(tr)
            fprintf(fid,'%s %s %f %f %s %s *',tr(j).netwk,tr(j).stnm,tr(j).stla,tr(j).stlo,tr(j).quality,tr(j).group);
            for k = 1:length(tr(j).cmps)
                fprintf(fid,' %s %s',tr(j).cmps{k},tr(j).filename{k});
            end
            fprintf(fid,' *');
            
            for k = 1:length(ind)
                fprintf(fid,' %f',tr(j).t(ind(k)));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);  
        fprintf('Output %s saved to %s\n',para.outlist,para.events{para.ievent});
        uicontrol(handle.h_hot);      
    end

    function CT_callback_savefig(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        set(handle.f1,'PaperPositionMode','auto');        
        ptmp = get(handle.f1,'PaperPosition');
        set(handle.f1,'Papersize',[ ptmp(3) ptmp(4)]);
        
        [a,b,c] = fileparts(para.events{para.ievent});
        eventname_tmp = [b,c];
        figname_base = ['Fig_CT_',eventname_tmp,'_'];
        tmp = dir(fullfile(para.events{para.ievent},[figname_base,'*.pdf']));
        if isempty(tmp)
            nnum = 1;
        else
            fignames = sort({tmp.name});
            figend = fignames{end};
            nnum = str2num(figend(length(figname_base)+1:end-4)) + 1;
            if isempty(nnum)
                nnum = 1;
            end
        end
        print(handle.f1,'-r300','-dpdf',fullfile(para.events{para.ievent},[figname_base,num2str(nnum,'%02d'),'.pdf']));
        % print(handle.f1,'-painters','-fillpage','-r300','-dpdf',fullfile(para.events{para.ievent},[figname_base,num2str(nnum,'%02d'),'.pdf']));
        fprintf('Save figure to %s\n',fullfile(para.events{para.ievent},[figname_base,num2str(nnum,'%02d'),'.pdf']));
       
        uicontrol(handle.h_hot);        
    end

%% Listbox
    function CT_callback_listbox(h,dummy)              
        sind = get(handle.h_listbox,'Value');
        if isempty(sind)
            return;
        end
        
        % set quality and group - show the first selected trace
        set(handle.h_quality,'Value',tr(para.idx_frame(sind(1))).quality_ind);
        set(handle.h_group,'Value',tr(para.idx_frame(sind(1))).group_ind);
        
        for k = 1:length(para.idx_frame)
            ind = find(tr(para.idx_frame(k)).cmps_isshow);
            for m = 1:length(handle.h_trace_plot1{k})
                set( handle.h_trace_plot1{k}(m),'linewidth',para.linewidth_show,'color',para.cmps_colors{ind(m)});
            end
        end
        for k = 1:length(sind)
            ind = find(tr(para.idx_frame(sind(k))).cmps_isshow);
            for m = 1:length(handle.h_trace_plot1{sind(k)})
                set( handle.h_trace_plot1{sind(k)}(m),'linewidth',para.linewidth_selected,'color',para.cmps_colors_sel{ind(m)});
            end            
        end
        
        % for stations in map
        for k = 1:length(para.idx)
            set(handle.h_map_st(para.idx(k)),'markersize',para.markersize_box,'linewidth',para.markerlinewidth_box);
        end
        for k = 1:length(sind)
            set(handle.h_map_st(para.idx_frame(sind(k))),'markersize',para.markersize_sel,'linewidth',para.markerlinewidth_sel);
        end
    end

    function CT_callback_listkey(src,evnt)        
        key = evnt.Key;
        sind = get(handle.h_listbox,'Value');
        tind = get(handle.h_tlistbox,'Value');        
        if strcmpi(key,'delete')
            for k = 1: length(sind)
                tr(para.idx_frame(sind(k))).quality = para.quality_list{end}; 
                tr(para.idx_frame(sind(k))).quality_ind = length(para.quality_list); 
                set(handle.h_quality,'Value',length(para.quality_list));
                set(handle.h_trace_plot1{sind(k)},'visible','off');
                set(handle.h_map_st(para.idx_frame(sind(k))),'marker',para.marker_quality{tr(para.idx_frame(sind(k))).quality_ind});        
            end
        elseif strcmpi(key,'backspace')
            for k = 1: length(sind)
                if strcmpi(tr(para.idx_frame(sind(k))).quality,para.quality_list{end});
                    tr(para.idx_frame(sind(k))).quality = para.quality_list{1};
                    tr(para.idx_frame(sind(k))).quality_ind = 1;
                    set(handle.h_quality,'Value',1);
                    set(handle.h_trace_plot1{sind(k)},'visible','on');
                    set(handle.h_map_st(para.idx_frame(sind(k))),'marker',para.marker_quality{1});
                end
            end
        elseif strcmpi(key,'return')
            para.idx_one = para.idx_frame(sind(1));
            CT_callback_plottrace2(src,evnt);
            uicontrol(handle.h_listbox);            
        elseif strcmpi(key, 'escape') || strcmpi(key, 'alt')
            uicontrol(handle.h_hot);
            return;
        else
            if ~isfield(para,'tname') || max(tind)>length(para.tname)
                return;
            end
                        
            if ~isempty(evnt.Modifier) && strcmpi(evnt.Modifier{:},'control') && strcmpi(key, 'leftarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) - 10*para.delta;
                end
            elseif ~isempty(evnt.Modifier) && strcmpi(evnt.Modifier{:},'control') && strcmpi(key, 'rightarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) + 10*para.delta;
                end
            elseif ~isempty(evnt.Modifier) && strcmpi(evnt.Modifier{:},'shift') && strcmpi(key, 'leftarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) - 100*para.delta;
                end
            elseif ~isempty(evnt.Modifier) && strcmpi(evnt.Modifier{:},'shift') && strcmpi(key, 'rightarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) + 100*para.delta;
                end                   
            elseif (isempty(evnt.Modifier)  || ~strcmpi(evnt.Modifier{:},'control')) && strcmpi(key, 'leftarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) - para.delta;
                end
            elseif (isempty(evnt.Modifier) || ~strcmpi(evnt.Modifier{:},'control')) && strcmpi(key, 'rightarrow')
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = tr(para.idx_frame(sind(k))).t(tind) + para.delta;
                end
%             elseif strcmpi(key,'0') || strcmpi(key,'hyphen') % delete selected time picks
            elseif strcmpi(key,'hyphen') % delete selected time picks
                for k = 1:length(sind)
                    tr(para.idx_frame(sind(k))).t(tind) = nan;
                end
            else
                return;
            end
            CT_callback_plottmark(src,evnt);
            uicontrol(handle.h_listbox);
        end
    end

    function CT_callback_quality(h,dummy)
        uicontrol(handle.h_hot);
        ind = get(handle.h_quality,'Value');
        sind = get(handle.h_listbox,'Value');
        if ~isfield(para,'idx_frame')
            return;
        end
        for k = 1:length(sind)
            tr(para.idx_frame(sind(k))).quality = para.quality_list{ind};
            tr(para.idx_frame(sind(k))).quality_ind = ind;
            if ind == length(para.quality_list)
                set(handle.h_trace_plot1{sind(k)},'visible','off');
            else
                set(handle.h_trace_plot1{sind(k)},'visible','on');
            end
            set(handle.h_map_st(para.idx_frame(sind(k))),'marker',para.marker_quality{tr(para.idx_frame(sind(k))).quality_ind});
        end
    end

    function CT_callback_group(h,dummy)
        uicontrol(handle.h_hot);
        ind = get(handle.h_group,'Value');
        sind = get(handle.h_listbox,'Value');
        if ~isfield(para,'idx_frame')
            return;
        end
        for k = 1:length(sind)
            tr(para.idx_frame(sind(k))).group = para.group_list{ind};
            tr(para.idx_frame(sind(k))).group_ind = ind;
            set(handle.h_map_st(para.idx_frame(sind(k))),'color',para.markercolor_group{tr(para.idx_frame(sind(k))).group_ind});
        end
    end
    
%% Arrival time picking
    function CT_callback_tlistbox(h,dummy)
        
        sind = get(handle.h_tlistbox,'Value');
        if isempty(sind)
            return;
        end
          
        CT_callback_plottmark(h,dummy);
        
    end

    function CT_callback_pickt(h,dummy)
        
        if ~isfield(para,'tname') || isempty(para.tname) || ~isfield(handle,'h_trace_plot1')
            uicontrol(handle.h_hot);
            return;
        end
        
        tind = get(handle.h_tlistbox,'Value');
        
        while 1
            [x,y,but] = myginput(1,'crosshair');
            if but ~= 1 && but ~= 3
                break;
            end
            y = floor(y);
            ind = find(para.offset_frame - y == 0);
            if isempty(ind)
                continue;
            end
            
            if but == 1
                tr(para.idx_frame(ind)).t(tind(1)) = x;
            elseif but == 3
                tr(para.idx_frame(ind)).t(tind(1)) = nan;
            end            
            CT_callback_plottmark(h,dummy);
        end
                    
    end

    function CT_callback_addt(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        if isfield(para,'tname')
            para.tname{length(para.tname)+1} = ['T',num2str(length(para.tname)+1)];
        else
            para.tname = {'T1'};
        end
        
        for i = 1:length(tr)
            tr(i).nt = length(para.tname);
            tr(i).t(tr(i).nt) = nan;
        end        
        set(handle.h_tlistbox,'String',para.tname);
        set(handle.h_tlistbox,'Value',length(para.tname));
        CT_callback_plottmark(h,dummy);
    end

    function CT_callback_delt(h,dummy)
        if ~isfield(para,'idx_frame')
            return;
        end
        tind = get(handle.h_tlistbox,'Value');   
        for i = 1:length(para.idx_frame)
            tr(para.idx_frame(i)).t(tind) = nan;
        end
        list_entry = get(handle.h_tlistbox,'String');
        nl = length(list_entry);
        nlkeep = setdiff(1:nl,tind);
        if isempty(nlkeep)
            tr = rmfield(tr,'t');
        else
            for i = 1:length(tr)
                tr(i).t = tr(i).t(nlkeep);
            end
        end
        for i = 1:length(tr)
            tr(i).nt = length(nlkeep);
        end
        para = rmfield(para,'tname');
        for i = 1:length(nlkeep)
            para.tname{i} = ['T',num2str(i)];
        end
        if isfield(para,'tname')
            set(handle.h_tlistbox,'String',para.tname);
            set(handle.h_tlistbox,'Value',min([tind(1) length(para.tname)]));
        else
            set(handle.h_tlistbox,'String',[]);  
            set(handle.h_tlistbox,'Value',1);          
        end
        CT_callback_plottmark(h,dummy);
    end
  
%% set, clear, reset, change figure focus
    % set parameters to GUI
    function CT_callback_para2GUI(h,dummy)
        if isfield(para,'n_per_frame')
            set(handle.h_ntrace_num,'String',num2str(para.n_per_frame));
        end        
        if isfield(para,'fl1')
            set(handle.h_filter_fl1,'String',num2str(para.fl1));
        end   
        if isfield(para,'fh1')
            set(handle.h_filter_fh1,'String',num2str(para.fh1));
        end 
        if isfield(para,'order1')
            set(handle.h_filter_order1,'String',num2str(para.order1));
        end
        if isfield(para,'fl2')
            set(handle.h_filter_fl2,'String',num2str(para.fl2));
        end   
        if isfield(para,'fh2')
            set(handle.h_filter_fh2,'String',num2str(para.fh2));
        end 
        if isfield(para,'order2')
            set(handle.h_filter_order2,'String',num2str(para.order2));
        end
        if isfield(para,'timewin')
            set(handle.h_timewin_L,'String',num2str(para.timewin(1)));
            set(handle.h_timewin_R,'String',num2str(para.timewin(2)));
        end    
        if isfield(para,'spec_twinl')
            set(handle.h_spec_twinl,'String',num2str(para.spec_twinl));
        end  
        if isfield(para,'spec_tstep')
            set(handle.h_spec_tstep,'String',num2str(para.spec_tstep));
        end
        if isfield(para,'ievent')
            set(handle.h_ievent,'String',num2str(para.ievent));
        end
        if isfield(para,'iframe')
            set(handle.h_iframe,'String',num2str(para.iframe));  
        end
        if isfield(para,'idx_frame')
            set(handle.h_quality,'Value',tr(para.idx_frame(1)).quality_ind);
            set(handle.h_group,'Value',tr(para.idx_frame(1)).group_ind);   
        end
    end

    % get parameters from GUI
    function CT_callback_GUI2para(h,dummy)
        h_cmp_Z = get(handle.h_cmp_Z,'Value');
        h_cmp_N = get(handle.h_cmp_N,'Value');
        h_cmp_E = get(handle.h_cmp_E,'Value');
        h_cmp_R = get(handle.h_cmp_R,'Value');
        h_cmp_T = get(handle.h_cmp_T,'Value');
        para.cmps_isshow = [h_cmp_Z h_cmp_N h_cmp_E h_cmp_R h_cmp_T];
        para.cmps_showind = find(para.cmps_isshow);
        para.ncmps_show = length(para.cmps_showind);
        para.cmps_colors_show = para.cmps_colors(para.cmps_showind);
        para.n_per_frame = str2num(get(handle.h_ntrace_num,'String'));
        if isempty(para.n_per_frame),   para.n_per_frame = 10;  end
        para.fl1 = str2num(get(handle.h_filter_fl1,'String'));
        para.fh1 = str2num(get(handle.h_filter_fh1,'String'));
        para.order1 = str2num(get(handle.h_filter_order1,'String'));
        para.filterlist_ind1 = get(handle.h_filter_list1,'Value');
        if isempty(para.fl1),	para.fl1 = 0;	end
        if isempty(para.fh1),	para.fh1 = 0;	end
        if isempty(para.order1),	para.order1 = 2;	end
        para.fl2 = str2num(get(handle.h_filter_fl2,'String'));
        para.fh2 = str2num(get(handle.h_filter_fh2,'String'));
        para.order2 = str2num(get(handle.h_filter_order2,'String'));
        para.filterlist_ind2 = get(handle.h_filter_list2,'Value');
        if isempty(para.fl2),	para.fl2 = 0;	end
        if isempty(para.fh2),	para.fh2 = 0;	end
        if isempty(para.order2),	para.order2 = 2;	end
        para.delta = str2num(get(handle.h_delta,'String'));
        para.timewin(1) = str2num(get(handle.h_timewin_L,'String'));
        para.timewin(2) = str2num(get(handle.h_timewin_R,'String')); 
        para.timewin = sort(para.timewin);
        para.isenv = get(handle.h_envelope,'Value');    
        para.issemilogy = get(handle.h_semilogy,'Value');  
        para.issmooth = get(handle.h_smooth,'Value'); 
        para.ievent = str2num(get(handle.h_ievent,'String'));
        para.iframe = str2num(get(handle.h_iframe,'String'));
        para.qshow_ind = get(handle.h_qshow,'Value');
        para.gshow_ind = get(handle.h_gshow,'Value');
    end
    
    % reset handles, parameters, and tr
    function CT_callback_reset(h,dummy)

        if isfield(handle,'h_map_base')
            delete(handle.h_map_base);
            handle = rmfield(handle,'h_map_base');
        end
        if isfield(handle,'h_map_st')
            delete(handle.h_map_st);
            handle = rmfield(handle,'h_map_st');
        end
        if isfield(handle,'h_map_sel')
            delete(handle.h_map_sel);
            handle = rmfield(handle,'h_map_sel');
        end
        if isfield(handle,'h_trace_plot1')
            for kk = 1:length(handle.h_trace_plot1)
                delete(handle.h_trace_plot1{kk});
            end
            handle = rmfield(handle,'h_trace_plot1');
        end
        if isfield(handle,'h_trace_plot2')
            delete(handle.h_trace_plot2);
            handle = rmfield(handle,'h_trace_plot2');
        end
        if isfield(handle,'h_plot_tmark')
            delete(handle.h_plot_tmark);
            handle = rmfield(handle,'h_plot_tmark');
            delete(handle.h_text_tmark);
            handle = rmfield(handle,'h_text_tmark');
        end
        if isfield(handle,'h_plot_loc')
            delete(handle.h_plot_loc);
            delete(handle.h_text_loc);
            handle = rmfield(handle,'h_plot_loc');
            handle = rmfield(handle,'h_text_loc');
        end
        if isfield(handle,'h_plot_ref')
            delete(handle.h_plot_ref);
            handle = rmfield(handle,'h_plot_ref');
        end
        cla(handle.hax1,'reset');
        cla(handle.hax2,'reset'); 
        
        set(handle.h_listbox,'String',[]);
        set(handle.h_listbox,'Value',1);
        set(handle.h_tlistbox,'String',[]);
        set(handle.h_tlistbox,'Value',1);
        set(handle.h_quality,'Value',1);
        set(handle.h_group,'Value',1);
        set(handle.h_sort_list,'Value',1);
        
        % clear variables
        if isfield(para,'idx')
            para = rmfield(para,'idx');
        end
        if isfield(para,'idx_one')
            para = rmfield(para,'idx_one');
        end
        if isfield(para,'idx_frame')
            para = rmfield(para,'idx_frame');
        end
        if isfield(para,'offset')
            para = rmfield(para,'offset');
        end
        if isfield(para,'offset_frame')
            para = rmfield(para,'offset_frame');
        end
        if isfield(para,'timewin_full')
            para = rmfield(para,'timewin_full');
        end
        if isfield(para,'slat')
            para = rmfield(para,'slat');
            para = rmfield(para,'slon');
        end
        if isfield(para,'x0')
            para = rmfield(para,'x0');
            para = rmfield(para,'y0');
        end
                    
        para.iframe = 1;
        para.nframes = 1;
        para.timewin = [0 0];
        
        CT_callback_para2GUI(h,dummy);
    end

%% Hot keys
    function CT_short_cut(src, evnt)
        if strcmpi(evnt.Key,'m')
            CT_callback_loadmap(src,evnt)
        elseif strcmpi(evnt.Key,'i')
            CT_callback_iniplot(src,evnt)
        elseif strcmpi(evnt.Key,'b') || strcmpi(evnt.Key,'backspace')
            CT_callback_preevent(src,evnt)
        elseif strcmpi(evnt.Key,'n') || strcmpi(evnt.Key,'space')
            CT_callback_nextevent(src,evnt)
        elseif strcmpi(evnt.Key,'comma')
            CT_callback_prepage(src,evnt)
        elseif strcmpi(evnt.Key,'period')
            CT_callback_nextpage(src,evnt)            
        elseif strcmpi(evnt.Key,'s')
            CT_callback_save(src,evnt)
        elseif ~isempty(evnt.Modifier) && strcmp(evnt.Modifier{:},'control') && strcmpi(evnt.Key,'d')
            CT_callback_del_event(src,evnt)
        elseif strcmpi(evnt.Key,'r')
            CT_callback_reset_event(src,evnt)
        elseif strcmpi(evnt.Key,'p') || strcmpi(evnt.Key,'t')
            CT_callback_pickt(src,evnt)
        elseif strcmpi(evnt.Key,'e')
            CT_callback_envelope(src,evnt)
        elseif strcmpi(evnt.Key,'z')
            CT_callback_windowzoom(src,evnt);
        elseif strcmpi(evnt.Key,'f')
            CT_callback_windowfull(src,evnt);
        elseif strcmpi(evnt.Key,'equal')
            CT_callback_ampup(src,evnt)
        elseif strcmpi(evnt.Key,'hyphen')
            CT_callback_ampdown(src,evnt)
        elseif strcmpi(evnt.Key,'alt') || (~isempty(evnt.Modifier) && ( strcmpi(evnt.Key,'0') && strcmpi(evnt.Modifier{:},'command')))
            uicontrol(handle.h_listbox)
        else
        end
    end

    function CT_copyright(h,dummy)
        
        msgbox('Copyright @Chunquan Yu; Reference: Chao, K. and C. Yu, A MATLAB GUI-based Package for Examining Triggered Tremor: A Case Study in New Zealand (2018), Seismol. Res. Lett., submitted.','Message','modal')
    end
end

