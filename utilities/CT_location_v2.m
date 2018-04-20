function [ loc ] = CT_location_v2( input )
% [ loc ] = CT_location( 'picked_time.txt' )

% v2: modified by Chunquan Yu, @ MIT-EAPS, Sep. 2015
% change input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: The sub-routine ' makeTTgrid1', 'trace_ray', 'delaz', and 'coortr'
%        are part of the CORAL Tools by Dr. Kenneth C. Creager
%   http://earthweb.ess.washington.edu/creager/
%   https://www.ess.washington.edu/dwp/people/profile.php?name=creager--ken
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

staLat = input.stla;
staLon = input.stlo;
TT = input.tt;
model = input.vmodel;

% search grid
opt.grd.lat= input.xlat;
opt.grd.lon= input.xlon;
opt.grd.depth = input.xdep;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% % Initial inout 1: arrival time
% %%%%%%%%%%%%%%%%%%%%%%%%
% % input='picked_time.txt';
% [stn_name,staLon,staLat,TT]=textread(input,'%s %f %f %f');

% calculate dT of all station-pair
num_stn=length(staLon); % number of staions
pa=combnk(1:num_stn,2); % list all possible pairs

% %%%%%%%%%%%%%%%%%%%%%%%%
% % Initial inout 2: searching range, velocity model
% %%%%%%%%%%%%%%%%%%%%%%%%
% lon_min = 132.1;
% lon_max = 133.6;
% lat_min =  32.7;
% lat_max =  34.1;
% dep     =  35.0; % km
% 
% opt.grd.depth = [dep];
% opt.grd.lon=[ lon_min : 0.01 : lon_max ];
% opt.grd.lat=[ lat_min : 0.01 : lat_max ];


% opt.model=0;
% opt.model=[
% 2.20    0.0
% 2.93    2.0
% 3.24    4.0
% 3.31    6.0
% 3.31    9.0
% 3.42    13.0
% 3.46    17.0
% 3.50    21.0
% 3.65    25.0
% 3.78    30.0
% 4.02    35.0
% 4.60    50.0
% 4.82    70.0
% ];
%[opt.TTgrid,opt.LAT,opt.LON,opt.DEP] = makeTTgrid1([staLat],[staLon],opt.grd,opt.model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make predicted differential travel-time grid
% [opt.TTgrid,opt.LAT,opt.LON,opt.DEP] = makeTTgrid1([staLat],[staLon],opt.grd);
[opt.TTgrid,opt.LAT,opt.LON,opt.DEP] = makeTTgrid1(staLat,staLon,opt.grd,model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  caculate RMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RMS =  sqrt [ ( (theo1-obs1)^2 + (theo2-obs2)^2 +.... (theo_n-obs_n)^2 )/n ] %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate td for all station pairs
num_pair=length(pa);
for p = 1:num_pair
    
        sta_a(p) = pa(p,1);
        sta_b(p) = pa(p,2);

        % theo time-different between 1-2
        aa = opt.TTgrid(1,sta_a(p)).TIM; % dt of 1 in pair 1-2
        bb = opt.TTgrid(1,sta_b(p)).TIM; % dt of 2 in pair 1-2
        dd(1,p).TTgrid.TIM = bb-aa;
        
        % observed time-different between 1-2
        TT_a(p) = TT(pa(p,1));
        TT_b(p) = TT(pa(p,2));
        dd_obs = TT_b(p) - TT_a(p); % obs tD
        
	% (theo-obs)^2
	sqr(p).theo_obs = times( (dd(1,p).TTgrid.TIM - dd_obs), (dd(1,p).TTgrid.TIM - dd_obs) );
        
end

% sum all (theo-obs)^2
sumtheo_obs = 0;
for q = 1:num_pair
    sumtheo_obs = sumtheo_obs + sqr(q).theo_obs; 
end

% take suqare root
RMS=sqrt( sumtheo_obs/num_pair );

%%%%%%%%%%%%%%%%%%%  find minimum RMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[out_row out_col] = find(RMS==min(min(RMS)));

lon = opt.LON(out_row,out_col); % find lon
lat = opt.LAT(out_row,out_col);
dep = opt.DEP(out_row,out_col);
RMS_min=min(min(RMS));

% calculate Chi-square
% see Peter Shearer, p.130-131
[lon_err,lat_err] = loc_err(opt,RMS,out_row,out_col,num_pair);

loc.lon     = lon;
loc.lat     = lat;
loc.dep     = dep;
loc.lon_err = lon_err;
loc.lat_err = lat_err;

%%%%%%%%%%%%%%%%%%% End of the Main Location Script %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sub-routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TTgrid,LAT,LON,DEP] = makeTTgrid1(staLat, staLon, grd, model);
%  makeTTGrid      calculate travel time grid
% USAGE: [TTgrid,LAT,LON,DEP] = makeTTgrid(staLat, staLon, grd, model);
%
% Input:
%  staLat   vector of station latitudes (geographic degrees)
%  staLon   vector of station longitudes (geographic degrees)
%  grd      structure containing 3-D source grid
%           grd.lon = vector of source longitudes (geographic degrees)
%           grd.lat = vector of source latitudes (geographic degrees)
%           grd.depth = vector of source depths (km positive down)
%
% Output:
%  TTgrid   vector of structures for each station containing containing fields
%           TTgrid(ksta).TIM  3-D array of travel time (s)
%           TTgrid(ksta).XI   2-D array of azimuths (deg)
%           TTgrid(ksta).DELT 3-D array of horizontal distancd (km)
%  LAT      3-D array of latitudes (geographic degrees)
%  LON      3-D array of longitudes (geographic degrees)
%  DEP      3-D array of depths (km)
% modified 7/20/07 to optionaly accept a new earth model
% staLat=[D.staLat]; staLon=[D.staLon]; grd=opt.grd; model=opt.model;

if nargin<4; model=''; end; % use default (PNSN P2 S-wave model) if no model is passed in

if min(grd.depth)<1
	d_depth=1-min(grd.depth);
	grd.depth=grd.depth+d_depth;
	disp(sprintf('Warning! Top of grid cannot be < 1 km. Shifting grid down %.1f km...',d_depth));
end

[LON,LAT,DEP]          = meshgrid(grd.lon,grd.lat,grd.depth);  % latitude and longitude of grid points
LAT2  = squeeze(LAT(:,:,1));
LON2  = squeeze(LON(:,:,1));
[ngrd1,ngrd2,ngrd3]= size(LON);                  % dimensions of the grids
ngrd               = ngrd1*ngrd2*ngrd3;          % total number of grid points

for kdepth=1:ngrd3;   % loop over depths

  [Xray,Tray,Pray]=trace_ray(grd.depth(kdepth),model);  % calculate distance, time, ray parameter for 'model' and source depth.

  for ksta=1:length(staLat);  % loop over stations

    [DELTS,AZIMS] = delaz(staLat(ksta),staLon(ksta),LAT2(:),LON2(:),0);  % distance and azimuth from station to each grid point
    DELTS=DELTS*111.1;                           % station to source epicentral distance (km) for lat, lon each grid point
    
    tmp  = interp1(Xray,Tray,DELTS);
    TIM  = reshape(tmp,  ngrd1,ngrd2);           % travel time (s) to center of array from each grid point
    XI   = reshape(AZIMS,ngrd1,ngrd2);           % station to source azimuth (deg) at each lat,lon grid point
    DELT = reshape(DELTS,ngrd1,ngrd2);           % station to source azimuth (deg) at each grid point
    
    if kdepth==1;
      TTgrid(ksta).XI   = XI;
      TTgrid(ksta).DELT = DELT;
    end
    TTgrid(ksta).TIM(:,:,kdepth)=TIM;
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,T,P,vs,vr]=trace_ray(depth, model);
%   trace_ray    calculate travel time curve by tracing rays
% USAGE: [X,T,P,vs,vr]=trace_ray(depth, model);
%
%  Input: 
%  depth is earthquake depth (km)
%  model is an optional velocity model (Nx2 array) with
%  velocity (km/s) in first column and depth (km) in second
%  velocity is a layered model, each line of the model gives the
%  velocity of a layer and depth at the top of the layer
%  use default (PNSN Puget Sound (P2)) model if model='' or if 
%  there is only one input argument
% 
% Output:
%  X horizontal distance (km)
%  T travel time (s)
%  P ray parameter (dT/dX s/km)
%  vs velocity at source (km/s)
%  vr velocity at receiver (surface) (km/s)
%
% example:
% clear; figure(20);clf; depth=[25 35 45 55]; colr = {'k' 'm' 'r' 'b'} ;
% for k=1:4; [X,T,P]=trace_ray(depth(k));hold on;h(k)=plot(X,T-X/4.5033,['.-' colr{k}]);end;axis([0 300 3 13]); grid on;
% legend(h,'25 km','35 km','45 km','55 km');xlabel('Range (km)'); ylabel('Travel Time - Range/4.5033 (s)')

%clear;close all;  % clear matlab workspace and close all figures

h0=depth;   % source depth

if nargin<2 || length(model)==0;  % use default model (PNSN model P2 converted to a shear wave model)

  vel_P2=[
   5.40  0.0
   6.38  4.0
   6.59  9.0
   6.73 16.0
   6.86 20.0
   6.95 25.0
   7.80 51.0;
   8.00 81.0];
 
  model = [vel_P2(:,1)/sqrt(3) ,  vel_P2(:,2)];
  
end

h  = model(:,2);    % column vector of depths corresponding to the tops of layers (km)
vl = model(:,1);    % velocity within each layer (km/s)

% add a new layer boundary at the source depth without changing the model
if length(find(h0==h))==0;
  i=max(find(h0>=h));
  h = [h(1:i); h0; h(i+1:end)];
  vl= vl([1:i i:end]);
else
  i=find(h0==h);
end
i0=i;
v0= vl(i0);
vs=v0;
vr=vl(1); %velocity at surface (receiver).

%   take_off_layer{1}  defines rays going up
%   take_off_layer{k}  (k>1) defines rays going down and reflecting critically off the layer boundary that is k-1 interfaces below the source

k0=1;
take_off_layer{k0} = -[0.01:5:80 81:2:87 88 89 89.9]*pi/180; % column vector of desired take-off angles (radians) < 0 means up-going ray

for k=i0+1:length(vl)-1;   % define rays reflecting critically off each interface below the source
  k0=k0+1;
  pmin = 1/vl(k+1);
  pmax = 1/vl(k);
  prange = pmin+(pmax-pmin)*[0:.05:.95 .96:.01:.99 .993:.002:.999];  
  prange = pmin+(pmax-pmin)*[.9:-.1:.1 .01];    
  take_off_layer{k0} = asin(v0*prange);
end

% combine the take-off angles into on vector 
take_off=[];
for k=1:1; % end loop at 1 for upgoing rays only and at k0 for all rays
  take_off=[take_off take_off_layer{k}]; 
end

k=find(take_off>0);
if length(k)>0;
  take_off(k) = max(take_off(k),asin(v0/vl(end))+eps );   % force all rays to turn with the model.
end


%h=[0:5.5:6310]';       % column vector of depths in the earth (km)
%r=6371-h;              % vector of equivalent radii
%v=prem(r);             % column vector of P-wave velocities (km/s) at depths h
%vl=(v(1:end-1)+v(2:end))/2; % mean velocity within each layer
dh=diff(h);                  % thickness of each layer (km)

%take_off=[30:1:85]'*pi/180;% column vector of desired take-off angles (radians)
%take_off=[26:.1:32]'*pi/180% column vector of desired take-off angles (radians)
%figure(1);clf;hold on     % open a new plot window
for l=1:length(take_off)   % loop over each ray
  clear dz dx dt x t i interface; % clear vectors containing ray increments
  k=1;
  interface(k) = find(h0==h);
  i(k)         = take_off(l);      % start a new ray
  x(k)         = 0;
  t(k)         = 0;
  while ~(interface(k)==1 & i(k)<0) % stop if interface is 1 and ray is going up.
    k=k+1;
    if i(k-1)<0;   % if ray going up then layer=interface-1
      layer = interface(k-1)-1; 
      layer_next = layer-1;
      interface(k) = interface(k-1)-1;
    else         
      layer = interface(k-1); % if ray going down that layer=interface
      layer_next = layer+1;
      interface(k) = interface(k-1)+1;
    end
    dz = h(interface(k)) - h(interface(k-1));
    dx = dz*tan(i(k-1)); % determine horizontal distance traveled in this layer
    dt = abs(dz/(vl(layer)*cos(i(k-1)))); % calculte travel time through this layer
    x(k)=x(k-1)+dx;
    t(k)=t(k-1)+dt;
    if layer_next == 0; break;end
    sini  = vl(layer_next)*sin(i(k-1))/vl(layer);% solve Snell's Law for ray angle in next layer
    if sini>1;                  % ray must reflect, change sign or ray angle 
      i(k) = -i(k-1);
    else
      i(k) = asin(sini);        % calculate ray angle in next layer
    end
    if k>20; break;end
  end
  %plot(x,-h(interface),'.-');axis('equal');hold on
  T(l,1)=t(end);
  X(l,1)=x(end);
end
P=sin(abs(take_off(:)))/v0; % horizontal slowness
%figure(2);subplot(211);plot(X,T,'.-r');axis([0 600 0 100]); subplot(212);plot(X,P,'.-r');axis([0 600 0 .3])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta,azeqst,azsteq]=delaz(eqlat,eqlon,stlat,stlon,flag);
%   delaz         compute earthquake/station distance and azimuth
% usage: [delta,azeqst,azsteq]=delaz(eqlat,eqlon,stlat,stlon,flag);
%
%     compute distance and azimuth from earthquake (eq) to station (st)
%     delta  = distance between (eq) and (st) in degrees
%     azeqst = azimuth from (eq) to (st) clockwise from north in degrees
%     azsteq = azimuth from (st) to (eq) clockwise from north in degrees
%
%     if input coordinates are geographic degrees   flag=0
%     if input coordinates are geocentric radians   flag=1
%
%     input latitudes and longitudes can be scalars or vectors
%     acceptable combinations are one earthquake with many stations,
%     one station and many earthquakes, or the same number of earthquakes
%     and stations
%     output vectors will have same dimensions as input vectors
%
%     calls coortr.m

% convert from geographic degrees to geocentric radians if necessary
% convert to spherical polar coordinates in radians (lat -> colatitude)

if flag==0,   % convert geographic degrees to geocentric radians
  [eqlat,eqlon]=coortr(eqlat,eqlon,flag); 
  [stlat,stlon]=coortr(stlat,stlon,flag); 
end
eqcolat=pi/2-eqlat;
stcolat=pi/2-stlat;

cos_eq=cos(eqcolat);
sin_eq=sin(eqcolat);
cos_st=cos(stcolat);
sin_st=sin(stcolat);
cos_eqst=cos(stlon-eqlon);
sin_eqst=sin(stlon-eqlon);

cos_delta=cos_eq.*cos_st + sin_eq.*sin_st.*cos_eqst;
sin_delta=real(sqrt(1-cos_delta.*cos_delta));
delta=atan2(sin_delta,cos_delta);
% if sin(delta)=0, set sin(delta)=eps=10**-16
sin_delta = sin_delta + (sin_delta==0)*eps;

% index is zero if expression is false, 1 if true; 
% if false, leave unchanged, if true azeqst=pi-azeqst
% this puts azeqst into the correct quadrant
azeqst=asin(sin_st.*sin_eqst./sin_delta);
index=(sin_eq.*cos_st - cos_eq.*sin_st.*cos_eqst < 0);
azeqst=azeqst + index.*(pi-2*azeqst);
azeqst=azeqst + (azeqst<0)*2*pi;

azsteq=asin(-sin_eq.*sin_eqst./sin_delta);
index=(cos_eq.*sin_st - sin_eq.*cos_st.*cos_eqst < 0);
azsteq=azsteq + index.*(pi-2*azsteq);
azsteq=azsteq + (azsteq<0)*2*pi;

% convert to degrees
delta=delta*180/pi;
azeqst=azeqst*180/pi;
azsteq=azsteq*180/pi;

% force answers to be real to fix problem at delta=180 degrees
delta=real(delta);
azeqst=real(azeqst);
azsteq=real(azsteq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [latout,lonout]=coortr(latin,lonin,flag);
%   coortr        geocentric/geographic coordinate transformation
% usage: [latout,lonout]=coortr(latin,lonin,flag);
%
%  purpose: transform between geographic and geocentric coordinates
%           geographic degrees to geocentric radians ( if flag=0 )
%        or geocentric radians to geographic degrees ( if flag=1 )
%           earthquake and station locations are typically
%           given in geographic coordinates, whereas most calculations
%           such as epicentral distance are given in geocentric coordinates
%           latin and lonin can be scalars or vectors, 
%           latout and lonout will match the dimensions of latin and lonin
%
%  if flag==0,
%     latin,  lonin  are latitude and longitude in geographic degrees
%     latout, lonout are latitude and longitude in geocentric radians
%  if flag==1,
%     latin,  lonin  are latitude and longitude in geocentric radians
%     latout, lonout are latitude and longitude in geographic degrees

if flag==0,
  latout=atan(tan(latin*pi/180)*0.9933056);
  lonout=lonin*pi/180;
elseif flag==1,
  latout=atan(tan(latin)/0.9933056)*180/pi;
  lonout=lonin*180/pi;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lon_err,lat_err] = loc_err(opt,RMS,out_row,out_col,num_pair);

% calculate Chi-square
% see Peter Shearer, p.130.
if num_pair<5
        ndf=num_pair;
else
        ndf=num_pair-4; % degree-of-freedom = n -4; n travel times are used;
end

%%%%%%%%%%%%%%%%%
% eq.5.30, p.130
%%%%%%%%%%%%%%%%%
% sum all variance
sumvar_best = (RMS(out_row,out_col)^2)*num_pair;
var_best = sumvar_best/ndf;

%sumvar_best = 0;
%for pp = 1:num_pair
%    sumvar_best = sumvar_best + sqr(pp).theo_obs(out_row,out_col);
%end
%var_best = sumvar_best/ndf;


%%%%%%%%%%%%%%%%%
% eq.5.31, p.131
%%%%%%%%%%%%%%%%%
sumvar_all = (RMS.^2)*num_pair;
chisqr = sumvar_all/var_best;

%sumvar_all = 0;
%for qq = 1:num_pair
%    sumvar_all = sumvar_all + sqr(qq).theo_obs;
%end
%chisqr = sumvar_all/var_best;

%%%%%%%%%%%%%%%%%%
% find 95% and 68% value for Chi-Sqr distribution
%RMS_chi=chi2pdf(RMS,ndf);
%chi_95=chi2inv(0.95,ndf); %chi95 = (ChiSqr<chi_95);
chi_68=chi2inv(0.68,ndf);

[err_row err_col] = find(chisqr<chi_68);
max_r=max(err_row);
max_c=max(err_col);
lon_max=opt.LON(max_r,max_c);
lat_max=opt.LAT(max_r,max_c);

min_r=min(err_row);
min_c=min(err_col);
lon_min=opt.LON(min_r,min_c);
lat_min=opt.LAT(min_r,min_c);

lon_err=(lon_max-lon_min)/2;
lat_err=(lat_max-lat_min)/2;
%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
