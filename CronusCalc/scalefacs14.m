%
% sf=scalefacs14(sample,scaling_model)
%
function sf=scalefacs14(sp14,scaling_model)
%
% Check if scaling_model was specified, otherwise set it to default
%
if (~exist('scaling_model','var')), scaling_model = 'all'; end

%
% Setup the scale factors.
%
sf.ST=sp14.ST;
sf.SLth=1;
sf.SLeth=1;
sf.P=sp14.P;
sf.elevation=sp14.elevation;

%
% Use Greg/Nat's code to compute the time dependent scale factors.
% The reason that we have a separate tdsfsample here is that the
% get_tdsf() function wasn't written by us.
%

load pmag_consts
tdsfsample.lat=sp14.latitude;
tdsfsample.long=sp14.longitude;
tdsfsample.pressure=sp14.P;
tdsfsample.elevation=sp14.elevation;
tdsfsample.scaling=scaling_model;
sf.tdsf=get_tdsf(tdsfsample,pmag_consts);
