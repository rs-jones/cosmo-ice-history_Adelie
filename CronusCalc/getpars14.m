%
%  [pp,sp14,sf14,cp14]=getpars14(sampledata14, scaling_model, maxdepth);  
%
%
% The sampledata14 vector contains the following information:
%
%1. Latitude (decimal degrees)
%2. Longitude (decimal degrees)
%3. Elevation (meters)
%4. Pressure (hPa)       
%5. sample thickness (cm)
%6. bulk density (g/cm^3)
%7. Shielding factor for terrain, snow, etc. (unitless)
%8. Erosion-rate epsilon (mm/kyr)
%9. Sample 10-Be concentration (atoms of 10-Be/g of target)
%10. Inheritance for Be (atoms 10-Be/g of target)
%11. Lambda effective (attenuation length in g/cm2)
%12. Lambdafe Effective neutron attenuation length (g/cm^2)
%13. Depth to top of sample (g/cm^2)
%
% The maxdepth parameter is optional if it is given, it specifies
% the maximum depth (in g/cm^2) for which production rates will be
% computed.  If not given, the default value (currently 2500
% g/cm^2) is used instead. 
%
function [pp,sp14,sf14,cp14]=getpars14(sampledata14, scaling_model, maxdepth)
%
% If maxdepth is not specified, use a default.
%
if (nargin < 3)
  maxdepth=2500;
end
%
% Get the physical parameters.
%
pp=physpars();
%
% Extract the sample parameters from the sampledatavector.
%
sp14=samppars14(sampledata14);
%
% Get the scale factors.
%
sf14=scalefacs14(sp14);
%
% Computed parameters.
%
cp14=comppars14(pp, sp14, sf14, maxdepth);
%
% Go ahead and produce contemporary scaling factors.
%
sf14.currentsf=getcurrentsf(sf14, 0, scaling_model, 'c');


