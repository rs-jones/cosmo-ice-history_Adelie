%
% N14=predN14depth(pp,sp,sf,cp,age,depths,scaling_model)
%
% Predicts the 14-C concentration for a depth sample.
%
%  pp,sp,sf,cp            physics, site, sample, and
%                         computed parameters.
%  age                    Exposure age in kyr.
%  depths                 Depths in g/cm^2.
%  deltat                 Time step in ka (deltat=0.1 is 100 yrs)
%
%
function N14=predN14depth(pp,sp,sf,cp,age,depths,scaling_model,deltat)
%
% Figure out how many depths we're doing this for.
%
ndepths=length(depths);
%
% Get the erosion rate in gm/(cm^2*yr) (sp.epsilon is in gm/cm^2/kyr)  
%
erosionrate=sp.epsilon/1000;
%
% Note that from now on in this function, the erosion rate is in g/cm^2/yr.
%
%
% Adjust contemporary depth to depth at start time.
%
currentdepths=depths+erosionrate*age*1000;
%
% We use a time step of 100 years for the integration in time.  
%
%deltat=0.1;
%
% We integrate until the point in time where the sample was collected.
%
tfinal=sp.tfinal;
%
% Figure out the depth spacing.
%
deltadepth=deltat*erosionrate*1000;
%
% Now, the main integration loop.
%
N14=zeros(ndepths,1);
t=tfinal-age;
while (t+deltat < tfinal)
%
% Update the elevation/latitude scaling factor to the current
% time.  Note that our t is in kyr, with negative values in the
% past, while the TDSF stuff is done in terms of years before present. 
%
  interptime=t+deltat/2;
  sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'c');
%
% Compute the production rate.  We use the mid point of the range
% of depths corresponding to the current time interval.  
%
  if (min(currentdepths-deltadepth/2)<0)
    fprintf(1,'currentdepth-deltadepth/2 is %f\n',currentdepths- ...
	    deltadepth/2);
    warning('Negative depth!');
  end
  pz14=prodz14(currentdepths-deltadepth/2,pp,sf,cp);

%
% There are two terms here.  The first term is the old inventory of cosmogenic
% nuclide, decayed for a time period of deltat.  The second term represents the
% newly generated cosmogenic nuclide, including the radioactive decay of 
% some of it that was generated early in the time period.  The 
% radioactive decay factor is:
%
% f=int(P*exp(-r(deltat-t)),t=0..deltat)=P*(1-exp(-r*deltat))/r
%
% Note that for extremely small values of deltat you could get roundoff 
% errors when taking 1-exp(-r*deltat).  This hasn't been a problem for 
% our range of deltat values.
%
% The effect of using this term is that predNXX's results are essentially
% indepedendent of deltat if the production rates are constant in time and
% the erosion rate is 0.   If the erosion rate is nonzero, then deltat must
% be small enough that very little erosion occurs during a time period, or
% else predNXX's result will depend on deltat.  
% 
%
% Update N14
%
  f14=(1.0-exp(-pp.lambda14C*deltat*1000))/pp.lambda14C;
  N14=N14*exp(-pp.lambda14C*deltat*1000)+...
      pz14*f14;
%
% Update t.
%
  t=t+deltat;
%
% Update depth
%
  currentdepths=currentdepths-deltat*1000*erosionrate;
end
%
% One final fractional time step.  deltatf is the length of this
% fractional time step.
%
deltatf=tfinal-t;
deltadepth=deltatf*erosionrate*1000;
%
% Update the elevation/latitude scaling factor to the current
% time.  Note that our t is in kyr, with negative values in the
% past, while the TDSF stuff is done in terms of years before present. 
%
interptime=t+deltatf/2;
sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'c');
%
% Compute the production rates.
%
pz14=prodz14(currentdepths-deltadepth/2,pp,sf,cp);
%
% Update N14
%
f14=(1.0-exp(-pp.lambda14C*deltatf*1000))/pp.lambda14C;
N14=N14*exp(-pp.lambda14C*deltatf*1000)+...
    pz14*f14;
