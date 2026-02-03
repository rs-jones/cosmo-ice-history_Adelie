%
% [ProdtotalC,ProdsC,ProdmuC]=prodz14(z,pp,sf,cp)
%
%
% The prodz function calculates the individual contribution 
% to the production of 14C from spallation and muon pathways 
% at a particular given depth. 

%  pp,sf,cp               physics, site, sample, and
%                         computed parameters.
%  z                      Depths in g/cm^2.
% 
% Get current scaling factors.
% Then call this function to compute total production for a vector of
% depths.  
% Note: This function has been updated to work with nuclide-dependent
% scaling factors (only applicable to the Sato/Lifton scaling scheme).  
%
function [ProdtotalC,ProdsC,ProdmuC]=prodz14(z,pp,sf,cp)
%
% Make sure that we don't have any depths that are too big.
%
if (max(z) > cp.maxdepth)
  error('Prodz called for depth greater than maxdepth.');
end
%
% Make sure that sf.Sel14 is a number.
%
if (isnan(sf.currentsf.Sel14))
    error('sf.currentsf.Sel14 is a NaN');
end
%
% We'll get some slightly negative depths due to roundoff errors.
% Fix them.
%
for i=1:length(z)
  if (z(i)<-1.0e-4)
    error('z(i) is negative!');
  end
  if (z(i)<0.0)
    z(i)=0.0;
  end
end

% %
% %find the number of depths given in the input vector z
% numberdepths=length(z);
% 
% 
% %For each z, find the appropriate negative muon flux and total flux terms by
% %interpolation
% %assume depths are always given in an already-sorted vector
% 
% negfluxdepth=interpolate(cp.depthvector,cp.negflux,z);
% totalfluxdepth=interpolate(cp.depthvector,cp.totalflux,z);


% Ps(z)=Sel*ST*Ps_0*exp(-z/BigLambda_f)
%
% Introduce a new term to speed up computation
%
expfactor=exp(-z/cp.Lambdafe);
%
% New productions for each spallation pathway separately, then deal with
% muons

ProdsC=sf.currentsf.Sel14*sf.ST*pp.PsC*expfactor;
%
% Muons are now being multiplied by the terrain shielding factor.  
%
ProdmuC=sf.ST*interpolate(cp.muon14(1,:),cp.muon14(2,:),z);

%
% Now, compute the total production.
%
ProdtotalC=ProdsC+ProdmuC;
%
% Look for NaN's in the results.
%
if ((sum(isnan(ProdtotalC)) > 0) )
  ProdsC
  ProdmuC
  error('Prodz14 produced NaN!');
end


