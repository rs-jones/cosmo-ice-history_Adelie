%
% cp=comppars14_mod(physpars,samppars,scalefactors,maxdepth)
%
% Generates a complete set of parameters for carbon-14 dating
% from the basic physics parameters, sample specific parameters, 
% and scale factors.
%
% maxdepth is an optional parameter.  The default value is 2500
% g/cm^2, or about 10 meters.  Computing muon production profiles
% down this far is extremely expensive.  For surface samples, a max
% depth of 25 g/cm^2 (or about 10cm)  is more appropriate.
%
% MODIFICATION to use Balco's model for muon production.
%
function cp=comppars14_mod(pp,sp14,sf14,maxdepth)
%
% First, set maxdepth to a default value if not supplied.
%
if (nargin < 4)
  maxdepth=2500;
end
%
% Record the maximum depth so that we know how deep it is safe to go.
%
cp.maxdepth=maxdepth;
%
% Setup Lambdafe.  
%
cp.Lambdafe=sp14.Lambdafe;
%

cp.Zs=sp14.ls*sp14.rb;


% % The RcEst above was originally used and has been replaced by Nat Lifton's
% % new model.  He  fit the function below to trajectory-traced cutoffs for the
% % modern dipolar field strength, so it should be more accurate overall.
% 
% dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
% %   Trajectory-traced dipolar estimate for these purposes
%    RcEst = (dd(1)*cos(d2r(sp14.latitude)) + ...
%        dd(2)*(cos(d2r(sp14.latitude))).^2 + ...
%        dd(3)*(cos(d2r(sp14.latitude))).^3 + ...
%        dd(4)*(cos(d2r(sp14.latitude))).^4 + ...
%        dd(5)*(cos(d2r(sp14.latitude))).^5 + ...
%        dd(6)*(cos(d2r(sp14.latitude))).^6);
   
%Use mean Solar Modulation Parameter (SPhiInf)
nnn=ceil(maxdepth/5);
%
deltadepth=maxdepth/nnn;    % this gives the step for each

cp.depthvector=[0:deltadepth:maxdepth];

% Store the muons production rates from the code
%
% Preallocate space for cp.muon14.  
%
cp.muon14=zeros(3,length(cp.depthvector));
%
% Store the depth vector.
%
cp.muon14(1,:)=cp.depthvector;  
%

%
% In the following, we'll
% Generate Carbon-14 production rates.
%
if ~isnan(sp14.concentration14)
    
    % Set the constants to C-14
    consts_14.Natoms = pp.Natoms14;
    consts_14.sigma190 = pp.sigma190_14;
    % consts_14.delsigma190 = pp.delsigma190_14; % not used
    consts_14.sigma0 = pp.sigma014;
    consts_14.k_neg = pp.k_neg14;
    % consts_14.fstar = pp.fstar14;

    % muon14=muonfluxsato(cp.depthvector,sf14.P,RcEst,pp.SPhiInf,pp,'yes');
    % 
    % % Now calculate the production rates.
    % z=cp.depthvector;
    % 
    % % Store muon scaling factor
    % cp.SFmufast=muon14.SFmufast;
    % cp.SFmuslow=muon14.SFmuslow;
    % 
    % % Depth-dependent parts of the fast muon reaction cross-section
    % % updated to have various possibilities for alpha
    % % Balco original - from Heisinger fast muon paper Sea Level
    % % Beta = 0.846 - 0.015 .* log((z./100)+1) + 0.003139 .* (log((z./100)+1).^2);
    % 
    % % % For Beacon heights
    % % aalpha = 0.75;
    % % Beta =  0.842344 - 0.0107398 log((z./100)+1) + 0.00240182 log((z./100)+1)^2
    % %
    % % aalpha = 0.85;
    % % Beta =  0.888695 - 0.00716992 log((z./100)+1) + 0.00169676 log((z./100)+1)^2
    % %
    % aalpha = 1.0;
    % Beta = 1.0;
    % %
    % % alpha = 1.15;
    % % Beta =  1.18129 + 0.00903804 log((z./100)+1) - 0.00273586 lnlog((z./100)+1)^2
    % %
    % % aalpha = 1.30;
    % % Beta =  1.45283+0.04615 lnlog((z./100)+1) - 0.0153481 lnlog((z./100)+1)^2 + 0.000672339 lnlog((z./100)+1)^3
    % % %
    % 
    % Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));
    % 
    % %sigma0 = consts.sigma190./(190.^aalpha);
    % 
    % % fast muon production
    % 
    % P_fast = muon14.phi.*Beta.*(Ebar.^aalpha).*sigma0.*Natoms;
    % 
    % % negative muon capture
    % P_neg = muon14.R.*k_negpartial.*fstar;
    % 
    % cp.P_total14=P_fast+P_neg;

    z = cp.depthvector;
    for i = 1:length(z)
        cp.P_total14(i) = P_mu_total_alpha1(z(i), sf14.P, consts_14); % Balco's muon model
    end

    cp.muon14(2,:) = cp.P_total14;
    % cp.negflux=muon14.R;
    % cp.totalflux=muon14.phi;

end

