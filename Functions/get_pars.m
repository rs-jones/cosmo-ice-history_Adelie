%
% sample_data = get_pars(sample_data,scaling_model,max_depth)  
%
% Gets parameters for each sample using CRONUScalc functions:
% pp=physical constants
% sp=sample parameters
% sf=scaling factors
% cp=computed parameters
% Exports to sample_data struct.
%
% Requires sample data to be loaded and sorted with get_data.m, and a
% scaling model to be set (i.e. 'DE','DU','LI','ST','LM','LSD'/'SF',
% 'LSDn'/'SA'). Optionally specify the max depth to compute production (in
% g cm^2).
%
% Written by Richard Selwyn Jones, Monash University
% richard.s.jones@monash.edu
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function sample_data = get_pars(sample_data,scaling_model,max_depth)

  warning('off','MATLAB:quad:ImproperFcnValue');

  disp('Generating nuclide parameters...');

  if strcmpi(scaling_model,'LSD')
      scaling = 'SF';
  elseif strcmpi(scaling_model,'LSDn')
      scaling = 'SA';
  else
      scaling = scaling_model;
  end
  if nargin < 3 || isempty(max_depth)
      maxdepth_gcm2 = 2500; % Specifies the maximum depth (in g/cm^2) for which production rates will be computed
  else
      maxdepth_gcm2 = max_depth;
  end

  
  % Find number of samples
  n_samples = numel(sample_data.CC.Be10Al26(:,1));

  % Get the physical parameters
  sample_data.pp = physpars(scaling);
  
  % Update Cronus 14C (and 10Be) data
  sample_data.pp.lambda10Be = -log(0.5)./1.387e6; % Chmeleff/Korschinek value (from Balco) = 4.9975e-07
  sample_data.pp.lambda14C = -log(0.5)./5730; % Based on Balco = 1.2097e-04
  sample_data.pp.lambda26Al = 9.83e-07; % Compatible with Nishiizumi standards (from Balco) - t(1/2) = 0.705e6 yr
  sample_data.pp.sigma010 = 2.8000e-31; % Based on Balco
  sample_data.pp.sigma014 = 2.3684e-30; % Based on Balco
  sample_data.pp.sigma026 = 3.8900e-30; % Based on Balco
  sample_data.pp.PsBe = 4.086; % From recalibrating with CRONUS primary (20161204), after Balco - ST ref at SLHL
  sample_data.pp.PsC = 12.1; % Based on Tulane CRONUS A after Balco (Hippe, 2017 = 12.24) - ST ref at SLHL
  sample_data.pp.PsAl = 28.535; % From recalibrating with CRONUS primary (20161204), after Balco - ST ref at SLHL
  sample_data.pp.sigmaPsBe = 4.086.*0.079;%0.322794; % Bootstrap from secondary cal data set, after Balco - ST ref at SLHL
  sample_data.pp.sigmaPsC = 12.1 * 0.05;%0.605; % Based on Balco (5% approx.) - ST ref at SLHL
  sample_data.pp.sigmaPsAl = 28.535.*0.104; % Based on Balco - ST ref at SLHL
  sample_data.pp.k_neg10 = 0.00191 .* 0.704 .* 0.1828; % From BCO fit, after Balco
  sample_data.pp.k_neg14 = 0.116 .* 0.704 .* 0.1828; % From Leymon High fit, after Balco
  sample_data.pp.k_neg26 = 0.0133 .* 0.296 .* 0.6559; % From BCO fit, after Balco
  
  % Get sample-specific parameters
  for s = 1:n_samples
  
      % Extract the sample parameters from a sample_data vector
      sample_data.sp1026{s} = samppars1026(sample_data.CC.Be10Al26(s,:));
      sample_data.sp14{s} = samppars14(sample_data.CC.C14(s,:));
      
      % Get the scale factors
      sample_data.sf1026{s} = scalefacs1026(sample_data.sp1026{s});
      sample_data.sf14{s} = scalefacs14(sample_data.sp14{s});
      
      % Computed parameters
      sample_data.cp1026{s} = comppars1026_mod(sample_data.pp,sample_data.sp1026{s},sample_data.sf1026{s},maxdepth_gcm2); % MODIFIED to use Balco's muon model
      sample_data.cp14{s} = comppars14_mod(sample_data.pp,sample_data.sp14{s},sample_data.sf14{s},maxdepth_gcm2); % MODIFIED to use Balco's muon model
      
      % Go ahead and produce contemporary scaling factors
      sample_data.sf1026{s}.currentsf = getcurrentsf(sample_data.sf1026{s},0,scaling,'albe');
      sample_data.sf14{s}.currentsf = getcurrentsf(sample_data.sf14{s},0,scaling,'c');
      
  end
  
  % Append scaling model
  sample_data.scaling_model = scaling;
  
  % Check given scaling models are consistent
  if ~strcmpi(scaling,sample_data.scaling_model)
      warning('specified scaling model is not the same as that given in the input data.')
  end

  disp('done.');
      
end
