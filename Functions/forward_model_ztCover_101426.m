%
% predN = forward_model_ztCover_101426(exposed_or_not,data,cover)
% predN = forward_model_ztCover_101426(exposed_or_not,data,cover,ee)
%
% Calculates predicted nuclide concentrations for each sample and nuclide
% by integrating over time and depth. Includes time-varying burial.
%
% It uses generated time intervals of exposure/burial (exposed_or_not).
%
% data is a required struct, created using get_data.m.
%
% cover is a struct containing surface cover depth (cover.z_cm) and 
% density (cover.density). Default is zero cover.
%
% ee is an optional struct of erosion rates, for periods of exposure
% (ee.e_expo) and burial (ee.e_bur). Defaults are zero.
%
% Output is a struct of predicted exposure ages for all time 
% intervals and summed intervals, and then for each nuclide.
%
%
%%

function predN = forward_model_ztCover_101426(exposed_or_not,data,cover,ee)

  % Suppress warnings
  warning('off','MATLAB:integral2:funVectorization');

  % Check inputs
  if (nargin < 2 || nargin > 4)
      error('forward_model_ztCover_101426 has wrong number of inputs!');
  end
  
  if (nargin < 3) || isempty(ee)
      ee.e_expo = 0; % If ee is missing, set erosion rate during exposure to zero
      ee.e_bur = 0;  % If ee is missing, set erosion rate during burial to zero
  end

  % Do calculations at model interval resolution if variable ice cover is used
  % (i.e. if more than three cover depth values (zero, snow value, ice value))
  if nargin > 3 && isscalar(cover.z_cm) && unique(cover.z_cm) > 3
      interval_time = exposed_or_not.interval_time;
      exp_logical = exposed_or_not.logical;
      cum_t_expo = interval_time;
      timeBP = exposed_or_not.time;
  else % Otherwise, use reduced intervals for faster computation
      interval_time = exposed_or_not.reduced_interval_time;
      exp_logical = exposed_or_not.reduced_logical;
      cum_t_expo = exposed_or_not.cum_t_expo;
      timeBP = exposed_or_not.reduced_time;
  end

  % Calculate mass depth of cover above rock sample surface
  if nargin < 3 || isempty(cover)
      mass_depth_cover = zeros(size(interval_time)); % Assume no cover
  elseif isscalar(cover.z_cm)
      mass_depth_cover = cover.density * cover.z_cm * ones(size(interval_time)); % Constant cover depth during burial
  else
      mass_depth_cover = cover.density .* cover.z_cm; % Time-varying cover during burial
  end
  mass_depth_cover = max(0, mass_depth_cover); % Ensure no negatives


  n_intervals = numel(exp_logical);
  n_samples = nnz(~cellfun(@isempty, data));

  N_10 = zeros(n_intervals,n_samples); % Row for each time interval, column for each sample
  N_14 = N_10; N_26 = N_10;  

  % Calculate for each interval
  for a = 1:n_intervals
      
      % Exposure interval
      if exp_logical(a) == 1
          
          for b = 1:n_samples

              % Compute depths at the end of this exposure period
              if a < n_intervals % Account for any subaerial erosion and snow cover, unless final interval
                  this_top_z = data{b}.top_z_gcm2 + ee.e_expo.*cum_t_expo(a+1) + mass_depth_cover(a);
                  this_bottom_z = data{b}.bottom_z_gcm2 + ee.e_expo.*cum_t_expo(a+1) + mass_depth_cover(a);
              else
                  this_top_z = data{b}.top_z_gcm2 + mass_depth_cover(a);
                  this_bottom_z = data{b}.bottom_z_gcm2 + mass_depth_cover(a);
              end             
                         
              if data{b}.nuclide10 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda10Be;
                  pr_func = @(z,t) PR_Z((z + ee.e_expo.*t),data{b}.pp,data{b}.sf10{1},data{b}.cp10{1},10) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight10;
                  this_N_10 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_10(a,b) = this_N_10 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_10(a,b) = this_N_10;
                  end

              end
              
              if data{b}.nuclide14 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda14C;
                  pr_func = @(z,t) PR_Z((z + ee.e_expo.*t),data{b}.pp,data{b}.sf14{1},data{b}.cp14{1},14) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight14;
                  this_N_14 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_14(a,b) = this_N_14 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_14(a,b) = this_N_14;
                  end

              end

              if data{b}.nuclide26 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda26Al;
                  pr_func = @(z,t) PR_Z((z + ee.e_expo.*t),data{b}.pp,data{b}.sf26{1},data{b}.cp26{1},26) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight26;
                  this_N_26 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_26(a,b) = this_N_26 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_26(a,b) = this_N_26;
                  end

              end
              
              clear Integrated_N;              
          end
      
      
      % Burial interval
      else
          
          for b = 1:n_samples

              % Compute depths at the end of this burial period
              if a < n_intervals % Account for ice cover and any subglacial erosion, unless final interval
                  this_top_z = data{b}.top_z_gcm2 + ee.e_bur.*cum_t_expo(a+1) + mass_depth_cover(a);
                  this_bottom_z = data{b}.bottom_z_gcm2 + ee.e_bur.*cum_t_expo(a+1) + mass_depth_cover(a);
              else
                  this_top_z = data{b}.top_z_gcm2 + mass_depth_cover(a);
                  this_bottom_z = data{b}.bottom_z_gcm2 + mass_depth_cover(a);
              end             
                         
              if data{b}.nuclide10 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda10Be;
                  pr_func = @(z,t) PR_Z((z + ee.e_bur.*t),data{b}.pp,data{b}.sf10{1},data{b}.cp10{1},10) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight10;
                  this_N_10 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_10(a,b) = this_N_10 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_10(a,b) = this_N_10;
                  end
              end
              
              if data{b}.nuclide14 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda14C;
                  pr_func = @(z,t) PR_Z((z + ee.e_bur.*t),data{b}.pp,data{b}.sf14{1},data{b}.cp14{1},14) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight14;
                  this_N_14 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_14(a,b) = this_N_14 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_14(a,b) = this_N_14;
                  end

              end

              if data{b}.nuclide26 == 1

                  % Define production rate functions
                  l = data{b}.pp.lambda26Al;
                  pr_func = @(z,t) PR_Z((z + ee.e_bur.*t),data{b}.pp,data{b}.sf26{1},data{b}.cp26{1},26) .* exp(-l.*t);

                  % Integrate over time and depth for each sample
                  Integrated_N = zeros(size(this_top_z));
                  for c = 1:length(this_top_z)
                      % Integrate each sample from top depth to bottom depth and from zero to the exposure time in this time step.
                      % Divide integral by sample thickness.
                      Integrated_N(c) = integral2(pr_func,this_top_z(c),this_bottom_z(c),0,interval_time(a),...
                          'Method','iterated', 'RelTol',1e-6,'AbsTol',1e-6) ./ (this_bottom_z(c) - this_top_z(c));
                  end

                  % Average by mineral weight (only necessary if multiple weights for a sample)
                  this_weight = data{b}.weight26;
                  this_N_26 = ( (sum(Integrated_N.*this_weight))./sum(this_weight) );

                  % Account for decay from end of interval to present 
                  % (only if there are more than one total interval + 0)
                  if (n_intervals > 2) && (a < n_intervals)
                      timeBP_intervalEnd = timeBP(a) - interval_time(a);
                      N_26(a,b) = this_N_26 .* exp(-l .* timeBP_intervalEnd);
                  else
                      N_26(a,b) = this_N_26;
                  end

              end
              
              clear Integrated_N;              
          end
       
      end
  end

  
  predN.all.N10 = N_10;
  predN.all.N14 = N_14;
  predN.all.N26 = N_26;
    
  % Sum all intervals
  predN.sum.N10 = sum(N_10);
  predN.sum.N14 = sum(N_14);
  predN.sum.N26 = sum(N_26);

  % Remove zeros (nuclides not measured)
  predN.sum.N10(predN.sum.N10==0) = NaN;
  predN.sum.N14(predN.sum.N14==0) = NaN;
  predN.sum.N26(predN.sum.N26==0) = NaN;

  % Add time and cover depth
  predN.interval_times = interval_time;
  predN.mass_depth_cover = mass_depth_cover;
  
end
