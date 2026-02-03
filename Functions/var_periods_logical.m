%
% exposed_or_not = var_periods_logical(t_start_a,mids,duration,mids_type)
% exposed_or_not = var_periods_logical(t_start_a,mids,duration,mids_type,expo_initial,expo_final,model_interval,frac_buried)
%
% Creates a logical of time intervals that are either exposed or buried
% based on specified mid-points of exposure/burial periods and 
% exposure/burial durations.
%
% Output is to be used to model nuclide concentrations through time.
%
% t_start is the model start time (in ka before present).
%
% mids is an array of mid-points for exposure/burial periods (in ka before 
% present).
%
% duration is the duration of exposure/burial (in ka) for the given
% periods.
%
% mids_type is a character input to identify whether mids and duration
% inputs are for 'exposure' or 'burial' periods.
%
% expo_initial is an optional input to specify duration of an initial
% exposure period from t_start (in ka).
%
% expo_final is an optional input to specify duration of a final
% exposure period (until present) (in ka).
%
% model_interval (in years) is optional. This should be 10, 100 or 1000 
% (default is 1000 years).
%
% frac_buried is the fraction of an exposure period that a sample is buried
% (optional).
%
%
%%

function exposed_or_not = var_periods_logical(t_start,mids,duration,mids_type,expo_initial,expo_final,model_interval,frac_buried)
  
  % Check inputs
  if (nargin < 4 || nargin > 8)
      error('var_periods_logical has the wrong number of inputs!');
  end
  if (~strcmp(mids_type,'exposure') && ~strcmp(mids_type,'burial'))
      error('mids_type should be "exposure" or "burial"!');
  end
  if (nargin < 5)
      expo_initial = [];
  end
  if (nargin < 6)
      expo_final = [];
  end
  if (nargin < 7) || isempty(model_interval)
      model_interval = 1000;
  end
  if (nargin < 8) || isempty(frac_buried)
      frac_buried = 0;
  end
  if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
      error('model_interval must be 10, 100 or 1000 years!');
  end
  
  
  if model_interval == 1000
      r_fac = -3;
  elseif model_interval == 100
      r_fac = -2;
  elseif model_interval == 10
      r_fac = -1;
  end
  
  t_start_a = round(t_start*1000,r_fac);
  dur_a = round(duration*1000,r_fac);
  mids_a = mids*1000;
  
  
  % Compute logical
  time = 0:model_interval:t_start_a;
  
  if (strcmp(mids_type,'exposure'))
      % Compute exposure periods
      exp_or_not = zeros(1,length(time));
      for a = 1:length(mids)
          period = mids_a(a)-(dur_a/2):model_interval:mids_a(a)+(dur_a/2);
          period_log = ismember(time,round(period,r_fac));
          exp_or_not(period_log) = 1;
      end
  else
      % Compute burial periods
      exp_or_not = ones(1,length(time));
      for a = 1:length(mids)
          period = mids_a(a)-(dur_a/2):model_interval:mids_a(a)+(dur_a/2);
          period_log = ismember(time,round(period,r_fac));
          exp_or_not(period_log) = 0;
      end
  end
  
  if ~isempty(expo_initial) % Add on any initial exposure
      initial = t_start_a-((expo_initial-1)*1000):model_interval:t_start_a;
      initial_log = ismember(time,initial);
      exp_or_not(initial_log) = 1;
  end
  
  if ~isempty(expo_final) % Add on any final exposure
      final = 0:model_interval:0+(expo_final*1000);
      final_log = ismember(time,final);
      exp_or_not(final_log) = 1;
  end

  exposed_or_not.logical = logical(fliplr(exp_or_not));
  
    
  % Append time and interval
  exposed_or_not.time = fliplr(time);
  exposed_or_not.values = fliplr(exp_or_not);
  exposed_or_not.interval_time = zeros(1,numel(exposed_or_not.logical)) + model_interval;
  
  
  % Account for fraction of exposure period that is buried, and
  % reduce exposed_or_not to single intervals of alternate exposure/burial
  if frac_buried>0
      expbur_change = find(diff([-1 exposed_or_not.logical -1]) ~= 0); % Find points along logical at either end of exposure periods
      interval_time = diff(expbur_change) .* model_interval; % Calculate time difference between these points
      expbur_change = expbur_change(1:end-1);
      n_intervals = numel(interval_time);
      for ii = 1:n_intervals
          reduced_time(ii) = exposed_or_not.time(expbur_change(ii)); % Determine new times
          reduced_logical(ii) = (exposed_or_not.logical(expbur_change(ii)) == 1); % Determine new logical
      end
      exposed_or_not.reduced_time_zeroBurFrac = reduced_time;
      exposed_or_not.reduced_logical_zeroBurFrac = reduced_logical;
      exp_times_buried = round(interval_time * frac_buried); % Calculate times buried for each interval
      exp_times_buried(~reduced_logical) = 0; % Keep only exposure interval times
      for iii = 1:n_intervals
          if reduced_logical(iii)
              interval_time(iii) = interval_time(iii) - exp_times_buried(iii); % Remove burial time from exposure period
              if iii<n_intervals
                  interval_time(iii+1) = interval_time(iii+1) + exp_times_buried(iii);  % Add burial time to burial period
                  reduced_time(iii+1) = reduced_time(iii+1) + exp_times_buried(iii); % Adjust start of burial period
              end
          end
          reduced_values(iii) = exposed_or_not.values(expbur_change(iii)); % Determine new values
      end
      t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
      cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  else
      % Reduce exposed_or_not to single intervals of alternate exposure/burial for improved efficiency
      expbur_change = find(diff([-1 exposed_or_not.logical -1]) ~= 0); % Find points along logical at either end of exposure periods
      interval_time = diff(expbur_change) .* model_interval; % Calculate time difference between these points
      expbur_change = expbur_change(1:end-1);
      n_intervals = numel(interval_time);
      for ii = 1:n_intervals
          reduced_time(ii) = exposed_or_not.time(expbur_change(ii)); % Determine new times
          reduced_logical(ii) = (exposed_or_not.logical(expbur_change(ii)) == 1); % Determine new logical
          reduced_values(ii) = exposed_or_not.values(expbur_change(ii)); % Determine new values
      end
      t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
      cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  end
  
  
%   % Reduce exposed_or_not to single intervals of exposure/burial for improved efficiency
%   expbur_change = find(diff([-1 exposed_or_not.logical -1]) ~= 0); % Find points along logical at either end of exposure periods
%   interval_time = diff(expbur_change) .* model_interval; % Calculate time difference between these points
%   expbur_change = expbur_change(1:end-1);
%   n_intervals = numel(interval_time); 
%   for r = 1:n_intervals
%       reduced_time(r) = exposed_or_not.time(expbur_change(r)); % Determine new times
%       reduced_logical(r) = (exposed_or_not.logical(expbur_change(r)) == 1); % Determine new logical
%       reduced_values(r) = exposed_or_not.values(expbur_change(r)); % Determine new values
%   end
%   t_expo = interval_time; t_expo(~reduced_logical) = 0; % Convert to time exposed (make burial periods zero)
%   cum_t_expo = fliplr(cumsum(fliplr(t_expo))); % Calculate cumulative exposure time (from present)
  
  
  % Append reduced time, logical, interval time and cumulative exposure time
  exposed_or_not.reduced_logical = [reduced_logical 0];
  exposed_or_not.reduced_time = [reduced_time 0];
  exposed_or_not.reduced_values = [reduced_values 0];
  exposed_or_not.reduced_interval_time = [interval_time 0];
  exposed_or_not.cum_t_expo = [cum_t_expo 0];
  
  % Get time of initial exposure
  initial_idx = find(exposed_or_not.reduced_logical,1);
  exposed_or_not.time_initial = exposed_or_not.reduced_time(initial_idx);
  
end
    