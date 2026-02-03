%
% out = periods_core(sample_data,model_inputs)
% out = periods_core(sample_data,model_inputs,model_interval)
% out = periods_core(sample_data,model_inputs,model_interval,plot_fig)
%
% Computes the predicted nuclide concentration for each sample and nuclide 
% in a bedrock core, based on specified times of exposure and floating
% parameters (including period duration). Then calculates the misfit with 
% measured nuclide concentrations.
%
% sample_data is a required struct, created using get_data.m.
%
% model_inputs are the starting values used for modelling the periods (with
% dynamic/variable duration) exposure history.
%
% model_interval is an optional input to specify the desired interval for
% nuclide calculation. This can be 10, 100 or 1000 years (default is 1000).
%
% n_iterations is the number of model iterations.
%
% plot_fig is an optional input of the figure handle generated with
% plot_core_concs.m. If included, the time series, exposure/burial periods
% and corresponding predicted nuclide concentrations are plotted for each
% modelled scenario.
%
% Output are the misfits between predicted and measured nuclide
% concentrations, for each sample.
%
%
%%

function out = periods_core(data,model_inputs,model_interval,n_iterations,plot_fig)


% If model time interval (in years) is not specified, use default
if (nargin < 4) || isempty(model_interval)
    model_interval = 1000;
end
if (model_interval ~= 1000 && model_interval ~= 100 && model_interval ~= 10)
    error('model_interval must be 10, 100 or 1000 years!');
end


% Get inputs
model_fields = fieldnames(model_inputs);
n_params = sum(contains(model_fields,'rand')); % Number of model parameters

if length(data) > length(n_params)
    DOF = length(data)-length(n_params); % Calculate degrees of freedom
else
    DOF = [];
end


% Calculate for each scenario
disp(strcat("Performing ", num2str(n_iterations)," simulations..."));

for n = 1:n_iterations

    % Get input values for each iteration
    this_model_time = model_inputs.rand_start_t(n);
    this_expo_dur = model_inputs.rand_expo_duration(n);
    if isfield(model_inputs,'rand_startexpo')
        this_start_exp = model_inputs.rand_startexpo(n);
    else
        this_start_exp = [];
    end

    % Generate exposure and burial periods
    exposed_or_not = var_periods_logical(this_model_time,model_inputs.expo_mids,this_expo_dur,'exposure',this_start_exp,model_inputs.recent_exp,model_interval);
    % exposed_or_not = var_periods_logical(this_model_time,model_inputs.expo_mids,this_expo_dur,'exposure',this_start_exp,recent_exp,model_interval,this_frac_buried);
    %%% CHECK - needs to account for sample position (relative burial when transect is included

    % Build time-varying cover structure
    icecover_logical = ~exposed_or_not.reduced_logical;
    ice_depth_m = icecover_logical .* 250; % Assume 250 m ice cover during burial (more than enough to effectively block production)
    these_inputs.cover_depth = ice_depth_m .* 100; % Convert m to cm
    these_inputs.cover_density = ones(size(these_inputs.cover_depth)) .* model_inputs.density_ice; % Array of corresponding ice density
    
    % Add snow cover if applicable
    if isfield(model_inputs,'rand_snow_depth')
        icefree_logical = logical(exposed_or_not.reduced_logical); % Find ice-free times to apply snow cover
        these_inputs.cover_depth(icefree_logical) = model_inputs.rand_snow_depth(n) * 100; % Convert m to cm
        these_inputs.cover_density(icefree_logical) = model_inputs.density_snow;
    end

    % Add subaerial erosion if applicable
    if isfield(model_inputs,'rand_subaerial_ero')
        these_inputs.subaerial_ero_rate = model_inputs.rand_subaerial_ero(n);
        these_inputs.rock_density = model_inputs.core_density;
    end

    % Add subglacial erosion if applicable
    if isfield(model_inputs,'rand_subglacial_ero')
        these_inputs.subglacial_ero_rate = model_inputs.rand_subglacial_ero(n);
        these_inputs.rock_density = model_inputs.core_density;
    end

    % Calculate predicted concentrations and misfit
    these_misfits = fit_core(exposed_or_not,these_inputs,data,plot_fig);

    % Calculate reduced chi squared from mean misfits
    if ~isempty(DOF)
        out.rChi2_1014(n) = these_misfits.misfit_1014 ./ DOF;
        out.rChi2_1026(n) = these_misfits.misfit_1026 ./ DOF;
        out.rChi2_10(n) = these_misfits.misfit_10 ./ DOF;
        out.rChi2_14(n) = these_misfits.misfit_14 ./ DOF;
        out.rChi2_26(n) = these_misfits.misfit_26 ./ DOF;
        out.rChi2_1014ratio(n) = these_misfits.misfit_1014ratio ./ DOF;
        out.rChi2_1026ratio(n) = these_misfits.misfit_1026ratio ./ DOF;
    else
        out.misfit_1014(n) = these_misfits.misfit_1014;
        out.misfit_1026(n) = these_misfits.misfit_1026;
        out.misfit_10(n) = these_misfits.misfit_10;
        out.misfit_14(n) = these_misfits.misfit_14;
        out.misfit_26(n) = these_misfits.misfit_26;
        out.misfit_1014ratio(n) = these_misfits.misfit_1014ratio;
        out.misfit_1026ratio(n) = these_misfits.misfit_1026ratio;
    end
    out.consistent_N10(n) = these_misfits.consistent_N10;
    out.consistent_N14(n) = these_misfits.consistent_N14;
    out.consistent_N26(n) = these_misfits.consistent_N26;
    
    out.predN_10(n,:) = these_misfits.predN_10;
    out.predN_14(n,:) = these_misfits.predN_14;
    out.predN_26(n,:) = these_misfits.predN_26;
    
    out.time_exposed(n) = sum(exposed_or_not.logical) .* model_interval;
    out.time_buried(n) = sum(~exposed_or_not.logical) .* model_interval;

end

out.data = data;
disp('finished.');


%%%%%%%%%%%%%%%%% Function to predict core concentrations %%%%%%%%%%%%%%%%%
  
    function out = fit_core(exposed_or_not,inputs,data,plot_fig)
        
        % Get inputs
        ee = [];
        cover.density = inputs.cover_density;
        cover.z_cm = inputs.cover_depth;
        if isfield(inputs,'subaerial_ero_rate')
            ee.e_expo = inputs.subaerial_ero_rate * inputs.rock_density;
        end
        if isfield(inputs,'subglacial_ero_rate')
            ee.e_bur = inputs.subglacial_ero_rate * inputs.rock_density;
        end

        % Calculate predicted concentrations
        predN = forward_model_ztCover_101426(exposed_or_not,data,cover,ee);


        % Determine whether samples have multiple nuclide measurements
        if isfield(predN.sum,'N10')
            logical_N10 = ~isnan(predN.sum.N10);
            nuclides = 10;
        end
        if isfield(predN.sum,'N14')
            logical_N14 = ~isnan(predN.sum.N14);
            nuclides = [nuclides,14];
        end
        if isfield(predN.sum,'N26')
            logical_N26 = ~isnan(predN.sum.N26);
            nuclides = [nuclides,26];
        end
        logical_101426 = [logical_N10;logical_N14;logical_N26];
        

        % Calculate misfit
        misfits_10 = NaN(size(data)); misfits_14 = NaN(size(data)); misfits_26 = NaN(size(data)); misfits_1014ratio = NaN(size(data)); misfits_1026ratio = NaN(size(data));
        consistent_N10 = false(size(data)); consistent_N14 = false(size(data)); consistent_N26 = false(size(data));
        for i = 1:numel(data)
            this_log = logical_101426(:,i);
            if all(this_log)
                misfits_10(i) = ((predN.sum.N10(i) - data{i}.N10)./data{i}.dN10) .^2;
                misfits_14(i) = ((predN.sum.N14(i) - data{i}.N14)./data{i}.dN14) .^2;
                misfits_26(i) = ((predN.sum.N26(i) - data{i}.N26)./data{i}.dN26) .^2;

                predN_1014ratio = predN.sum.N14(i)./predN.sum.N10(i);
                meas_1014ratio = data{i}.N14./data{i}.N10;
                misfits_1014ratio(i) = predN_1014ratio - meas_1014ratio;

                predN_1026ratio = predN.sum.N26(i)./predN.sum.N10(i);
                meas_1026ratio = data{i}.N26./data{i}.N10;
                misfits_1026ratio(i) = predN_1026ratio - meas_1026ratio;

                consistent_N10(i) = (predN.sum.N10(i) >= (data{i}.N10 - data{i}.dN10)) && ...
                    (predN.sum.N10(i) <= (data{i}.N10 + data{i}.dN10));
                consistent_N14(i) = (predN.sum.N14(i) >= (data{i}.N14 - data{i}.dN14)) && ...
                    (predN.sum.N14(i) <= (data{i}.N14 + data{i}.dN14));
                consistent_N26(i) = (predN.sum.N26(i) >= (data{i}.N26 - data{i}.dN26)) && ...
                    (predN.sum.N26(i) <= (data{i}.N26 + data{i}.dN26));
            else
                if logical_N10(i)
                    misfits_10(i) = ((predN.sum.N10(i) - data{i}.N10)./data{i}.dN10) .^2;
                    consistent_N10(i) = (predN.sum.N10(i) >= (data{i}.N10 - data{i}.dN10)) && ...
                                (predN.sum.N10(i) <= (data{i}.N10 + data{i}.dN10));
                end
                if logical_N14(i)
                    misfits_14(i) = ((predN.sum.N14(i) - data{i}.N14)./data{i}.dN14) .^2;
                    consistent_N14(i) = (predN.sum.N14(i) >= (data{i}.N14 - data{i}.dN14)) && ...
                                (predN.sum.N14(i) <= (data{i}.N14 + data{i}.dN14));
                end
                if logical_N26(i)
                    misfits_26(i) = ((predN.sum.N26(i) - data{i}.N26)./data{i}.dN26) .^2;
                    consistent_N26(i) = (predN.sum.N26(i) >= (data{i}.N26 - data{i}.dN26)) && ...
                                (predN.sum.N26(i) <= (data{i}.N26 + data{i}.dN26));
                end
                if logical_N10(i) && logical_N14(i)
                    predN_1014ratio = predN.sum.N14(i)./predN.sum.N10(i);
                    meas_1014ratio = data{i}.N14./data{i}.N10;
                    misfits_1014ratio(i) = predN_1014ratio - meas_1014ratio;
                end
                if logical_N10(i) && logical_N26(i)
                    predN_1026ratio = predN.sum.N26(i)./predN.sum.N10(i);
                    meas_1026ratio = data{i}.N26./data{i}.N10;
                    misfits_1026ratio(i) = predN_1026ratio - meas_1026ratio;
                end
            end
        end
       
        out.misfit_1014 = mean([misfits_10,misfits_14],'omitnan');
        out.misfit_1026 = mean([misfits_10,misfits_26],'omitnan');
        out.misfit_10 = mean(misfits_10,'omitnan');
        out.misfit_14 = mean(misfits_14,'omitnan');
        out.misfit_26 = mean(misfits_26,'omitnan');
        out.misfit_1014ratio = mean(misfits_1014ratio,'omitnan');
        out.misfit_1026ratio = mean(misfits_1026ratio,'omitnan');
        out.consistent_N10 = any(consistent_N10) && numel(find(consistent_N10))>1;
        out.consistent_N14 = any(consistent_N14) && numel(find(consistent_N14))>1;
        out.consistent_N26 = any(consistent_N26) && numel(find(consistent_N26))>1;

        out.predN_10 = predN.sum.N10;
        out.predN_14 = predN.sum.N14;
        out.predN_26 = predN.sum.N26;


        if ~isempty(plot_fig)
            
            % Plot time-series data
            ts_plot = plot_expbur_time(exposed_or_not,plot_fig);
            
            
            % Plot predicted concentrations for samples
            plot_core_concs_pred(data,plot_fig,predN.sum,nuclides,logical_101426,[]);
            
            
            % Plot exposure/burial periods
            plot_ExpBur(exposed_or_not,ts_plot,plot_fig);
            
            
            drawnow;
            
        end
    end

end
