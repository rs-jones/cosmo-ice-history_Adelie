%
% out = run_periods(sample_type,sample_data_core,sample_data_transect,model_inputs,model_interval,n_iterations,plotfig)
% out = run_periods(sample_type,sample_data_core,sample_data_transect,model_inputs,model_interval,n_iterations,plotfig,erratics_height)
%
% Computes the predicted nuclide concentration for each sample and nuclide,
% based on specified times of exposure and floating parameters (including
% period duration). Then calculates the misfit with measured nuclide
% concentrations.
%
% sample_type is a required input specifying 'core', 'transect' or
% ['core','transect'].
%
% sample_data_core and sample_data_transect are required structs for the
% bedrock core and elevation transect samples, created using get_data.m.
% Use empty field input (i.e. []) if not including core/transect data.
%
% model_inputs are the starting values used for modelling the periods
% exposure history.
%
% model_interval is the input to specify the desired interval for
% nuclide calculation. This can be 10, 100 or 1000 years (default is 1000).
%
% n_iterations is the number of model iterations.
%
% plotfig is a true/false logical to specify if a figure should be plotted
% (updated each model iteration).
%
% erratics_height is an optional input to specify the surface height of any
% erratics in the transect relative to the ground (bedrock core), in
% metres. Default is zero.
%
% Output are the misfits between predicted and measured nuclide
% concentrations, for each sample.
%
%
%%

function out = run_periods(sample_type,sample_data_core,sample_data_transect,model_inputs,model_interval,n_iterations,plotfig,erratics_height)

if any(ismember(sample_type,"core")) && isempty(sample_data_core)
    error('"core" specified but sample data not provided.')
end
if any(ismember(sample_type,"transect")) && isempty(sample_data_transect)
    error('"transect" specified but sample data not provided.')
end

if exist("erratics_height",'var') && isempty(erratics_height)
    erratics_height = 0;
end

% Check that nuclides exist to model
nuclides = model_inputs.nuclides; % Nuclides to model

has10 = (any(ismember(sample_type,"transect")) && any(sample_data_transect.logical_10)) || ...
    (any(ismember(sample_type,"core")) && any(sample_data_core.logical_10));
has14 = (any(ismember(sample_type,"transect")) && any(sample_data_transect.logical_14)) || ...
    (any(ismember(sample_type,"core")) && any(sample_data_core.logical_14));
has26 = (any(ismember(sample_type,"transect")) && any(sample_data_transect.logical_26)) || ...
    (any(ismember(sample_type,"core")) && any(sample_data_core.logical_26));

for i = 1:length(nuclides)
    this_nuclide = nuclides(i);

    if any(ismember(["10","10Be","Be10","Be-10"],this_nuclide)) && ~has10
        error('Be10 data do not exist to model.')
    end
    if any(ismember(["14","14C","C14","C-14"],this_nuclide)) && ~has14
        error('C14 data do not exist to model.')
    end
    if any(ismember(["26","26Al","Al26","Al-26"],this_nuclide)) && ~has26
        error('Al26 data do not exist to model.')
    end
end

% Re-organise core data for each sample
if any(ismember(sample_type,"core"))
    core_data = cell(1,length(sample_data_core.s));
    sample_data_core.logical_1014 = any(sample_data_core.logical_10 & sample_data_core.logical_14);
    for a = 1:length(sample_data_core.s)
        core_data{a}.name = sample_data_core.s{a}.name{1};
        core_data{a}.pp = sample_data_core.pp;
        core_data{a}.sf10 = sample_data_core.sf1026(a);
        core_data{a}.cp10 = sample_data_core.cp1026(a);
        core_data{a}.sf26 = sample_data_core.sf1026(a);
        core_data{a}.cp26 = sample_data_core.cp1026(a);
        core_data{a}.sf14 = sample_data_core.sf14(a);
        core_data{a}.cp14 = sample_data_core.cp14(a);
        core_data{a}.top_z_cm = sample_data_core.s{a}.top_z;
        core_data{a}.bottom_z_cm = sample_data_core.s{a}.bottom_z;
        core_data{a}.top_z_gcm2 = sample_data_core.s{a}.top_z_gcm2;
        core_data{a}.bottom_z_gcm2 = sample_data_core.s{a}.bottom_z_gcm2;
        if isfield(sample_data_core.s{a},'weight10')
            core_data{a}.weight10 = sample_data_core.s{a}.weight10;
        end
        if isfield(sample_data_core.s{a},'weight14')
            core_data{a}.weight14 = sample_data_core.s{a}.weight14;
        end
        if isfield(sample_data_core.s{a},'weight26')
            core_data{a}.weight26 = sample_data_core.s{a}.weight26;
        end
        if any(ismember(["10","10Be","Be10","Be-10"],nuclides))
            core_data{a}.nuclide10 = sample_data_core.s{a}.nuclide10;
        else
            core_data{a}.nuclide10 = 0;
        end
        if any(ismember(["14","14C","C14","C-14"],nuclides))
            core_data{a}.nuclide14 = sample_data_core.s{a}.nuclide14;
        else
            core_data{a}.nuclide14 = 0;
        end
        if any(ismember(["26","26Al","Al26","Al-26"],nuclides))
            core_data{a}.nuclide26 = sample_data_core.s{a}.nuclide26;
        else
            core_data{a}.nuclide26 = 0;
        end
        if core_data{a}.nuclide10 == 1
            core_data{a}.N10 = sample_data_core.s{a}.N10;
            core_data{a}.dN10 = sample_data_core.s{a}.dN10;
        end
        if  core_data{a}.nuclide14 == 1
            core_data{a}.N14 = sample_data_core.s{a}.N14;
            core_data{a}.dN14 = sample_data_core.s{a}.dN14;
        end
        if  core_data{a}.nuclide26 == 1
            core_data{a}.N26 = sample_data_core.s{a}.N26;
            core_data{a}.dN26 = sample_data_core.s{a}.dN26;
        end
    end
    out.core_data = core_data;
    if isfield(sample_data_core.CC,"Be10")
        out.core_density = mean(sample_data_core.CC.Be10(:,6),'omitmissing');
    elseif isfield(sample_data_core.CC,"C14")
        out.core_density = mean(sample_data_core.CC.C14(:,6),'omitmissing');
    else
        out.core_density = mean(sample_data_core.CC.Al26(:,6),'omitmissing');
    end
end

% Re-organise transect data for each sample
if any(ismember(sample_type,"transect"))
    transect_data = cell(1,length(sample_data_transect.s));
    sample_data_transect.logical_1014 = any(sample_data_transect.logical_10 & sample_data_transect.logical_14);
    for a = 1:length(sample_data_transect.s)
        transect_data{a}.name = sample_data_transect.s{a}.name{1};
        transect_data{a}.pp = sample_data_transect.pp;
        transect_data{a}.sf10 = sample_data_transect.sf1026(a);
        transect_data{a}.cp10 = sample_data_transect.cp1026(a);
        transect_data{a}.sf26 = sample_data_transect.sf1026(a);
        transect_data{a}.cp26 = sample_data_transect.cp1026(a);
        transect_data{a}.sf14 = sample_data_transect.sf14(a);
        transect_data{a}.cp14 = sample_data_transect.cp14(a);
        transect_data{a}.top_z_cm = sample_data_transect.s{a}.top_z;
        transect_data{a}.bottom_z_cm = sample_data_transect.s{a}.bottom_z;
        transect_data{a}.top_z_gcm2 = sample_data_transect.s{a}.top_z_gcm2;
        transect_data{a}.bottom_z_gcm2 = sample_data_transect.s{a}.bottom_z_gcm2;
        if isfield(sample_data_transect.s{a},'weight10')
            transect_data{a}.weight10 = sample_data_transect.s{a}.weight10;
        end
        if isfield(sample_data_transect.s{a},'weight14')
            transect_data{a}.weight14 = sample_data_transect.s{a}.weight14;
        end
        if isfield(sample_data_transect.s{a},'weight26')
            transect_data{a}.weight26 = sample_data_transect.s{a}.weight26;
        end
        if any(ismember(["10","10Be","Be10","Be-10"],nuclides))
            transect_data{a}.nuclide10 = sample_data_transect.s{a}.nuclide10;
        else
            transect_data{a}.nuclide10 = 0;
        end
        if any(ismember(["14","14C","C14","C-14"],nuclides))
            transect_data{a}.nuclide14 = sample_data_transect.s{a}.nuclide14;
        else
            transect_data{a}.nuclide14 = 0;
        end
        if any(ismember(["26","26Al","Al26","Al-26"],nuclides))
            transect_data{a}.nuclide26 = sample_data_transect.s{a}.nuclide26;
        else
            transect_data{a}.nuclide26 = 0;
        end
        if transect_data{a}.nuclide10 == 1
            transect_data{a}.N10 = sample_data_transect.s{a}.N10;
            transect_data{a}.dN10 = sample_data_transect.s{a}.dN10;
        end
        if transect_data{a}.nuclide14 == 1
            transect_data{a}.N14 = sample_data_transect.s{a}.N14;
            transect_data{a}.dN14 = sample_data_transect.s{a}.dN14;
        end
        if transect_data{a}.nuclide26 == 1
            transect_data{a}.N26 = sample_data_transect.s{a}.N26;
            transect_data{a}.dN26 = sample_data_transect.s{a}.dN26;
        end
    end
    out.transect_data = transect_data;
    if isfield(sample_data_transect.CC,"Be10")
        out.transect_density = mean(sample_data_transect.CC.Be10(:,6),'omitmissing');
    elseif isfield(sample_data_transect.CC,"C14")
        out.transect_density = mean(sample_data_transect.CC.C14(:,6),'omitmissing');
    else
        out.transect_density = mean(sample_data_transect.CC.Al26(:,6),'omitmissing');
    end
end

% Add age of recent exposure if present
if ~isempty(model_inputs.recent_expo_sample)

    if (any(ismember(sample_type,"core")) && ~isfield(sample_data_core,'ages')) || (any(ismember(sample_type,"transect")) && ~isfield(sample_data_transect,'ages'))
        error('Exposure ages are not provided in the input data.')
    end

    % Find sample and get ages
    if any(ismember(sample_type,"core")) && ~any(ismember(sample_type,"transect"))
        recent_samp_idx = find(sample_data_core.names==model_inputs.recent_expo_sample);
        sample_ages = sample_data_core.ages;
    elseif ~any(ismember(sample_type,"core")) && any(ismember(sample_type,"transect"))
        recent_samp_idx = find(sample_data_transect.names==model_inputs.recent_expo_sample);
        sample_ages = sample_data_transect.ages;
    else
        recent_samp_idx(1) = find(sample_data_core.names==model_inputs.recent_expo_sample);
        sample_ages{1} = sample_data_core.ages;
        recent_samp_idx(2) = find(sample_data_transect.names==model_inputs.recent_expo_sample);
        sample_ages{2} = sample_data_transect.ages;
    end

    if has14 % Take 14C age as default
        if numel(recent_samp_idx) == 1
            out.recent_exp = sample_ages.C14(recent_samp_idx,1)/1000;
        else
            try out.recent_exp = sample_ages{1}.C14(recent_samp_idx(1),1)/1000;
            catch
                out.recent_exp = sample_ages{2}.C14(recent_samp_idx(2),1)/1000;
            end
        end
    elseif has10 % Otherwise, use 10Be age
        if numel(recent_samp_idx) == 1
            out.recent_exp = sample_ages.Be10(recent_samp_idx,1)/1000;
        else
            try out.recent_exp = sample_ages{1}.Be10(recent_samp_idx(1),1)/1000;
            catch
                out.recent_exp = sample_ages{2}.Be10(recent_samp_idx(2),1)/1000;
            end
        end
    else % Otherwise, use 26Al age
        if numel(recent_samp_idx) == 1
            out.recent_exp = sample_ages.Al26(recent_samp_idx,1)/1000;
        else
            try out.recent_exp = sample_ages{1}.Al26(recent_samp_idx(1),1)/1000;
            catch
                out.recent_exp = sample_ages{2}.Al26(recent_samp_idx(2),1)/1000;
            end
        end
    end
else
    out.recent_exp = [];
end


% Create random model inputs
out.rand_start_t = param_rand_gen(model_inputs.start_time_bnds(1),model_inputs.start_time_bnds(2),n_iterations);
out.rand_expo_duration = param_rand_gen(model_inputs.expo_dur_bnds(1),model_inputs.expo_dur_bnds(2),n_iterations);
if ~isempty(model_inputs.startexpo_bnds)
    out.rand_startexpo = param_rand_gen(model_inputs.startexpo_bnds(1),model_inputs.startexpo_bnds(2),n_iterations);
end
if ~isempty(model_inputs.snow_depth_bnds)
    out.rand_snow_depth = param_rand_gen(model_inputs.snow_depth_bnds(1),model_inputs.snow_depth_bnds(2),n_iterations);
end
if ~isempty(model_inputs.subaerial_erosion_bnds)
    out.rand_subaerial_ero = param_rand_gen(model_inputs.subaerial_erosion_bnds(1),model_inputs.subaerial_erosion_bnds(2),n_iterations);
end
if ~isempty(model_inputs.subglacial_erosion_bnds)
    out.rand_subglacial_ero = param_rand_gen(model_inputs.subglacial_erosion_bnds(1),model_inputs.subglacial_erosion_bnds(2),n_iterations);
end

out.density_snow = 0.27; % g cm-3
out.density_ice = 0.917; % g cm-3
out.expo_mids = model_inputs.expo_mids;

% Determine burial fractions for each sample depending on relative elevation difference (not relevant if analysing a single core)
if isfield(model_inputs,'burialfrac_bnds') && ~isempty(model_inputs.burialfrac_bnds)
    out.rand_burialfrac = param_rand_gen(model_inputs.burialfrac_bnds(1),model_inputs.burialfrac_bnds(2),n_iterations);

    sample_positions = sample_data_transect.Position;
    if any(ismember(sample_type,"core"))
        sample_positions = [sample_positions,sample_data_core.Position];
    end
    pos_scaled = rescale(sample_positions);
    for i = 1:length(out.rand_burialfrac)
        out.sample_burial_fracs(i,:) = out.rand_burialfrac(i)*pos_scaled;
    end
end


% Run model for core samples
if any(ismember(sample_type,"core"))
    disp(' '); disp('Periods history (variable duration) for bedrock core');

    if plotfig
        core_fig = plot_core_concs(sample_data_core,1);
        drawnow
    end

    core_inputs = out;
    out.bestfits_core = periods_core(core_data,core_inputs,model_interval,n_iterations,core_fig);

    out.bestfits_core.consistent_1014 = out.bestfits_core.consistent_N10 & out.bestfits_core.consistent_N14;
    out.bestfits_core.consistent_1026 = out.bestfits_core.consistent_N10 & out.bestfits_core.consistent_N26;
end


% % Run model for transect samples
% if any(ismember(sample_type,"transect"))
%     disp(' '); disp('Periods history (variable duration) for elevation transect');
%
%     transect_inputs = out;
%     transect_inputs.sample_elevations = transect_elevations;
%     sample_height_above_bedrock = erratics_height;
%     if ~isempty(model_inputs.snow_depth_bnds)
%         transect_inputs.rand_snow_depth = out.rand_snow_depth - sample_height_above_bedrock; % Reduce snow cover for erratics relative to bedrock
%         transect_inputs.rand_snow_depth(transect_inputs.rand_snow_depth < 0) = 0;
%     end
%     out.bestfits_transect = transect_multistage_MC(transect_data,transect_inputs,model_interval,n_iterations,[]);
%
%     out.bestfits_transect.consistent_1014 = out.bestfits_transect.consistent_N10 & out.bestfits_transect.consistent_N14;
% end
%
%
% % Evaluate consistency between core and transect results
% if any(ismember(sample_type,"core")) && any(ismember(sample_type,"transect"))
%     out.bestfit_coreANDtransect.consistent_N10 = out.bestfits_core.consistent_N10 & out.bestfits_transect.consistent_N10;
%     out.bestfit_coreANDtransect.consistent_N14 = out.bestfits_core.consistent_N14 & out.bestfits_transect.consistent_N14;
%     out.bestfit_coreANDtransect.consistent_1014 = out.bestfits_core.consistent_1014 & out.bestfits_transect.consistent_1014;
%     out.bestfit_coreANDtransect.consistent_core14_transect1014 = out.bestfits_core.consistent_N14 & out.bestfits_transect.consistent_1014;
% end


% Export nuclides modelled
out.nuclides_modelled = nuclides;


end