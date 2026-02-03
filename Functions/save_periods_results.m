%
% save_periods_results(base_save_name,model_inputs,n_iterations,MC_results)
%
% Builds a descriptive filename from model input bounds and saves the
% Monte Carlo output structure. This keeps the main script clean while
% ensuring saved files encode the model configuration used.
%
% base_save_name is a user-defined base name for the output file.
%
% model_inputs is the inputs structure used to run run_periods.
%
% n_iterations is the number of Monte Carlo iterations.
%
% MC_results is the output structure returned by run_periods.
%
% Optional parameter bounds are appended only if present.
%
% Output file format:
%   <base>_iterX_startBndsA-B_expoDurBndsC-D_[optional fields].mat
%
% Saved variable name inside the file remains MC_results.
%
%%

function save_periods_results(base_save_name,model_inputs,n_iterations,MC_results)

% Required bounds
start_str = sprintf('%g-%g', ...
    model_inputs.start_time_bnds(1), ...
    model_inputs.start_time_bnds(2));

expo_str = sprintf('%g-%g', ...
    model_inputs.expo_dur_bnds(1), ...
    model_inputs.expo_dur_bnds(2));

% Start assembling filename
save_str = strcat( ...
    base_save_name, ...
    '_iter',num2str(n_iterations), ...
    '_startBnds',start_str, ...
    '_expoDurBnds',expo_str);

% Optional initial exposure duration
if isfield(model_inputs,'startexpo_bnds') && ...
        ~isempty(model_inputs.startexpo_bnds)
    startexpo_str = sprintf('%g-%g', ...
        model_inputs.startexpo_bnds(1), ...
        model_inputs.startexpo_bnds(2));
    save_str = strcat(save_str,'_startExpoBnds',startexpo_str);
end

% Optional snow cover
if isfield(model_inputs,'snow_depth_bnds') && ...
        ~isempty(model_inputs.snow_depth_bnds)
    snow_str = sprintf('%g-%g', ...
        model_inputs.snow_depth_bnds(1), ...
        model_inputs.snow_depth_bnds(2));
    save_str = strcat(save_str,'_snowBnds',snow_str);
end

% Optional subaerial erosion
if isfield(model_inputs,'subaerial_erosion_bnds') && ...
        ~isempty(model_inputs.subaerial_erosion_bnds)
    aero_str = sprintf('%g-%g', ...
        model_inputs.subaerial_erosion_bnds(1), ...
        model_inputs.subaerial_erosion_bnds(2));
    save_str = strcat(save_str,'_subaeroBnds',aero_str);
end

% Optional subglacial erosion
if isfield(model_inputs,'subglacial_erosion_bnds') && ...
        ~isempty(model_inputs.subglacial_erosion_bnds)
    gero_str = sprintf('%g-%g', ...
        model_inputs.subglacial_erosion_bnds(1), ...
        model_inputs.subglacial_erosion_bnds(2));
    save_str = strcat(save_str,'_subgeroBnds',gero_str);
end

% Optional burial fraction bounds
if isfield(model_inputs,'burialfrac_bnds') && ...
        ~isempty(model_inputs.burialfrac_bnds)
    burial_str = sprintf('%g-%g', ...
        model_inputs.burialfrac_bnds(1), ...
        model_inputs.burialfrac_bnds(2));
    save_str = strcat(save_str,'_burialFracBnds',burial_str);
end

% Final save
save(strcat(save_str,'.mat'),'MC_results');

disp(['Results saved to: ', save_str, '.mat'])

end