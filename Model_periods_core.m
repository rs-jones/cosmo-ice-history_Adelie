%% SCRIPT TO MODEL MULTIPLE NUCLIDE CONCENTRATIONS FROM SURFACE SAMPLES USING DYNAMIC PERIODS OF EXPOSURE AND BURIAL
%
% Performs a suite of Monte Carlo simulations to determine the best-fit
% exposure/burial scenario to explain the measured nuclide concentations. 
% Plots the corresponding history, and predicted concentrations on a 
% two-isotope diagram.
%
% Uses the floating parameters of total exposure-burial time and exposure
% period duration (in years) to compute a multi-stage exposure-burial 
% history.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT DATA (one row per rock sample):
%
%  1. Sample name
%  2. Latitude (decimal degrees)
%  3. Longitude (decimal degrees)
%  4. Elevation (m a.s.l.)
%  5. Atmospheric pressure (hPa; set to 0 if unknown)
%  6. Relative position (e.g. distance from ice margin or elevation above ice)
%  7. Sample thickness (cm)
%  8. Bulk density (g cm^-3)
%  9. Shielding factor (terrain, snow, etc.; unitless)
% 10. Top depth of sample (cm)
% 11. Bottom depth of sample (cm)
% 12. Year the sample was collected (calendar year)
%
% ---
% Nuclide-specific data blocks (NaN if not measured):
%
% 10Be:
% 13. 10Be concentration (atoms g^-1)
% 14. 10Be concentration uncertainty (1 sigma; atoms g^-1)
% 15. 10Be final mineral weight (g)
% 16. 10Be exposure age (mean; years) [optional]
% 17. 10Be exposure age uncertainty (internal 1 sigma; years) [optional]
% 18. 10Be exposure age uncertainty (external 1 sigma; years) [optional]
%
% 14C:
% 19. 14C concentration (atoms g^-1)
% 20. 14C concentration uncertainty (1 sigma; atoms g^-1)
% 21. 14C final mineral weight (g)
% 22. 14C exposure age (mean; years) [optional]
% 23. 14C exposure age uncertainty (internal 1 sigma; years) [optional]
% 24. 14C exposure age uncertainty (external 1 sigma; years) [optional]
%
% 26Al:
% 25. 26Al concentration (atoms g^-1)
% 26. 26Al concentration uncertainty (1 sigma; atoms g^-1)
% 27. 26Al final mineral weight (g)
% 28. 26Al exposure age (mean; years) [optional]
% 29. 26Al exposure age uncertainty (internal 1 sigma; years) [optional]
% 30. 26Al exposure age uncertainty (external 1 sigma; years) [optional]
%
% ---
% Optional metadata:
%
% 31. Scaling model used (e.g. 'DE','DU','LI','ST','LM','LSD','LSDn').
%     If provided, the same scaling model is assumed for all samples.
%
% Exposure ages are not required for forward modelling, but may be used
% by specific model configurations for comparison, calibration, or
% sensitivity analyses.
%
%
%% Initialise model workspace (USER ACTION REQUIRED)

clear % Start fresh
addpath(genpath(pwd))

% SET inputs
input_name = 'exampleCoreData_101426'; % Name used for sample data .xlsx
scaling_model = 'LSDn'; % Scaling model - 'ST','LM','LSD','LSDn'
model_interval = 100; % Optionally set model interval used in calculations - 10, 100 or 1000 years (default is 1000)


%% Get Data

% Load sample data
core_data = get_data_101426(input_name);

% Get parameters for samples
prod_max_depth = 1000; % OPTIONAL USER MODIFICATION: Provide maximum depth (m) if considering production at large depths (e.g. for subglacial erosion)
    % leave empty to use default; increase depth if model gives the error "Prodz called for depth greater than maxdepth."

core_data = get_pars(core_data,scaling_model,prod_max_depth*100*2.7);


% Generate figure of core concentrations and plot handles
close; % Close existing figure, if open
plot_core_concs(core_data,0);


%% Model inputs (USER ACTION REQUIRED)

inputsMC.nuclides = ["10Be","14C","26Al"]; % Specify nuclides to model
inputsMC.recent_expo_sample = []; % Specify ID of sample to use as the age of most recent exposure period (e.g. "Sample1-10") - leave empty to exclude

inputsMC.expo_mids = [4,8,15,22,25,28,33,35,44,53,67,86,92,104,165]; % Specify mid-points of exposure periods (ka before present)

inputsMC.start_time_bnds = [20 100]; % Set start time for the model in ka before present (lower and upper)
inputsMC.expo_dur_bnds = [1 15]; % Set exposure time duration during periods in kyr (lower and upper)

%inputsMC.burialfrac_bnds = [0 .9]; % Fraction of exposure period that a sample could be buried (lower and upper)
inputsMC.startexpo_bnds = []; % Optionally set exposure duration at the start of the model, prior to the dynamic periods (kyr; lower and upper) - leave empty to exclude

inputsMC.snow_depth_bnds = []; % Set surface cover depth in metres (lower and upper) - leave empty to exclude
inputsMC.subaerial_erosion_bnds = []; % Set subaerial erosion rate in mm/ka (lower and upper) - leave empty to exclude - see Marrero et al. (2018)
inputsMC.subglacial_erosion_bnds = []; % Set subglacial erosion rate in mm/ka (lower and upper) - leave empty to exclude - see Hallet et al. (1996)

MC_iterations = 50; % Set number of Monte Carlo iterations - initially test with just a few iterations, and increase once happy with setup

plotfig = true; % Set true to plot a figure (updated each model iteration); otherwise, set false.


%% Run model

close all % Close any existing figures in the workspace
MC_results = run_periods("core",core_data,[],inputsMC,model_interval,MC_iterations,plotfig); % Perform Monte Carlo analysis


% Save output
save_name = 'MC_example-core_dynamicPeriods_bestfits'; % OPTIONAL USER MODIFICATION: Set base save name
save_periods_results(save_name,inputsMC,MC_iterations,MC_results);


%% Process results 

save_name = 'MC_example-core_dynamicPeriods_bestfits'; % OPTIONAL USER MODIFICATION: Set base save name

close all % Close any existing figures in the workspace


% Plot all scenarios where the predicted concentrations are consistent with measurements
consistentField_toPlot = 'bestfits_core.consistent_1026'; % USER ACTION REQUIRED: Specify consistent results field to plot 
    % (see MC_results struct) - e.g. consistent_1014, consistent_1026
plot_consistent_scenarios(MC_results,consistentField_toPlot,save_name);


% Plot histograms of the parameters
consistentField_toPlot = 'bestfits_core.consistent_1026'; % USER ACTION REQUIRED: Specify consistent results field to plot 
    % (see MC_results struct) - e.g. consistent_1014, consistent_1026
plot_parameter_histograms(MC_results,consistentField_toPlot,save_name);


% Plot measured with predicted concentrations 
nuclides_toPlot = [10,26]; % USER ACTION REQUIRED: Set nuclides to plot - 10=Be10, 14=C14, 26=Al26
n_bestfits = 5; % USER ACTION REQUIRED: Specify number of bestfits to plot (lowest misfit predictions)
plot_concs_measWithPred(MC_results,'core',core_data,nuclides_toPlot,n_bestfits,save_name);
