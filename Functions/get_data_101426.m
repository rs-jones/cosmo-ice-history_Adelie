%
% sample_data = get_data_101426(input_name)
%
% Loads measured cosmogenic nuclide data from an input .xlsx file,
% extracts sample metadata and nuclide concentrations, sorts and combines
% samples with multiple nuclide measurements, and exports all information
% to a single structured output suitable for modelling.
%
% This function supports arbitrary combinations of in situ cosmogenic
% 10Be, 14C, and/or 26Al measured in each sample. Samples may contain one,
% two, or all three nuclides. Missing nuclide data must be entered as NaN
% in the spreadsheet and are automatically ignored.
%
% The function also prepares nuclide-specific input structures compatible
% with CRONUScalc-style production rate calculations, and sorts data for
% plotting and multi-nuclide comparison.
%
% For CRONUScalc calculations, erosion rates and inheritance are assumed
% to be zero. Attenuation lengths are determined automatically using the
% Sato scaling scheme, based on the geographic location of each sample.
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
% -------------------------------------------------------------------------
% NUCLIDE-SPECIFIC DATA BLOCKS (NaN if not measured):
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
% -------------------------------------------------------------------------
% OPTIONAL METADATA:
%
% 31. Scaling model used (e.g. 'DE','DU','LI','ST','LM','LSD','LSDn').
%     If provided, the same scaling model is assumed for all samples.
%
%
% Exposure ages are not required for forward modelling, but may be used
% by specific model configurations for comparison, calibration, or
% sensitivity analyses.
%
% -------------------------------------------------------------------------
% OUTPUT:
%
% sample_data – structure containing:
%   - sample-level metadata
%   - nuclide-specific concentrations, uncertainties, and mineral weights
%   - logical flags indicating which nuclides are present
%   - depth information in cm and g cm^-2
%   - CRONUScalc-compatible input arrays for each nuclide
%
% -------------------------------------------------------------------------
% Written by Richard Selwyn Jones
% Monash University
% richard.s.jones@monash.edu
%
% Associated with the iceTEA tools suite, which is built on CRONUS-Earth 
% and CRONUScalc code.
%
%%

function sample_data = get_data_101426(input_name)

disp('Loading measured nuclide data...');

% Load data
[in_data,in_txt,~] = xlsread([input_name '.xlsx']);

n_cols = size(in_txt,2);

% Input checks
if n_cols > 31
    error('sample input data has too many fields!');
end
if n_cols > 18 && n_cols < 31
    error('sample input data seems to have exposure ages but fields are missing');
end

% Define columns for text fields
col_names = 1;
col_scaling = 31;

% Define columns for data fields
col_latlon = 1:2;
col_elevpress = 3:4;
col_pos = 5;
col_thk = 6;
col_dens = 7;
col_shield = 8;
col_depths = 9:10;
col_year = 11;
col_10concs = 12:13;
col_10weight = 14;
col_10ages = 15:17;
col_14concs = 18:19;
col_14weight = 20;
col_14ages = 21:23;
col_26concs = 24:25;
col_26weight = 26;
col_26ages = 27:29;


n_samples = size(in_data,1);
names = in_txt(2:end,col_names)';

% Replace NaN nuclide concentrations with zero
for i = 1:n_samples
    if isnan(in_data(i,col_10concs(1))), in_data(i,col_10concs) = 0; end % 10Be
    if isnan(in_data(i,col_14concs(1))), in_data(i,col_14concs) = 0; end % 14C
    if isnan(in_data(i,col_26concs(1))), in_data(i,col_26concs) = 0; end % 26Al
end

% Determine thickness or sample depths if missing
for i = 1:n_samples
    if isempty(in_data(i,col_thk)) || in_data(i,col_thk)==0 || isnan(in_data(i,col_thk))
        in_data(i,col_thk) = in_data(i,col_depths(2)) - in_data(i,col_depths(1));
    end
    if isempty(in_data(i,col_depths(1))) || isnan(in_data(i,col_depths(1)))
        in_data(i,col_depths(1)) = in_data(i,col_depths(2)) - in_data(i,col_thk);
    end
    if isempty(in_data(i,col_depths(2))) || isnan(in_data(i,col_depths(2)))
        in_data(i,col_depths(2)) = in_data(i,col_depths(1)) + in_data(i,col_thk);
    end
end

% Logical flags for nuclides
logical_10 = any(in_data(:,col_10concs(1)),2)';
logical_14 = any(in_data(:,col_14concs(1)),2)';
logical_26 = any(in_data(:,col_26concs(1)),2)';
NN = [any(logical_10),any(logical_14),any(logical_26)];

% Sort per-nuclide (duplicate AMS handling)
if any(logical_10)
    sorted_10 = sort_data_nuclide(in_data(logical_10,:), ...
        in_data(logical_10,col_10concs), names(logical_10),'Be10');
else, sorted_10 = [];
end

if any(logical_14)
    sorted_14 = sort_data_nuclide(in_data(logical_14,:), ...
        in_data(logical_14,col_14concs), names(logical_14),'C14');
else, sorted_14 = [];
end

if any(logical_26)
    sorted_26 = sort_data_nuclide(in_data(logical_26,:), ...
        in_data(logical_26,col_26concs), names(logical_26),'Al26');
else, sorted_26 = [];
end

% Merge datasets
logical_1014 = logical_10 & logical_14;
logical_1026 = logical_10 & logical_26;
%logical_101426 = logical_10 & logical_14 & logical_26;

logical_any_multi = any(logical_1014 | logical_1026);
logical_all_multi = all(logical_1014 | logical_1026);

if logical_all_multi % All samples have 10Be + at least one other nuclide

    sorted_names = sorted_10.names;
    sorted_data  = sorted_10.data;
    duplicate_logical = sorted_10.duplicates;
    duplicates_10 = sorted_10.duplicates;
    if NN(2)
        duplicates_14 = sorted_14.duplicates;
    else
        duplicates_14 = zeros(size(duplicates_10));
    end
    if NN(3)
        duplicates_26 = sorted_26.duplicates;
    else
        duplicates_26 = zeros(size(duplicates_10));
    end

elseif logical_any_multi % Some samples have multiple nuclides

    comb_names = sorted_10.names;
    if NN(2)
        comb_names = [comb_names; sorted_14.names];
    end
    if NN(3)
        comb_names = [comb_names; sorted_26.names];
    end
    [sorted_names,~,~] = unique(comb_names,'stable');
    sorted_data = sorted_10.data;
    duplicate_logical = sorted_10.duplicates;
    duplicates_10 = sorted_10.duplicates;

    if NN(2) % Append 14C-only samples
        extra_14 = ~ismember(sorted_14.names, sorted_10.names);
        sorted_data = [sorted_data; sorted_14.data(extra_14,:)];
        duplicate_logical = [duplicate_logical sorted_14.duplicates(extra_14)];
        duplicates_14 = [zeros(1,length(sorted_10.names)) sorted_14.duplicates(extra_14)];
    end
    if NN(3) % Append 26Al-only samples
        extra_26 = ~ismember(sorted_26.names, sorted_10.names);
        sorted_data = [sorted_data; sorted_26.data(extra_26,:)];
        duplicate_logical = [duplicate_logical sorted_26.duplicates(extra_26)];
        duplicates_26 = [zeros(1,length(sorted_10.names)) sorted_26.duplicates(extra_26)];
    end

% Only one nuclide exists
elseif NN(1) && ~NN(2) && ~NN(3)
    sorted_names = sorted_10.names;
    sorted_data  = sorted_10.data;
    duplicate_logical = sorted_10.duplicates;
    duplicates_10 = sorted_10.duplicates;

elseif ~NN(1) && NN(2) && ~NN(3)
    sorted_names = sorted_14.names;
    sorted_data  = sorted_14.data;
    duplicate_logical = sorted_14.duplicates;
    duplicates_14 = sorted_14.duplicates;

elseif ~NN(1) && ~NN(2) && NN(3)
    sorted_names = sorted_26.names;
    sorted_data  = sorted_26.data;
    duplicate_logical = sorted_26.duplicates;
    duplicates_26 = sorted_26.duplicates;
    
else % No sample has multiple nuclides

    sorted_names = [];
    sorted_data = [];
    duplicate_logical = [];
    if NN(1)
        sorted_names = [sorted_names; sorted_10.names];
        sorted_data  = [sorted_data;  sorted_10.data];
        duplicates_10 = sorted_10.duplicates;
    end
    if NN(2)
        sorted_names = [sorted_names; sorted_14.names];
        sorted_data  = [sorted_data;  sorted_14.data];
        duplicates_14 = sorted_14.duplicates;
    end
    if NN(3)
        sorted_names = [sorted_names; sorted_26.names];
        sorted_data  = [sorted_data;  sorted_26.data];
        duplicates_26 = sorted_26.duplicates;
    end
end


% Calculate number of sorted samples
n_sorted = size(sorted_data,1);

logical_sorted_10 = any(sorted_data(:,col_10concs(1)), 2)';
logical_sorted_14 = any(sorted_data(:,col_14concs(1)), 2)';
logical_sorted_26 = any(sorted_data(:,col_26concs(1)), 2)';


% Collate sample details
for c = 1:n_sorted

    sample_data.s{c}.name = sorted_names(c);

    % ---- 10Be ----
    if logical_sorted_10(c)
        sample_data.s{c}.nuclide10 = 1;
        sample_data.s{c}.N10 = sorted_data(c,col_10concs(1));
        sample_data.s{c}.dN10 = sorted_data(c,col_10concs(2));
        sample_data.s{c}.weight10 = sorted_data(c,col_10weight);
    else
        sample_data.s{c}.nuclide10 = 0;
    end

    % ---- 14C ----
    if logical_sorted_14(c)
        sample_data.s{c}.nuclide14 = 1;
        sample_data.s{c}.N14 = sorted_data(c,col_14concs(1));
        sample_data.s{c}.dN14 = sorted_data(c,col_14concs(2));
        sample_data.s{c}.weight14 = sorted_data(c,col_14weight);
    else
        sample_data.s{c}.nuclide14 = 0;
    end

    % ---- 26Al ----
    if logical_sorted_26(c)
        sample_data.s{c}.nuclide26 = 1;
        sample_data.s{c}.N26 = sorted_data(c,col_26concs(1));
        sample_data.s{c}.dN26 = sorted_data(c,col_26concs(2));
        sample_data.s{c}.weight26 = sorted_data(c,col_26weight);
    else
        sample_data.s{c}.nuclide26 = 0;
    end

    % ---- Geometry & density (nuclide-independent) ----
    sample_data.s{c}.top_z    = sorted_data(c,col_depths(1));
    sample_data.s{c}.bottom_z = sorted_data(c,col_depths(2));
    sample_data.s{c}.density  = sorted_data(c,col_dens);

end

% Add all details for each nuclide measurement (if samples were combined)
if any(duplicate_logical)

    % ---- 10Be duplicates ----
    if exist('duplicates_10','var') && any(duplicates_10)
        dup_s = find(duplicates_10);
        for d = 1:length(dup_s)
            s = dup_s(d);
            sample_data.s{s}.name      = sorted_10.dup_names{d};
            sample_data.s{s}.top_z     = sorted_10.dup_top_z{d};
            sample_data.s{s}.bottom_z  = sorted_10.dup_bottom_z{d};
            sample_data.s{s}.weight10  = sorted_10.dup_weight_Be10{d};
        end
    end

    % ---- 14C duplicates ----
    if exist('duplicates_14','var') && any(duplicates_14)
        dup_s = find(duplicates_14);
        for d = 1:length(dup_s)
            s = dup_s(d);
            sample_data.s{s}.name      = sorted_14.dup_names{d};
            sample_data.s{s}.top_z     = sorted_14.dup_top_z{d};
            sample_data.s{s}.bottom_z  = sorted_14.dup_bottom_z{d};
            sample_data.s{s}.weight14  = sorted_14.dup_weight_C14{d};
        end
    end

    % ---- 26Al duplicates ----
    if exist('duplicates_26','var') && any(duplicates_26)
        dup_s = find(duplicates_26);
        for d = 1:length(dup_s)
            s = dup_s(d);
            sample_data.s{s}.name      = sorted_26.dup_names{d};
            sample_data.s{s}.top_z     = sorted_26.dup_top_z{d};
            sample_data.s{s}.bottom_z  = sorted_26.dup_bottom_z{d};
            sample_data.s{s}.weight26  = sorted_26.dup_weight_Al26{d};
        end
    end
end

% Calculate mass depths
for e = 1:n_sorted
    sample_data.s{e}.top_z_gcm2 = sample_data.s{e}.top_z .* sorted_data(e,col_dens);
    sample_data.s{e}.bottom_z_gcm2 = sample_data.s{e}.bottom_z .* sorted_data(e,col_dens);
end


% CC matrices (CRONUScalc-compatible)

% Assume certain sample details are unknown or zero
e_rate = zeros(n_sorted,1);   % Erosion rate (mm/kyr)
inh10  = zeros(n_sorted,1);   % 10Be inheritance
inh14  = zeros(n_sorted,1);   % 14C inheritance
inh26  = zeros(n_sorted,1);   % 26Al inheritance

% Initially set attenuation length to zero
init_L = zeros(n_sorted,1);


if any(logical_sorted_10) && any(logical_sorted_26)

    % Build Be10–Al26 matrix
    sample_data.CC.Be10Al26 = ...
        [ sorted_data(:,[col_latlon,col_elevpress,col_thk,col_dens,col_shield]), ...
          e_rate, ...
          sorted_data(:,col_10concs(1)), ...   % Be10
          sorted_data(:,col_26concs(1)), ...   % Al26
          inh10, inh26, ...
          init_L, ...
          sorted_data(:,col_depths(1)), ...
          sorted_data(:,col_year) ];    % year

    % Determine atmospheric pressure
    sample_data.CC.Be10Al26 = atm_pressure( ...
        sorted_data(:,1), sorted_data(:,2), sample_data.CC.Be10Al26);

    % Determine attenuation lengths (Sato)
    for d = 1:n_sorted
        L(d,1) = attenuationlength( ...
            sample_data.CC.Be10Al26(d,1), ...
            sample_data.CC.Be10Al26(d,2), ...
            sample_data.CC.Be10Al26(d,3), ...
            sample_data.CC.Be10Al26(d,4));
    end
    sample_data.CC.Be10Al26(:,13) = L;

else
    % Build Be10-only matrix
    sample_data.CC.Be10Al26 = ...
        [ sorted_data(:,[col_latlon,col_elevpress,col_thk,col_dens,col_shield]), ...
          e_rate, ...
          sorted_data(:,col_10concs(1)), ...
          zeros(n_sorted,1), ...   % dummy Al column
          inh10, inh14, ...
          init_L, ...
          sorted_data(:,col_depths(1)), ...
          sorted_data(:,col_year) ];

    sample_data.CC.Be10Al26 = atm_pressure( ...
        sorted_data(:,1), sorted_data(:,2), sample_data.CC.Be10Al26);

    for d = 1:n_sorted
        L10(d,1) = attenuationlength( ...
            sample_data.CC.Be10Al26(d,1), ...
            sample_data.CC.Be10Al26(d,2), ...
            sample_data.CC.Be10Al26(d,3), ...
            sample_data.CC.Be10Al26(d,4));
    end
    sample_data.CC.Be10Al26(:,13) = L10;
end

if any(logical_sorted_14)

    % Build C14 matrix
    sample_data.CC.C14 = ...
        [ sorted_data(:,[col_latlon,col_elevpress,col_thk,col_dens,col_shield]), ...
          e_rate, ...
          sorted_data(:,col_14concs(1)), ...
          inh14, ...
          init_L, ...
          sorted_data(:,col_depths(1)), ...
          sorted_data(:,col_year) ];

    sample_data.CC.C14 = atm_pressure( ...
        sorted_data(:,1), sorted_data(:,2), sample_data.CC.C14);

    % Fixed attenuation length for C-14
    sample_data.CC.C14(:,11) = 140;
end


% Uncertainty matrices

uncert_1026 = zeros(size(sample_data.CC.Be10Al26));
if any(logical_sorted_10) && any(logical_sorted_26)    
    uncert_1026(:,9)  = sorted_data(:,col_10concs(2)); % Be uncertainty
    uncert_1026(:,10) = sorted_data(:,col_26concs(2)); % Al uncertainty
    uncert_1026(:,3)  = 5; % Elevation - use 5 m as default uncertainty
    sample_data.CC_uncert.Be10Al26 = uncert_1026;
else
    uncert_1026(:,9) = sorted_data(:,col_10concs(2)); % Be uncertainty
    uncert_1026(:,3) = 5; % Elevation - use 5 m as default uncertainty
    sample_data.CC_uncert.Be10Al26 = uncert_1026;
end
if isfield(sample_data.CC,'C14')
    uncert_14 = zeros(size(sample_data.CC.C14));
    uncert_14(:,9) = sorted_data(:,col_14concs(2));
    uncert_14(:,3) = 5; % Elevation - use 5 m as default uncertainty
    sample_data.CC_uncert.C14 = uncert_14;
end


% Ages struct (if ages are present)
if any(~isnan(sorted_data(:,col_10ages(1))))
    sample_data.ages.Be10 = sorted_data(:,col_10ages);
end
if any(~isnan(sorted_data(:,col_14ages(1))))
    sample_data.ages.C14  = sorted_data(:,col_14ages);
end
if size(in_data,2) == n_cols && any(~isnan(sorted_data(:,col_26ages(1))))
    sample_data.ages.Al26 = sorted_data(:,col_26ages);
end
if size(in_data,2) == n_cols && (any(~isnan(sorted_data(:,col_10ages(1)))) || any(~isnan(sorted_data(:,col_14ages(1)))) || any(~isnan(sorted_data(:,col_26ages(1)))))
    sample_data.ages.scaling = in_txt(2:end,col_scaling);
end

% Sample position
sample_data.position = sorted_data(:,col_pos)';

% Default cover depth
sample_data.cover.z = 0;

% Export nuclide logicals
sample_data.logical_10 = logical_10;
sample_data.logical_14 = logical_14;
sample_data.logical_26 = logical_26;

% Export sample names
sample_data.names = sorted_names;

disp('done.');

end
