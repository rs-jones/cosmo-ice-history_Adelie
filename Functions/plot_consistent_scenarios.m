function plot_consistent_scenarios(results, field, save_name)

min_marker_size = 5;
max_marker_size = 30;

field_struct = split(field,'.');
field_group = results.(field_struct{1});
is_consistent = field_group.(field_struct{2});

has_snow = isfield(results,'rand_snow_depth');
if has_snow
    snow_depth = results.rand_snow_depth;
end

% --- Select correct misfit set
if contains(field,'core')
    base = results.bestfits_core;
elseif contains(field,'transect')
    base = results.bestfits_transect;
else
    error('Cannot identify misfit field.');
end

if contains(field,'1014')
    rChi2 = base.rChi2_1014;
elseif contains(field,'1026')
    rChi2 = base.rChi2_1026;
elseif contains(field,'10') && ~contains(field,'14')
    rChi2 = base.rChi2_10;
elseif contains(field,'14')
    rChi2 = base.rChi2_14;
elseif contains(field,'26')
    rChi2 = base.rChi2_26;
else
    error('Cannot identify rChi2 field.');
end

% --- Durations
if isfield(base,'time_exposed')
    exposure_duration = base.time_exposed/1000;
    burial_duration = base.time_buried/1000;
else
    % fallback for thinning-style runs
    exposure_duration = results.rand_thinning_end/1000;
    burial_duration = ...
        (results.rand_start_t - results.rand_thinning_end)/1000;
end

figure; clf;

nRows = 2;
nCols = has_snow + 1;
cmap = flipud(parula);

% Burial vs exposure
ax1 = subplot(nRows, nCols, 1);
hold on;

% Inconsistent
h_incon = scatter(burial_duration(~is_consistent), ...
                  exposure_duration(~is_consistent), ...
                  min_marker_size, [.6,.6,.6], 'o');

% Consistent
h_con = scatter(burial_duration(is_consistent), ...
                exposure_duration(is_consistent), ...
                max_marker_size, 'r', 'filled');

xlabel('Burial duration (kyr)');
ylabel('Exposure duration (kyr)');
grid on;

% Add legend outside axes
lg = legend([h_incon h_con], ...
    {'Inconsistent','Consistent'}, ...
    'Location','eastoutside');
lg.Box = 'off';

% Display message if none consistent
if ~any(is_consistent)
    disp('No modelled scenarios are consistent with measurements.');
end


% Misfit surface
ax2 = subplot(nRows,nCols,1+nCols);
scatter3(burial_duration, exposure_duration, rChi2, ...
         max_marker_size, rChi2,'filled');

xlabel('Burial duration (kyr)');
ylabel('Exposure duration (kyr)');
zlabel('Reduced \chi^2');
colormap(ax2,cmap);
colorbar;
set(gca,'ZDir','reverse');
grid on;

% Snow panels
if has_snow
    ax3 = subplot(nRows,nCols,2); hold on;
    scatter(snow_depth(~is_consistent), exposure_duration(~is_consistent), ...
        min_marker_size,[.6 .6 .6],'o');
    scatter(snow_depth(is_consistent), exposure_duration(is_consistent), ...
        max_marker_size,'r','filled');
    xlabel('Snow depth (m)');
    ylabel('Exposure duration (kyr)');
    grid on;

    ax4 = subplot(nRows,nCols,2+nCols);
    scatter3(snow_depth, exposure_duration, rChi2, ...
        max_marker_size, rChi2,'filled');
    xlabel('Snow depth (m)');
    ylabel('Exposure duration (kyr)');
    zlabel('Reduced \chi^2');
    colormap(ax4,cmap);
    colorbar;
    set(gca,'ZDir','reverse');
    grid on;
end

linkaxes([ax1 ax2],'xy');
if has_snow
    linkaxes([ax3 ax4],'xy');
end

%set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]); % Make full screen

if nargin > 2 && ~isempty(save_name)
    output_dir = fullfile(pwd,'Figures');

    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end

    outfile = fullfile(output_dir, ...
        strcat(save_name,'_scatter-consistent.png'));

    saveas(gcf,outfile);
end

end