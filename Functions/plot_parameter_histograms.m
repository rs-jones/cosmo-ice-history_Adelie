function plot_parameter_histograms(results, field, save_name)

max_col = [0.85 0 0];
nbins = 20;

field_struct = split(field,'.');
field_group = results.(field_struct{1});
is_consistent = field_group.(field_struct{2});

if ~any(is_consistent)
    warning('No modelled scenarios are consistent with measurements - histograms empty.');
end

plot_fields = {'rand_start_t','rand_expo_duration'};
field_labels = {'Model start time (ka)', ...
                'Exposure duration (ka)'};
plot_data = {results.rand_start_t/1000, ...
             results.rand_expo_duration};

if isfield(results,'rand_snow_depth')
    plot_fields{end+1}='rand_snow_depth';
    field_labels{end+1}='Snow depth (m)';
    plot_data{end+1}=results.rand_snow_depth;
end

nPlots=numel(plot_fields);
nRows=ceil(sqrt(nPlots));
nCols=ceil(nPlots/nRows);

figure; clf;

for k=1:nPlots
    subplot(nRows,nCols,k); hold on;

    all_vals=plot_data{k};
    cons_vals=all_vals(is_consistent);

    bin_edges=linspace(min(all_vals),max(all_vals),nbins+1);

    histogram(cons_vals,bin_edges,...
        'FaceColor',max_col,...
        'EdgeColor','none',...
        'Normalization','probability');

    xlabel(field_labels{k});
    ylabel('Relative probability');
    grid on;
end

set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);

if nargin>2 && ~isempty(save_name)
    output_dir = fullfile(pwd,'Figures');

    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end

    outfile = fullfile(output_dir, ...
        strcat(save_name,'_parameter-histograms.png'));

    saveas(gcf,outfile);
end

end