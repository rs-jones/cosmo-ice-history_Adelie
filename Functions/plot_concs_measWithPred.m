function plot_concs_measWithPred(results, samples_type, ...
    sample_data, nuclide, n_bestfits, save_name)

if strcmpi(samples_type,'core')
    these_bestfits='bestfits_core';
elseif strcmpi(samples_type,'transect')
    these_bestfits='bestfits_transect';
else
    error('must specify "core" or "transect"')
end

output_dir = fullfile(pwd,'Figures');
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

base = results.(these_bestfits);

% --- Select misfit
if isequal(sort(nuclide),[10 14])
    [~,bestfit_idx]=mink(base.rChi2_1014,n_bestfits);
elseif isequal(sort(nuclide),[10 26])
    [~,bestfit_idx]=mink(base.rChi2_1026,n_bestfits);
elseif isequal(nuclide,10)
    [~,bestfit_idx]=mink(base.rChi2_10,n_bestfits);
elseif isequal(nuclide,14)
    [~,bestfit_idx]=mink(base.rChi2_14,n_bestfits);
else
    [~,bestfit_idx]=mink(base.rChi2_26,n_bestfits);
end

predN_masked.N10 = base.predN_10(bestfit_idx,:);
predN_masked.N14 = base.predN_14(bestfit_idx,:);
predN_masked.N26 = base.predN_26(bestfit_idx,:);


% Filter measured nuclides to match input
sample_data_plot = sample_data;

plot10 = any(nuclide == 10);
plot14 = any(nuclide == 14);
plot26 = any(nuclide == 26);

for i = 1:length(sample_data_plot.s)
    if isfield(sample_data_plot.s{i},'nuclide10') && sample_data_plot.s{i}.nuclide10==1
        sample_data_plot.s{i}.nuclide10 = plot10;
    end
    if isfield(sample_data_plot.s{i},'nuclide14') && sample_data_plot.s{i}.nuclide14==1
        sample_data_plot.s{i}.nuclide14 = plot14;
    end
    if isfield(sample_data_plot.s{i},'nuclide26') && sample_data_plot.s{i}.nuclide26==1
        sample_data_plot.s{i}.nuclide26 = plot26;
    end
end


%% CORE
if strcmpi(samples_type,'core')

    fig_core_h = plot_core_concs(sample_data_plot,0);

    logical_all = [
        ~isnan(predN_masked.N10(1,:));
        ~isnan(predN_masked.N14(1,:));
        ~isnan(predN_masked.N26(1,:))
    ];

    plot_core_concs_pred(results.core_data,...
        fig_core_h,predN_masked,nuclide,logical_all,[]);

    if ~isempty(save_name)
        outfile = fullfile(output_dir, ...
            strcat(save_name,'_measWithPred_vsDepth.png'));
        saveas(gcf,outfile);
    end

%% TRANSECT
elseif strcmpi(samples_type,'transect')

    % ---- Existing 10Be–14C plots retained ----
    if all(ismember([10 14],nuclide))

        x_lim = [500,6e3];
        y_lim = [500,5e3];
        add_names = 1;

        fig_1014_h = plot_normconcs_1014( ...
            sample_data_plot,2,0,0,0,x_lim,y_lim,add_names);

        plot_concs_1014_pred(sample_data_plot, ...
            fig_1014_h,predN_masked,[],0);

        if ~isempty(save_name)
            outfile = fullfile(output_dir, ...
                strcat(save_name,'_measWithPred_1014.png'));
            saveas(gcf,outfile);
        end
    end

    % ---- 10Be–26Al (if plotting tools exist) ----
    if all(ismember([10 26],nuclide)) && ...
            exist('plot_normconcs_1026','file')

        fig_1026_h = plot_normconcs_1026(sample_data_plot);

        plot_concs_1026_pred(sample_data_plot, ...
            fig_1026_h,predN_masked);

        if ~isempty(save_name)
            outfile = fullfile(output_dir, ...
                strcat(save_name,'_measWithPred_1026.png'));
            saveas(gcf,outfile);
        end
    end
end

end