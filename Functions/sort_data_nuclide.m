%
% sorted_data = sort_data_nuclide(in_data, N, names, nuclide_tag)
%
% Sorts in_data if multiple samples were combined for nuclide measurements.
% Duplicate measurements are identified when both nuclide concentration
% and 1-sigma uncertainty are identical.
%
% If such duplicates are found, a new data matrix is produced where each
% new row represents a unique nuclide measurement. Individual sample
% properties (depths, weights) are retained in cell arrays.
%
% INPUTS
%   in_data : numeric matrix of sample data (same format as get_data_101426)
%   N       : [N_i, dN_i] matrix for the nuclide being sorted
%   names   : cell array of sample names
%   nuclide_tag : string identifier: '10', '14', or '26'
%
% OUTPUT
%   sorted_data : struct with fields
%       .data
%       .names
%       .duplicates
%       .dup_names
%       .dup_top_z
%       .dup_bottom_z
%       .dup_weight
%
% Written by Richard Selwyn Jones
%
%%

function sorted_data = sort_data_nuclide(in_data, N, names, nuclide_tag)

  % Define columns
  col_topz = 9;
  col_bottz = 10;
  if strcmpi(nuclide_tag,'Be10')
      col_weight = 14;
  elseif strcmpi(nuclide_tag,'C14')
      col_weight = 20;
  elseif strcmpi(nuclide_tag,'Al26')
      col_weight = 26;
  end

  % Identify duplicate means
  [~,ia,ic] = unique(N(:,1),'stable');
  dup_N = ismember(N(:,1), N(ia(accumarray(ic,1)>1),1));

  % Identify duplicate uncertainties
  [~,ia,ic] = unique(N(:,2),'stable');
  dup_dN = ismember(N(:,2), N(ia(accumarray(ic,1)>1),2));

  true_dup = dup_N & dup_dN;

  if any(true_dup)

      in_data_N = in_data(~true_dup,:);
      in_names_N = names(~true_dup)';

      dup_ids = unique(find(true_dup));
      n_dup = numel(dup_ids);

      sorted = zeros(n_dup,size(in_data,2));
      sorted_names = cell(n_dup,1);

      for k = 1:n_dup
          ind = find(true_dup);
          ind = ind(k);

          sorted_names{k} = names{ind};

          sorted(k,:) = in_data(ind,:);

          % Retain per-sample values
          sorted_data.dup_names{k} = names(ind);
          sorted_data.dup_top_z{k} = in_data(ind,col_topz)';
          sorted_data.dup_bottom_z{k} = in_data(ind,col_bottz)';
          sorted_data.(['dup_weight_' nuclide_tag]){k} = in_data(ind,col_weight)';
      end

      sorted_data.names = [in_names_N; sorted_names];
      sorted_data.data = [in_data_N; sorted];
      sorted_data.duplicates = [false(1,size(in_data_N,1)) true(1,n_dup)];

  else
      sorted_data.names = names(~true_dup)';
      sorted_data.data = in_data(~true_dup,:);
      sorted_data.duplicates = false(1,size(sorted_data.data,1));
  end

end