%
% iso_plot = plot_normconcs_1014(sample_data,sigma,plot_ts)
% iso_plot = plot_normconcs_1014(sample_data,sigma,plot_ts,expo_intervals,bur_intervals)
% iso_plot = plot_normconcs_1014(sample_data,sigma,plot_ts,expo_intervals,bur_intervals,x_lim,y_lim,add_names)
%
% Plots normalised nuclide concentrations against each other (14C vs 10Be).
%
% The steady-state erosion island is plotted as a black line, while 
% exposure and burial isochrons are plotted as grey dashed and dot-dashed
% lines, respectively.
%
% sample_data is a required struct, created using get_data.m.
%
% sigma should be either 1 or 2. If '1' then only 1 sigma concentration
% uncertainty is plotted, and if '2' then both 1 and 2 sigma uncertainties
% are shown.
%
% plot_ts is a binary input for the plot type - '0' to plot just nuclide
% concentrations, or '1' to additionally plot the time-series data and
% corresponding periods of exposure/burial.
%
% expo_intervals and bur_intervals can be vectors of values in unit ka,
% which corresponds to the plotted and labelled isochrons. If empty, then
% defaults are used.
%
% x_lim and y_lim set the limits for the x-axis and y-axis ([lower upper]).
% If empty, then defaults are used.
%
% The figure (fig) and axes (ax) handles are exported as an output.
%
%
%%

function iso_plot = plot_normconcs_1014(sample_data,sigma,plot_ts,expo_intervals,bur_intervals,x_lim,y_lim,add_names)
  
  % Check inputs
  if (nargin < 3 || nargin > 8)
      error('plot_normconcs_1014 has wrong number of inputs!');
  end
  
  if (~sigma == 1 || ~sigma == 2)
      error('sigma should be 1 or 2!');
  end
  if plot_ts ~= 1 && plot_ts ~= 0
      error('plot_ts should be binary (1 or 0)!');
  end
  
  if (nargin > 3 && ~isempty(expo_intervals))
      expo_int = expo_intervals;
  else
      expo_int = [25,75,50,100,125]; % Default exposure intervals to be shown (ka)
  end
  if (nargin > 3 && ~isempty(bur_intervals))
      bur_int = bur_intervals;
  else
      bur_int = [10,20,30]; % Default burial intervals to be shown (ka)
  end
  if nargin < 8
      add_names = 0;
  end
  
  
  % Compute production rates for samples and normalise concentrations
  sample_data = norm_concs_1014(sample_data);
  
   
  % Get average parameters for samples
  for a = 1:length(sample_data.sp1026)
      s_rho(a) = sample_data.sp1026{a}.rb;
      s_L10(a) = sample_data.sp1026{a}.Lambdafe;
      s_L14(a) = sample_data.sp14{a}.Lambdafe;
  end
  rho = mean(s_rho);
  L10 = mean(s_L10);
  L14 = mean(s_L14);
  
  % Get average production rates
  P10_all = [];
  P14_all = [];
  for b = 1:length(sample_data.s)
      if isfield(sample_data.s{b}, 'PR_10')
          P10_all(end+1) = sample_data.s{b}.PR_10;
      end
      if isfield(sample_data.s{b}, 'PR_14')
          P14_all(end+1) = sample_data.s{b}.PR_14;
      end
  end
  P10 = mean(P10_all,'omitnan');
  P14 = mean(P14_all,'omitnan');
  
   
  % Continually exposed (non-eroding)
  A_10 = sample_data.pp.lambda10Be + rho / L10 * 0; % Zero erosion
  A_14 = sample_data.pp.lambda14C + rho / L14 * 0; % Zero erosion
  t_expos = linspace(1,10e6,10000+1); % Time exposed (a) - 1a to 10Ma
  
  be_conc_expo0ero = ((P10 ./ A_10) .* (1-exp(-1*A_10*t_expos))) ./ P10;
  c_conc_expo0ero = ((P14 ./ A_14) .* (1-exp(-1*A_14*t_expos))) ./ P14;
  
  
  % Continually exposed (steady-state erosion line)
  e_subaerial = linspace(1,0,20000); % Sub-aerial erosion (cm/a)
  A_10 = sample_data.pp.lambda10Be + rho / L10 * e_subaerial;
  A_14 = sample_data.pp.lambda14C + rho / L14 * e_subaerial;
  t_expos = 5e7; % Constant time exposed (50 Ma)
  
  be_conc_expoSSero = ((P10 ./ A_10) .* (1-exp(-1*A_10*t_expos))) ./ P10;
  c_conc_expoSSero = ((P14 ./ A_14) .* (1-exp(-1*A_14*t_expos))) ./ P14;
  
  
  % Burial isochrons
  A_10 = sample_data.pp.lambda10Be + rho / L10 * 0; % Zero erosion
  A_14 = sample_data.pp.lambda14C + rho / L14 * 0; % Zero erosion

  t_expos = linspace(0,10e6,10000+1); % Time exposed (a) - 1a to 10Ma
  t_expos(1) = 1; t_expos = [0,t_expos];
  t_bur = 0:1000:1e5; % Burial time (a) - 1a to 100 ka
  t_bur(1) = 1; t_bur = [0,t_bur];
  
  for k=1:length(t_expos) % Time consuming part
      be_conc(k) = (P10 ./ A_10) .* (1-exp(-1*A_10*t_expos(k)));
      c_conc(k) = (P14 ./ A_14) .* (1-exp(-1*A_14*t_expos(k)));
      for i=1:length(t_bur)
          be_conc_burial(k,i) = (be_conc(k)*exp(-1*A_10*t_bur(i))) ./ P10;
          c_conc_burial(k,i) = (c_conc(k)*exp(-1*A_14*t_bur(i))) ./ P14;
      end
  end
    
  
  % Determine isochron grid
  labels_expo0ero = strcat(num2str(expo_int'),' ka');
  expo0ero_col = [];
  for c = 1:length(expo_int)
      expo0ero_col(c) = find(t_expos==(expo_int(c)*1000));
  end
  labels_t_bur = strcat(num2str(bur_int'),' ka');
  burial_col = [];
  for c = 1:length(bur_int)
      burial_col(c) = find(t_bur==(bur_int(c)*1000));
  end
  
  
  % Sort nuclide data to plot; Mix if necessary
  plot_N = N_to_plot_1014(sample_data);
  
 
  % PLOT
  fig = figure;
  if plot_ts == 1
      a1 = subplot(3,1,1); a2 = subplot(3,1,2); a3 = subplot(3,1,3);
      set(a1,'pos',[0.14 0.40 0.78 0.54]);
      set(a2,'pos',[0.14 0.15 0.78 0.12]);
      set(a3,'pos',[0.14 0.11 0.78 0.03]);
  else
      a1 = subplot(1,1,1);
      a1.Position(2) = 0.13;
  end
  
  % Make colours
  grey = [.5 .5 .5];
  colour = [0.86,0.098,0.106; 0.98,0.703,0.68]; % Sample ellipses = Red
  %colour = [0.48,0.20,0.58; 0.76,0.64,0.81]; % Sample ellipses = Purple
  %colour = [0,0.53,0.216; 0.65,0.858,0.628]; % Sample ellipses = Green

  axes(a1);
  plot(be_conc_burial(:,burial_col),c_conc_burial(:,burial_col),'-.','Color',grey,'Linewidth',1.2);
  hold on;
  plot(be_conc_burial(expo0ero_col,:)',c_conc_burial(expo0ero_col,:)','--','Color',grey,'Linewidth',1.2);
  plot(be_conc_expo0ero,c_conc_expo0ero,'-k','Linewidth',1.4);
  plot(be_conc_expoSSero,c_conc_expoSSero,'-','Color','k','Linewidth',1.4);
  if ~expo_intervals == 0
      text(be_conc_expo0ero(1,expo0ero_col),c_conc_expo0ero(1,expo0ero_col),labels_expo0ero,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',grey)
      text(be_conc_burial(expo0ero_col(1)+1,burial_col),c_conc_burial(expo0ero_col(1)+1,burial_col),labels_t_bur,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',grey)
  end

  % Plot sample ellipses
  for el = 1:length(plot_N.norm_N10)
      N10 = plot_N.norm_N10(el);
      dN10 = plot_N.norm_dN10(el);
      N14 = plot_N.norm_N14(el);
      dN14 = plot_N.norm_dN14(el);
      if sigma == 1
          h_e = ellipse_concs(N10,dN10,N14,dN14,1);
          set(h_e,'Color',colour(1,:),'LineWidth',1,'LineStyle','-');
          plot(N10,N14,'.','Color',colour(1,:));
      else
          [h_e1,h_e2] = ellipse_concs(N10,dN10,N14,dN14,2);
          set(h_e1,'Color',colour(1,:),'LineWidth',1,'LineStyle','-');
          set(h_e2,'Color',colour(2,:),'LineWidth',1,'LineStyle','-');
          plot(N10,N14,'.','Color',colour(1,:));
      end
  end
 
  hold off;
  
  
  % Adjust axes
  if (nargin > 5 && ~isempty(y_lim))
      ylim(y_lim);
  else
      xlim([0,1.5e5]);
  end
  if (nargin > 5 && ~isempty(x_lim))
      xlim(x_lim);
  else
      xlim([0,1.5e5]);
  end
  %set(a1,'xscale','log');
  %xlabel('^{10}Be ^* (atoms g^{-1})')
  %xlabel('^{10}Be ^* (equivalent to years)')
  xlabel('^{10}Be ^*')
  %ylabel('^{14}C ^* (atoms g^{-1})')
  %ylabel('^{14}C ^* (equivalent to years)')
  ylabel('^{14}C ^*')
  
  box on;
  if expo_intervals == 0
      grid on;
  end
  
  % Add sample names to plot
  if add_names == 1
     add_names_concs1014(sample_data); 
  end
  
  % Export figure handles
  if plot_ts == 1
      iso_plot = [fig a1 a2 a3];
  else
      iso_plot = [fig a1];
  end
    

  %%%%%%%%%%%%%%%% Function to add sample names to the plot %%%%%%%%%%%%%%%
  
    function add_names_concs1014(norm_sample_data)

        if nnz(norm_sample_data.logical_10) ~= nnz(norm_sample_data.logical_14)

            warning('Names can only be added where samples have both Be-10 and C-14 measurements.')

        else

            % Get data for each sample
            names = cell(length(norm_sample_data.s),1);
            n10 = zeros(length(norm_sample_data.s),1);
            for n = 1:length(norm_sample_data.s)
                this_sample = norm_sample_data.s{n};
                names{n} = this_sample.name{1};
                n10(n) = this_sample.norm_N10;
                n14(n) = this_sample.norm_N14;
            end
            n10_5percent = (n10/100)*5;

            % Add text to plot
            h = text(n10+n10_5percent,n14,char(names));
            set(h,'FontSize',10);

        end

    end
  
end
