%% compute_brainmap_colorbar_bounds.m
%
% Gets statistics on group/individual level brainmaps. Takes a
% vertex_values file and computes things like the vertex average, vertex
% standard deviation, vertex variance, vertex maximum, vertex minimum

% Load vertex_values file of interest
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy 250ms\spectrogram_baseline5to0_common\Gamma_power_CP\vertex_values.mat')
% Compute vertex-wise sum across patients
vertex_values_sumL = squeeze(nansum(vertex_valuesL, 2));
divisorL = squeeze(sum(vertex_valuesL ~= 0, 2));
divisorL(divisorL == 0) = 1;
vertex_values_sumR = squeeze(nansum(vertex_valuesR, 2));
divisorR = squeeze(sum(vertex_valuesR ~= 0, 2));
divisorR(divisorR == 0) = 1;

% Compute scaled or normal averaging
scaled = 1;
if scaled
    vertex_values_sumL = vertex_values_sumL ./ sqrt(divisorL);
    vertex_values_sumR = vertex_values_sumR ./ sqrt(divisorR);
else
    vertex_values_sumL = vertex_values_sumL ./ divisorL;
    vertex_values_sumR = vertex_values_sumR ./ divisorR;
end

% Plot histogram of the vertex values
bins = 200;
figure
H = histogram([vertex_values_sumL; vertex_values_sumR], bins,...
    'Normalization', 'probability');
xlabel('Group Vertex Value (Scaled Average Gamma Power)')
ylabel('Normalized counts')
title('Mean Gamma Full Brain Distribution')

%% Compute the outer boundaries to retain some percentage of the data in the
% visualization
alpha = .05;
bounds = [1, bins];
loss = sum(H.Values([1 : bounds(1), bounds(2) : end]));
while loss < alpha
    upper_tail_loss = sum(H.Values(bounds(2) - 1 : end));
    lower_tail_loss = sum(H.Values(1 : bounds(1) + 1));
    if upper_tail_loss < lower_tail_loss
        bounds(2) = bounds(2) - 1;
    else
        bounds(1) = bounds(1) + 1;
    end
    loss = sum(H.Values([1 : bounds(1), bounds(2) : end]));
end
color_bar = H.BinEdges([bounds(1), bounds(2) + 1]); % a and b
tail_percentages = [sum(H.Values(1 : bounds(1))), sum(H.Values(bounds(2) : end))];

%% Compute the inner boundaries to retain some percentage of the data in the
% visualization
beta = .84;
[~, zero_init] = min(abs(H.BinEdges));
bounds = [zero_init, zero_init];
loss = sum(H.Values(bounds(1) : bounds(2)));
s = 1;
while loss < beta
    upper_trunk_loss = sum(H.Values(zero_init : bounds(2)));
    lower_trunk_loss = sum(H.Values(bounds(1) : zero_init));
    if s == 1
        
        bounds(2) = bounds(2) + 1;
        s = 0;
    else
        if bounds(1) > 1
            bounds(1) = bounds(1) - 1;
            s = 1;
        else
            s = 0;
        end
    end
    loss = sum(H.Values(bounds(1) : bounds(2)));
end
color_bar_inner = H.BinEdges([bounds(1), bounds(2) + 1]); % c and d
trunk_percentage = [sum(H.Values(bounds(1) : zero_init)), sum(H.Values(zero_init + 1 : bounds(2)))];

% create colorbar
mycmap = create_colorbar(color_bar, color_bar_inner);

% display colorbar
figure; colorbar; colormap(mycmap); caxis(color_bar)
