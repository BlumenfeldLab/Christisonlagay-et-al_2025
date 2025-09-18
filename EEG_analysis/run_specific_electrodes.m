%%%% grab thalamic contacts for timecourse:
%%% written by KChristisonLagay March 2021


mappeddatadrive='V';
relaventfile=[mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\451NH\meanpower_traces_77bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_from2to2_baseline5to0_commonr.mat']


load(relaventfile)

xtick_spacing_filename = [mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\filter_times_win100ms_shift25ms.mat'];
load (xtick_spacing_filename)


CPtraces=squeeze(meanpower_traces(43:45, 1, 3, :));

CNPtraces=squeeze(meanpower_traces(43:45, 2, 3, :));

figure; hold on

CPtraces_mean=nanmean(CPtraces,1);
CNPtraces_mean=nanmean(CNPtraces,1);

subtracted_mean=CPtraces_mean-CNPtraces_mean;

% plot(CPtraces(1,:), 'r')
% plot(CPtraces(2,:), 'b')
% plot(CPtraces(3,:), 'y')
plot(subtracted_mean, 'k')
plot(CPtraces_mean, 'r')
plot(CNPtraces_mean, 'b')


axis([19 77 -8 20])
times=(T-1)*1000;
xticks([19 29 39 49 59 69]);
xticklabels({'-500', '-250', '0', '250', '500', '750'});
xlabel('ms')
ylabel('weighted zscore')

xline(42, '--')
xline(39, '--')

