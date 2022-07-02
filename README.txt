### ROADMAP

* Simulation for full model (SMA):
- exec: run_SMA4e
- analysis: analyze_ACCRT

* Simulation for reduced models:
- exec: run_SMA4e
- exec: run_sameseed (linear) or run_sameseed_quad (quadratic)
- analysis: analyze_full_trials

* Simulation for time-dependent linear reduced model:
- exec: run_agvals

* Simulation for WWM:
- exec: run_2MA

* Simulation for full scaled neural network model:
- fit: spontaneous_NMDA_group or spontaneous_NMDA_group_suppression
- fit: test_selectivtiy or test_selectivity_suppression

### POTENTIAL ISSUES
- file path problems
- Unaccounted for files: plot_accrt_selectivity, plot_accrt_thre_nsig, plot_accrt_thre_phase (?),
plot_accrt_thre-2, plot_accrt_thre-2-extreme, plot_aQ_by_bm, plot_aQ_by_bm-2, plot_phase_space_selectivity, plot_phase_space_slope,
plot_psych_metric_compare, plot_psych_metric_UMSM, plot_psych_metric-2, plot_psych_metric, plot_sameseed_compare_error, test*
- Files that should be added: plot_contour, plot_lag_as_derivative, plot_phase_space_accrt, plot_psych_metric_agvals,
plot_psych_metric_nsig, plot_psych_metric-fast, plot_reward_rate, plot_sameseed_compare_meanacc, plot_SNMDA_traces_2d,
plot_SNMDA_traces, plot_steady_state_demonstration, plot_trace_gif
