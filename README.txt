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
- Move full scaled neural network model into exec/analysis to be consistent