<!--Source code: https://github.com/othneildrew/Best-README-Template/edit/master/README.md -->

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/L24358/SDM">
    <img src="https://github.com/L24358/SDM/blob/main/graphs/SDM.PNG" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">SDM</h3>

  <p align="center">
    Code for <strong><em>Attractor Decision Network with Selective Inhibition</em></strong> by Belle (Pei-Hsien) Liu, Chung-Chuan Lo and Kuo-An Wu. <br/><br />
    <a href="https://github.com/L24358/SDM"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/L24358/SDM">View Demo</a>
    ·
    <a href="https://github.com/L24358/SDM/issues">Report Bug</a>
    ·
    <a href="https://github.com/L24358/SDM/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#exec">exec</a></li>
        <li><a href="#analysis">analysis</a></li>
        <li><a href="#plot">plot</a></li>
        <li><a href="#fit">fit</a></li>
      </ul>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

The ability to decide swiftly and accurately in an urgent scenario is crucial for an organism's survival. The neural mechanisms underlying the perceptual decision and trade-off between speed and accuracy have been extensively studied in the past few decades. Among several theoretical models, the attractor neural network model has successfully captured both behavioral and neuronal data observed in many decision experiments. However, a recent experimental study revealed additional details that were not considered in the original attractor model. In particular, the study shows that the inhibitory neurons in the posterior parietal cortex of mice are as selective to decision making results as the excitatory neurons, whereas the original attractor model assumes the inhibitory neurons to be unselective. In this study, we investigate a revised attractor model with selective inhibition, and analyze in detail the differences in computational ability that it brings. We proposed a reduced model for both the unselective and selective models, and showed that the selective model alone produces a time-varying energy landscape. This time dependence of the energy landscape allows the selective model to integrate information carefully in initial stages, then quickly converge to an attractor once the choice is clear. This results in the selective model having a more efficient speed-accuracy trade-off that is unreachable by unselective models.

<!-- USAGE EXAMPLES -->
## Usage

### exec

Code that runs simulation of full models, reduced models, time-dependent reduced models, and full-scale neural network models.
- ``run_2MA``: simulates WWM.
- ``run_full_scale_model_wsuppression``: simulates full scale models, with suppression of activity after decision is made.
- ``run_full_scale_model``: simulates full scale models, without suppression of activity after decision is made.
- ``run_LRM``: simulates linear reduced models.
- ``run_QRM``: simulates quadratic reduced models.
- ``run_SM_v4e``: simulates selective models (version 4e).
- ``run_TDLRM``: simulates linear time-dependent reduced models.

### analysis

Code that analyzes the simulation results, mainly extracting information regarding accuracy and reaction time.
- ``analyze_full_scale_model_selectivity_wsuppression``: analyzes full scale models with suppression.
- ``analyze_full_scale_model_selectivity``: analyzes full scale models without suppression.
- ``analyze_SM_parametrize_accrt``: parametrizes the speed and accuracy of full models.
- ``compare_SM_RM_accrt``: compares the accuracy and reaction time of full and reduced models.

### plot

Code that plots results.

- ``compare_ab_absim``: compares theoretical (noiseless) vs averaged (noisy) a, b and aQ, bQ values.
- ``compare_different_WWM_by_nullcline``: plot difference of nullcline between two WWM.
- ``compare_psychometric_wWWM``: plot psychometric function in the style of Wong and Wang's paper for comparison purposes.
- ``compare_SM_QRM_meanrtc``: compare the mean reaction time of the full model and the quadratic reduced model.
- ``compare_SM_RM_meanacc``: compare the accuracy of the full model, the linear reduced model and the quadratic reduced model.
- ``compare_SM_RM_meanrtc``: compare the mean reaction time of the full model, the linear reduced model and the quadratic reduced model.
- ``compare_SM_RM_trialbytrial``: plot the trial-by-trial comparison of the full model and reduced model.
- ``plot_abL``: fits aL(t), bL(t) values.
- ``plot_accrt_thre``: plots speed-accuracy trade-off plane for different thresholds.
- ``plot_aQ_by_bQ``: plot aQ as a function of bQ for different SMs.
- ``plot_bL_by_splus``: calculate crossing time as a function of splus, conditioned on the benchmark value (=0.05).
- ``plot_coherence_gif``: plot gif that illustrates the 2AFC experiment considered, and the concept of coherence.
- ``plot_dRT_by_threshold``: plot Delta_RT as a function of threshold.
- ``plot_lag_by_aQ_by_cinp``: plot the relation between a_Q, c_inp for the minus direction, and deviation (lag).
- ``plot_lag_by_derivative``: plot deviation (lag) as a function of the derivative of splus, sminus.
- ``plot_lag``: plot deviation (lag) as a function of time.
- ``plot_landscape_1D``: plot deviation (lag) as a function of the derivative of splus, sminus.
- ``plot_linearity_single_trial``: plot linear relation for a single trial of UM, SM plus; and seemingly linear relation of a single trial of SM minus.
- ``plot_log_by_log``: plot log-log relation of a single trial for SM, minus direction.
- ``plot_multiple_trials_SNMDA_1D``: plot multiple trials of S_NMDA in 1D.
- ``plot_multiple_trials_SNMDA_2D``: plot multiple trials of S_NMDA in 2D.
- ``plot_nullcine_favor1``: plot nullcline in which population 1 is favored.
- ``plot_nullcline_favor2``: plot nullcline in which population 2 is favored.
- ``plot_nullcline_winp``: plot nullcline with input.
- ``plot_phase_space_accrt``: plot accuracy and reaction time parameters in the phase space of gamma2 and gamma3.
- ``plot_psychometric_bootstrap_sample``: plot averaged psychometric function with single bootstrap sample.
- ``plot_psychometric_different_nsig``: plot psychometric function for models with different amount of noise (nsig).
- ``plot_psychometric_fast``: plot psychometric function of UM with faster GABA dynamics.
- ``plot_psychometric_LTDRM``: plot the psychometric function for LTDRM.
- ``plot_regrets_by_benchmark``: plot the percentage of regret trials as a function of benchmark.
- ``plot_reward_rate``: plot reward rate as a function of selectivity.
- ``plot_roc_demonstration``: demonstrative plot of response distribution and roc.
- ``plot_rt_cdf``: fit reaction time to gamma distribution.
- ``plot_selectivity_gif``: plot gif to explain the concept of selectivity.
- ``plot_single_trial_SNMDA_SGABA``: plot S_NMDA, S_GABA as a function of time.
- ``plot_single_trial_SNMDA``: plot S_NMDA as a function of time.
- ``plot_steady_state_demonstration``: plot demonstration of how the quadratic relation is a subtraction of deviation (lag) from the steady-state solution.
- ``plot_trace_gif``: plot traces of simulations in gif.
- ``plot_velocity_magnitude_contourplot``: plot velocity magnitude as a contour plot.
- ``plot_velocity_magnitude``: plot velocity magnitude as a function of time.
- ``sweep_WWM_v2``: plot parameter sweep of WWM as 2D slices.

### test

Code that tests rudimentary functions and ideas. Includes ``test_auc_method`` and ``test_psychmetric_std``.

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/L24358/SDM/issues) for a list of proposed features (and known issues).


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Belle (Pei-Hsien) Liu - belle.l24358@gmail.com

Project Link: [https://www.biorxiv.org/content/10.1101/2021.10.05.463257v2](https://www.biorxiv.org/content/10.1101/2021.10.05.463257v2)

