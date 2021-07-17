<!--Source code: https://github.com/othneildrew/Best-README-Template/edit/master/README.md -->

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/L24358/SDM">
    <img src="https://github.com/L24358/SDM/blob/main/graphs/SDM.PNG" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">SDM</h3>

  <p align="center">
    Code for <strong><em>Attractor Decision Network with Selective Inhibition</em></strong> by Belle (Pei-Hsien) Liu, Chung-Chuan Lo and Kuo-An Wu. <br/> bioRxiv link: ____
    <br />
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

- ``ss_flysim_iter=0``: runs flysim in current directory.
- ``ss_flysim_iter=1``: runs flysim in all first-level sub-directories.
- ``ss_flysim_iter=2``: runs flysim in all second-level sub-directories.

### analysis

- ``gen_pro``: generates protocol file.
- ``gen_conf``: generates configuration files.
- ``classes``: includes some commonly used classes. In particular, class ``motif`` generates small neuronal circuits when the ID is given. See how to specify circuit ID here:

### plot

### fit

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/L24358/SDM/issues) for a list of proposed features (and known issues).


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Belle (Pei-Hsien) Liu - belle.l24358@gmail.com

Project Link: [https://github.com/L24358/SDM](https://github.com/L24358/SDM)

