# ASSC26 Tutorial

This repository hosts materials for the "A practical guide to EEG analysis tools used in the neuroscience of consciousness and cognition" tutorial, for [ASSC26](https://theassc.org/assc-26/#tutorials).

The slides for the course are below:

[Intro](https://docs.google.com/presentation/d/1X7yLfVzOVmVLF_CH6Phm_tp9EJp2D_Q4st-tBgOjoY8/edit?usp=sharing)

[Synchrony](https://docs.google.com/presentation/d/1NUFq5qjLdXKk3nqVy6pAoQ2IAdzKaC3PZqalbEzcrqg/edit?usp=sharing)

[NESS2023](https://docs.google.com/presentation/d/1WpJl8nP0swUqbdUwYPvmXec1HPUf0hXI2juuXgZ0AU0/edit?usp=sharing)

## Overview

In the neuroscience of consciousness, electroencephalography (EEG) is an invaluable tool to non-invasively study human brain function. Indeed, the EEG allows the study of mechanisms that may hamper specific brain functions – for example, slow waves in the awake EEG are often linked to pathological dynamics. EEG is currently used as a clinical tool to monitor levels of consciousness during anesthesia and in patients with disorders of consciousness. In the process of using EEG to understand conscious states, methodological choices in data analysis figure heavily in the development of scientific theories and vice versa. With recent developments in methods for analyzing EEG data, it is essential that practitioners appropriately apply analysis tools.

In this tutorial, the presenters will provide a hands-on introduction to key methodological concepts and popular tools used in the analysis of EEG brain signals. Working from interactive notebooks, this tutorial is designed as an approachable guide to neural data analysis for both investigators and trainees. Attendees will also receive a brief introduction to each method’s statistical and theoretical underpinnings. We will provide an overview of EEG data analysis techniques that target biomarkers of consciousness, grouped into four categories:

1. Slow waves: EEG slow waves are a common feature of all unconscious states, and are thought to reflect up- and down-states in the spiking of cortical neurons. Recent research suggests that different levels of unconsciousness are differentiated by the entrainment of high-frequency activity to the slow wave in frontal vs posterior EEG. This tutorial will cover time-domain methods for decomposing the EEG signal into slow waves and high-frequency activity, and relating them to each other (Hilbert-based methods, State space oscillators, phase-amplitude coupling).

2. Aperiodic neural activity: States of consciousness dramatically differ in overall signal power. This can be seen in the slope of the power spectrum, which reflects aperiodic (non-oscillatory) activity. Measuring aperiodic activity is important in order to detect band-limited oscillations, and as a feature of interest that relates to underlying spiking levels, and to the balance of excitation and inhibition. This tutorial will cover spectral parameterization (specparam, formerly fooof), a method for measuring aperiodic activity and neural oscillations.

3. Synchrony: Measures of coupling between EEG recording sites often invoke the notion of synchrony as coherence or phase locking. These measures are sometimes used, in turn, to infer graph edges in network analyses. This tutorial will define and compare pairwise coherence, phase locking indices, and global coherence, as well as demonstrate some of their uses in the inference of functional connectivity.

4. Complexity: Over the past years, complexity-related measures have been shown to be reliable biomarkers of consciousness across several conditions. Among them, the perturbational complexity index has shown unprecedented sensitivity and specificity in detecting capacity for consciousness in healthy controls and patients with disorders of consciousness. This tutorial will cover measures of complexity (Perturbational Complexity Index, Lempel-Ziv complexity).

The interactive portion of the tutorial will be followed by a panel discussion and Q&A session regarding analytical approaches and their applications to consciousness science.

## Repository Layout

This tutorial has two components, slides and hands-on code examples.

If you want to try out the hands-on materials you will need a working Python installation (see below) and MATLAB.

The hands-on materials are organized as follows:

- `data:` contains some example data files that can be used to explore the methods
- `connectivity:` contains materials related to applying functional connectivity measures to neural data
- `freq_domain:` contains materials related to using frequency-domain spectral decomposition
- `time_domain:` contains materials using state space models to capture time-domain rhythms and their interactions.

## Software Environment

Part of the materials in this tutorial require the Python programming language, as well as Jupyter notebooks. Participants will need a working install of Python (>=3.6), with standard scientific packages (numpy, scipy, matplotlib), as well as some neuroscience specific modules.

To assist with getting an environment setup, we have provided an environment file which can be initialized as:

    > git clone ***
    > cd ASSC26_EEGMethods
    > conda env create -f environment.yml
    > conda activate ASSC26
    > jupyter lab

The above assumes you have the [conda](https://docs.conda.io/en/latest/) package manager, which can be downloaded as part of the [Anaconda distribution](https://www.anaconda.com/download).

Some of the tutorial requires a version of MATLAB compatible with .mlx files.

## Additional References

We will also cover topics related to these other repositories:

### PCIst
[https://www.github.com/renzocom/PCIst](https://www.github.com/renzocom/PCIst)

Comolatti R, Pigorini A, Casarotto S, Fecchio M, Faria G, Sarasso S, Rosanova M, Gosseries O, Boly M, Bodart O, Ledoux D, Brichant JF, Nobili L, Laureys S, Tononi G, Massimini M, Casali AG. A fast and general method to empirically estimate the complexity of brain responses to transcranial and intracranial stimulations. Brain Stimul. 2019 Sep-Oct;12(5):1280-1289. doi: 10.1016/j.brs.2019.05.013. Epub 2019 May 15. PMID: 31133480.

### TAAC tool
[https://www.github.com/iTCf/TAAC](https://www.github.com/iTCf/TAAC)

Russo S, Sarasso S, Puglisi GE, Dal Palù D, Pigorini A, Casarotto S, D'Ambrosio S, Astolfi A, Massimini M, Rosanova M, Fecchio M. TAAC - TMS Adaptable Auditory Control: A universal tool to mask TMS clicks. J Neurosci Methods. 2022 Mar 15;370:109491. doi: 10.1016/j.jneumeth.2022.109491. Epub 2022 Jan 31. PMID: 35101524.

### rt-TEP
[https://www.github.com/iTCf/rt-TEP.git](https://www.github.com/iTCf/rt-TEP.git)

Casarotto S, Fecchio M, Rosanova M, Varone G, D'Ambrosio S, Sarasso S, Pigorini A, Russo S, Comanducci A, Ilmoniemi RJ, Massimini M. The rt-TEP tool: real-time visualization of TMS-Evoked Potentials to maximize cortical activation and minimize artifacts. J Neurosci Methods. 2022 Mar 15;370:109486. doi: 10.1016/j.jneumeth.2022.109486. Epub 2022 Jan 21. PMID: 35074394.

## Citing

If you use any of the following code or methods contained in this repository in your work, please cite the following papers:

### Global coherence
Weiner VS, Zhou DW, Kahali P, Stephen EP, Peterfreund RA, Aglio LS, Szabo MD, Eskandar EN, Salazar-Gomez AF, Sampson AL, Cash SS. Propofol disrupts alpha dynamics in functionally distinct thalamocortical networks during loss of consciousness. Proceedings of the National Academy of Sciences. 2023 Mar 14;120(11):e2207831120.
Please reach out to dwzhou@mit.edu about citing the subspace methods.

## Contributors

This Tutorial was created by:
- [David W. Zhou](https://scholar.google.com/citations?user=aiiwYt0AAAAJ&hl=en), PhD, Carney Institute for Brain Science, Brown University, Providence, RI, USA
- [Emily P. Stephen](http://emilystephen.com/), PhD, Dept. of Mathematics and Statistics, Boston University, Boston, MA, USA
- [Thomas Donoghue](https://tomdonoghue.github.io/), PhD, Dept. of Biomedical Engineering, Columbia University, New York City, NY, USA
- Matteo Fecchio, PhD, Center for Neurotechnology and Neurorecovery, Massachusetts General Hospital, Boston, MA, USA

## Acknowledgements

Some of these materials are adapted from the short course "Advanced topics in the analysis of neural electrophysiology data: decomposing rhythmic and broadband components", also available
[here](https://github.com/TomDonoghue/NESSTutorial2023).
