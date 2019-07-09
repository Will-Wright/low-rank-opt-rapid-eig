Description
-----------

**Summary**

This software package solves large-scale phase retrieval and blind deconvolution problems by framing them as low-rank semidefinite programs. 

This package is a modification of the original package [low-rank-opt](https://www.cs.ubc.ca/~mpf/pubs/low-rank-spectral-optimization-via-gauge-duality/) ([mirror](https://github.com/Will-Wright/low-rank-opt-orig)) by Michael Friedlander.  The method is based on the paper [Low-rank spectral optimization via gauge duality](https://arxiv.org/abs/1508.00315) by Michael Friedlander and Ives Macedo.  

Our software package uses the same main algorithm (`saga_sd.m`) and helper functions as the original `low-rank-opt` package.  Our modifications are summarized below and detailed in the **Contents** section.

**Contributions**

This fork optimizes the original `saga_sd` algorithm for handling noisy phase retrieval problems.  The original algorithm is very effective for noiseless models, but has difficulty handling noisy models.  The modifications focus on two goals:
   1. providing an optimal method for solving the eigenvalue subproblem in the main routine, and
   2. using reliable termination conditions to terminate when the algorithm stagnates.


Demo Tutorials
--------------

The following experiments demonstrate the effectiveness of the `saga_sd` algorithm for phase retrieval and the efficiency of the new `adaptive_eigs_params` algorithm (which we call `saga_sd rapid-eig`) for solving the sequence of eigenvalue problems in `saga_sd`.

### Comparison of `saga_sd` and `saga_sd rapid-eig`

* `experiments.figure.noisyimage_adaptive_eig_full_exp`: runs a series of experiments (with results cached in this repository) demonstrating the efficiency of `saga_sd rapid-eig` for a variety of phase retrieval problems. Generates the following results (along with many other plots).

This plot demonstrates that `saga_sd rapid-eig` (red dashed line) is more efficient than `saga_sd` (blue solid line) for a variety of random phase retrieval problems with varying size (left), oversampling rate (middle), and noise ratio (right).

<p align="center">
<img src="Numerics-ada_vs_orig_various_params.png">
</p>


This table shows `saga_sd rapid-eig` is nearly optimal in the sense of selecting the best parameter (number of requested eigenvalues `r` in the MATLAB function 
[eigs](https://www.mathworks.com/help/matlab/ref/eigs.html)) for each iterate in the sequence of eigenvalue problems.  Note that `m` refers to the block-size in `eigs`, `L` is the phase retrieval oversampling rate, and `epsilon_rel` is the noise ratio.  See chapters 6 and 7 of the dissertation [A rapid eigenvalue method for noisy phase retrieval](https://github.com/Will-Wright/dissertation-rapid_eigenvalue_method_for_noisy_phase_retrieval/blob/master/will_wright_dissertation.pdf) for details.


<p align="center">
<img src="saga_sd-vs-rapid-eig.png">
</p>


Note: you may also run `experiments.figure.noisyimage_adaptive_eig_demo` for a smaller set of experiments and results.



### Comparison of `saga_sd` and [`wflow`](https://arxiv.org/abs/1407.1065) for noisy phase retrieval

* `experiments.table.noisyrandom_saga_vs_wflow_rel_errs`: demonstrates `saga_sd` is generally more likely to recover an accurate signal than the competitor `wflow` algorithm for noisy phase retrieval problems.

<p align="center">
<img src="saga_sd-vs-wflow.png">
</p>


* `experiments.figure.noisyimage_signal_relerr_various_saga_iterates`: solves two noisy image recovery problems (at 10% and 30% noise to signal) and displays figures from various iterates, demonstrating that this optimization method makes significant progress in 20 iterations. 
 
![denoised_images](/cache/noisy_image_parrot_iterates_and_x_rel_err.png)



Contents
--------

All contents of this package are functionally identical to the original `low-rank-opt` unless otherwise noted as **new** or **modified**.

* `./setpath.m`: script to add dependencies to the path. *Run this first.*

* `./saga_sd.m`: dual-descent solver
   * **Mod**: adaptive eigenvalue strategy for rapid computation of the eigenvalue subproblem
   * **Mod**: new termination conditions for noisy phase retrieval to identify stagnation
   * **Mod**: decreasing steplength strategy to replace linesearch when objective function is nondifferentiable
   * **Mod**: other minor alterations which do not affect internal behavior (see solver internal comments for details)

* `./project.m`: projection onto the dual gauge constraint
   * **Mod**: projection onto constraint set for prototyping l1-PhaseLift method

* `./pfd.m`: primal-from-dual recovery procedure

* `./dfp.m`: dual-from-primal recovery procedure ([optionally] used as spacer steps in `saga_sd`)

* `./spg.m`: projected-gradient method with Barzilai-Borwein steplengths and a non-monotonic linesearch

* `./+basis`: two implementations of restricted bases used for blind deconvolution (Haar DWT and a subset of the canonical basis)

* `./+cdp`: functions to generate Coded Diffraction Patterns for phase-recovery experiments

* `./+experiments`: functions to reproduce experiments
   * **Mod**: `/+exputil`, `/+figure`, `/+table`, and `/+test` subdirectories to demonstrate behavior of `saga_sd.m` for noisy phase retrieval

* `./+hop`: linear operators for phase recovery (`pl.m`) and blind deconvolution (`bd.m`)

* `./+kernel`: motion blur kernels used for blind deconvolution. One to match the kernel used by [Ahmed, Recht and Romberg](http://dx.doi.org/10.1109/TIT.2013.2294644) (`motion.m`) and two other variants to play with: `cmotion.m` and `dcmotion.m`

* `./+rwt`: the [Rice Wavelet Toolbox](https://github.com/ricedsp/rwt)

* `./+util`: utility functions

* `./data`: images used in experiments

* `./cache`: outputs from the experiments used to generate the tables in the paper

* `./external`: other third-party codes used and (possibly) modified to run our experiments and comparisons. Essentially, codes made available by [Mahdi Soltanolkotabi](http://www-bcf.usc.edu/~soltanol/WFcode.html), [Ali Ahmed](aliahmed.org/code.html), a [TFOCS](http://cvxr.com/tfocs/)-based solver for phase recovery very graciously provided by [Thomas Strohmer](https://www.math.ucdavis.edu/~strohmer/), and [Mark Schmidt's minFunc](http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
   * **Mod**: [Mark Schmidt's minConf](https://www.cs.ubc.ca/~schmidtm/Software/minConf.html) for comparison with `saga_sd.m`
   * **Mod**: `./external/eigensolvers` for adaptive eigenvalue method

Generating problem instances and solver inputs
----------------------------------------------

## Example: Phase-recovery experiments

Create options for the problem-generation routine `experiments.gendatapl`:
* 1-dim signal size 128
* random complex Gaussian signal
* octanary masks

```
gOpts = struct('n',128,'m',1,'signal','gaussian','mask','octanary');
```

Apply these solvers to the phase-recovery problem:
* `SAGA` (the gauge-dual approach)
* `SAGA-feas` (feasibility variant of `saga`)
* `TFOCS`
* `WFLOW`

```
setpath;
sOpts = struct('verbosity',0); # quiet, please!
stats(1) = experiments.pl('solver','saga','solverOpts',sOpts);
stats(2) = experiments.pl('solver','saga-feas','solverOpts',sOpts);
stats(3) = experiments.pl('solver','tfocs','solverOpts',sOpts);
stats(4) = experiments.pl('solver','wflow','solverOpts',sOpts);

# report relative errors in recovery
stats.xErrorRel
```

The script `experiments.pl` is a convenient way to generate data and apply a solver. For example:

1. generate data:
   * a random complex Gaussian signal `x0` of length `n = 128`
   * a measurement operator `A` based on octanary masks
   * the measurement vector `b`
   ```
   [A, b, x0, info] = ...
   experiments.gendatapl('n',128,'m',1,'signal','gaussian','mask','octanary');
   ```

2. Solve the PhaseLift problem for this instance using the gauge-dual approach:
  ```
  [x, r, info_sol] = saga_sd(A, b);
  ```

3. Report on the recovery error:
```
xError = util.hermitianerror(x(:),x0(:),'fro');
xErrorRel = xError / norm(x0)^2; % = ||xx'-x0x0'||_F/||x0x0'||_F
```

## Example: blind deconvolution
```
solver = 'saga';       % other options: `saga-feas` and `aug-lag`
filtername = 'motion'; % other options: `cmotion` and `dcmotion`
genOpts = struct('filtername',filtername);
stats = experiments.bd('solver',solver,'genOpts',genOpts);
```

Generating the tables in the paper
----------------------------------

First step is making sure you are at the main folder and then setting the paths
by calling:
```
setpath
```

The tables in the paper can be generated with the commands below.

*Warning:* Some of the tests can take significant amount of time because
they run many thousands of examples. We suggest that you first download
the solution cache file (2.4GB):
```
unzip('http://www.math.ucdavis.edu/~mpf/low-rank-opt/cache.zip')
```

* For the full random Gaussian noiseless experiments:
```
experiments.noiselesspl
```

* For the full random instances with noise:
```
experiments.noisypl
```

* For the large nebula image:
```
experiments.naturalpl
```

* For the blind deconvolution problem:
```
experiments.shapesbd
```

It might be interesting to run smaller problems to get acquainted with the
outputs and have a faster run. For that, the following are appropriate:

* For a subset of the random Gaussian signals (without noise):
```
experiments.noiselesspl('test',true)
```

* For a subset of the random instances with noise:
```
experiments.noisypl('test',true)
```

* For smaller (rescaled) versions of the nebula image (feel free to substitute
  the 1/50 parameter by another value for different rescalings):
```
experiments.naturalpl('resizeImage',1/50)
```

Data directory
--------------
* `nebula.tif`: resized image from [Hubble Telescope][].

[Hubble Telescope]: http://hubblesite.org/newscenter/archive/releases/nebula/2015/12/image/a/

Cache directory
---------------
This directory contains outputs from the experiments used to generate the tables
in the paper.

* `.mat` files are generated with the output info for each of the solvers
and the main parameters used. For example:
    - `solution_saga-feas_11_1.mat` contains the output info from running the
      solver `SAGA-feas` on the noiseless phase recovery problem with 11 masks
      and 2 as seed for the random number generator;
    - `solution_wflow_9_0_050_1.mat` contains the output info from running the
      solver `WFLOW` on the noisy phase recovery problem with 9 masks, 0.050 as
      eta and 1 as seed for the random number generator;
    - `wusterland_saga_10.mat` contains the output info from running the solver
      `SAGA` on the noiseless phase recovery problem for the nebula image with
      10 masks;
    - `shapesbd_aug-lag.mat` contains the output info from running the augmented
      Lagrangian solver on the blind deconvolution problem.

* `.png` files are generated to depict the images used in the blind
deconvolution problem. In this case, no matter what solver is used, the files
terminating on: `*_signal.png`, `*_kernel.png` and `*_measur.png` are all equal
across solvers as they depict the original input image, the original motion
blur kernel and the exact measurements, respectively. For example:
    - `aug-lag_signalEst.png` depicts the estimated (deblurred) image computed
      by the augmented Lagrangian solver;
    - `saga-feas_kernelEst.png` depicts the estimated motion blur kernel
      computed by the `SAGA-feas` solver;
    - `saga_measurEst.png` depicts the measurements computed with the estimated
      signal and motion blur kernel by the `SAGA` solver.

Miscellaneous
-------------

Create a zip archive of the repo with the command:
```
git archive --format=zip --prefix=low-rank-opt/ HEAD > /tmp/low-rank-opt.zip
```

Create a zip archive of the solution cache:
```
zip -r cache.zip cache
```
