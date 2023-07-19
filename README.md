# bamoso
Bayesian modelling of animal social interactions

## Installation

You need [`cmdstanr`](https://mc-stan.org/cmdstanr/) in order to install and run `bamoso`. 
This in turn requires a working C++ toolchain first. 
Check out the [getting started-guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) from `cmdstanr` to see how to get everything set up properly. 
Also [this document](https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html#cpp-toolchain) might be helpful.

Then install `cmdstanr`.

```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

And you also need the `remotes` packages, which is easy to install from CRAN:

```
install.packages("remotes")
```

Then check whether things are set up correctly:

```
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
```

You may need to follow any instructions that the output of that call gives.
If you haven't used `cmdstanr` before, it's likely that you at least need to install `cmdstan` itself:

```
cmdstanr::install_cmdstan()
```

If this gives positive feedback, install `bamoso`:

```
library(remotes)
remotes::install_github("gobbios/bamoso", dependencies = TRUE, build_vignettes = FALSE)
```

If you want to install (recompile) the intro vignette, use:

```
install.packages("tinytex")
library(remotes)
remotes::install_github("gobbios/bamoso", dependencies = TRUE, build_vignettes = TRUE)
```

This will take a little more time to finish, but if successful gives you access to the vignette:

```
vignette("intro", package = "bamoso")
```



