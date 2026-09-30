# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`seirvax`: an R package with an SEIR model (hospitalisation, deaths, parallel
vaccinated stratum, age × region groups) written in odin2 and run with dust2,
plus wrappers to build, run, fit and plot it. See `README.md` for the user-facing
description and `docs/` for the model description and how to change it.

## Running R on this machine

R is not installed on the host; it runs in a local Docker image (`gunnar-r`,
built from `../docker/Dockerfile`). odin2, dust2 and monty are not in the image
— they live in the gitignored project-local library `.Rlib/`. Run R with:

```sh
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e R_LIBS=/work/.Rlib \
  -v "$PWD":/work -w /work gunnar-r:latest Rscript -e 'devtools::test()'
```

This is specific to this machine, so keep it out of the repo's own documentation:
`README.md` and `docs/` assume people have their own R setup.

Useful commands in that container (`Rscript -e '...'`):

- `odin2::odin_package(".")` — regenerate C++ after editing `inst/odin/*.R`
- `devtools::document()`, `devtools::test()`, `devtools::check(manual = FALSE)`
- `devtools::install(".", lib = "/work/.Rlib", upgrade = "never")`

`devtools::load_all()` compiles without optimisation, so the stochastic model and
particle MCMC run several times slower than in an installed build. Install before
timing anything.

## Conventions

- `inst/odin/seir_ode.R` and `inst/odin/seir_stoch.R` must keep identical
  parameters and states; the R layer treats them interchangeably.
- Generated files (`R/dust.R`, `R/cpp11.R`, `src/`, `inst/dust/`) are never
  edited by hand — regenerate them.
- Follow the checklist in `docs/modifying_the_model.md` for any model change, and
  add a test in `tests/testthat/`.
