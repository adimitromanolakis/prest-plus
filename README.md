# PREST-plus

> **Note:** This repository is a copy of the official PREST-plus distribution hosted at
> https://utstat.toronto.edu/sun/Software/Prest/. Source version v4.09 (released 2010-10-05).
>
> The method is described in the paper:
> Sun L, Dimitromanolakis A.
> **PREST-plus identifies pedigree errors and cryptic relatedness in the GAW18 sample using genome-wide SNP data.**
> *BMC Proceedings* 2014;8(Suppl 1):S23. [doi:10.1186/1753-6561-8-s1-s23](https://doi.org/10.1186/1753-6561-8-s1-s23)
>
> Mirror written 2026-11-13.

PREST-plus is a software package for the detection of pedigree errors, cryptic relatedness and relationship misspecification in GWAS or linkage data. Using an optimized MLE estimator for IBD probabilities, PREST-plus computes accurate estimates for IBD0/1/2 using any number and combination of SNP or microsatellite marker data. It works as efficiently and accurately with a microsatellite linkage panel or SNP data from a GWAS study.

## Building

Requires `g++` (or another standard C++ compiler). The `Makefile` builds two binaries.

```sh
make          # builds `prest` and `convert-geno`
```

- `prest` — the main PREST-plus program
- `convert-geno` — utility to convert genotype input formats

Use `make clean` to remove build artifacts.

## Documentation

Full documentation is included in this repository as [PREST-plus v4.09.pdf](PREST-plus%20v4.09.pdf). It describes how to install and use PREST-plus on your dataset.

## License

This project is released open-source. See the [project home](https://utstat.toronto.edu/sun/Software/Prest/) for details.
