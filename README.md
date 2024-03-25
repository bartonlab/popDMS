# popDMS

popDMS is a method developed by the [Barton lab](https://bartonlab.github.io) for inferring the functional effects of mutations from deep mutational scanning (DMS) experiments, also known as multiplexed assays for variant effects (MAVEs). Code for popDMS is written in Python3 and C++. For more information, see [this preprint](https://www.biorxiv.org/content/10.1101/2024.01.29.577759v1) describing popDMS and its application to simulations and a variety of DMS data sets.

### Software dependencies

Methods to infer epistasis are implemented in C++11 and make use of the [GNU Scientific Library](https://www.gnu.org/software/gsl/) and [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

Version 3.4.0 of Eigen that we use can be downloaded from this [link](https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip). For epistasis inference, this file should be unzipped into the `epistasis/` directory.


### Running popDMS

popDMS uses codon counts in [dms_tools format](http://jbloomlab.github.io/dms_tools/fileformats.html#deep-mutational-scanning-counts-file) or sequence counts in [MaveDB-HGVS format](https://www.mavedb.org/docs/mavehgvs/index.html) for input. For reference, [this link](https://github.com/bartonlab/paper-DMS-inference/blob/main/data/raw_data/FP16_DNA_codoncounts.csv) demonstrates the format for codon counts, and [this link](https://github.com/bartonlab/paper-DMS-inference/blob/main/data/raw_data/TpoR_nucleotide_counts.csv) shows an example file in MaveDB-HGVS format. 

Running popDMS differs slightly depending on the format of the input data.

#### Using codon counts

When using codon counts as input, we require three variables: `codon_counts_files`, `replicates`, and `times`. Here `codon_counts_files` contains a list of file paths to codon counts files. For each file, there must be a corresponding entry in the list `replicates` that identifies which replicate the file belongs to, and an entry in the list `times` that gives the time (in numbers of generations) that sequencing was performed to obtain this data. For examples, see [popDMS.ipynb](popDMS.ipynb).

#### Using sequence counts

When using sequence counts, we require five variables: `haplotype_counts_file`, `reference_sequence_file`, `n_replicates`, `time_points`, and `time_cols`. The variable `haplotype_counts_file` gives the path to a file containing the sequence counts. To normalize the selection coefficients relative to a reference sequence, `reference_sequence_file` should provide the path to a file storing the reference sequence in plain text (for an example, [see here](https://github.com/bartonlab/paper-DMS-inference/blob/main/data/raw_data/TpoR_reference_sequence.dat)). The total number of replicates is specified by `n_replicates`. For each replicate, the time(s) at which data was collected are given as a list in `time_points`. Finally, for each replicate, the variable `time_cols` points to the columns in the `haplotype_counts_file` that store sequence counts for each time. For examples, see [popDMS.ipynb](popDMS.ipynb).

#### Interpreting the output

For both approaches, popDMS will compute and save the variant frequencies needed to calculate selection coefficients. Using these files, the code to infer the selection coefficients can quickly be rerun using the `infer_independent` (for codon counts) or `infer_correlated` (for sequence counts) methods. Both methods will save a compressed comma separated values (CSV) file containing the inferred selection coefficients at the inferred optimal value regularization strength. The file can be unzipped to be viewed in plain text or with a program such as Microsoft Excel.

The columns of the selection coefficient are:
- `site`: Specifies the site at which the variant is observed, following the numbering of sites in the original input file
- `amino_acid`: Specifies the amino acid (including stops)
- `WT_indicator`: Set to `True` if the amino acid matches the reference at that site, and `False` otherwise
- `rep_x`: Values in these columns give the selection coefficients inferred for each replicate independently, with replicates numbered starting from 0 (`rep_0`)
- `joint`: Joint selection coefficients inferred across all replicates

This CSV file can be used for downstream analysis. We also provide a built-in plotting function `fig_dms` that takes a path to the CSV file as input and produces a heatmap of the inferred selection coefficients.


### Epistasis inference: (already merged in one bash file to run automatically)

The format of input data for epistasis inference is described in [this file](epistasis_inference/README_bash.txt). Once data has been stored in this format, inference of epistatic interactions proceeds by running the shell script `run_epistasis.sh` in the `epistasis` directory.


# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
