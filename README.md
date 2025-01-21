# BEEML2
BEEML2 is a python package to scan DNA sequences for transcription factor binding sites based on biophysical parameters. The name is subject to change, but is currently based on the name used by [Zhao et al., 2009](https://doi.org/10.1371/journal.pcbi.1000590). This package is based on code that I developed in graduate school and used in [Friedman et al., 2021](https://doi.org/10.7554/eLife.67403) and [Friedman et al., 2025](https://doi.org/10.1016/j.cels.2024.12.004).

## Installation
Create a conda environment, clone this repository, and then pip install it. For example:
```sh
conda create --name beeml2 python=3.11 -y
conda activate beeml2
git clone https://github.com/rfriedman22/beeml2.git
cd beeml2
pip install .
```

## Quickstart
Load in some motifs from a MEME-formatted file and some sequences from a FASTA file:
```python
from beeml2.io import read_meme, read_fasta
pwms = read_meme("examples/motifs.meme")
sequences = read_fasta("examples/sequences.fasta")
```

First, we need to convert the motifs from position weight matrices (PWMs) to energy weight matrices (EWMs). If you don't do this, downstream functions won't work.
```python
from beeml2.motifs import pwms_to_ewms
ewms = pwms_to_ewms(pwms)
```

There are two main ways you can use EWMs to scan your sequence:
- To find individual motif instances
- To calculate total TF occupancy
In both cases, you need to specify the relative $K_D$ (i.e. binding affinity relative to the optimal binding site) you want to use as a cutoff. The default is `0.0273`, which corresponds to the $\mu = 9$ value used in publications.

If you want to look for individual motif instances, you can run the following:
```python
from beeml2.occupancy import find_motifs
motif_hits = find_motifs(sequences, ewms)
motif_hits["chr1-135793436-135793600_CPPE"]
```
The result is a pandas Series where the indices are sequence IDs and the values are DataFrames containing the result of the motif scan for each sequence. For example, the above line will show a table indicating that the sequence has a CRX motif on the plus strand on the (zero-based) interval \[78, 86\) with ~9\% relative affinity, and an NRL motif on the minus strand on the interval \[136, 147\) with ~6% relative affinity.

You can use the result of `find_motifs` to calculate the number of motifs in each sequence:
```python
from beeml2.occupancy import count_motifs
count_motifs(motif_hits)
```

You can also save the result of `find_motifs` to a flat text file:
```python
from beeml2.io import save_motif_hits
save_motif_hits(motif_hits, "motif_hits.tsv")
```

If you want to calculate the total occupancy of each TF on each sequence, you can run:
```python
from beeml2.occupancy import calculate_occupancy
occupancies = calculate_occupancy(sequences, ewms)
```
The difference with calculating the total occupancy is that it considers the probability of the TF binding at each position of the sequence and then sums it up to obtain the total number of predicted molecules bound to the sequence. In practice, this is similar to counting the number of motifs, but it allows for continuous values (e.g. a site can have a predicted occupancy of 0.75).

Once you have the predicted occupancies of each TF on each sequence, you can calculate metrics such as the information content of the sequence:
```python
from beeml2.occupancy import information_content
information_content(occupancies)
```

## Contributing
If you would like to contribute code, please create a branch and open a pull request. Things that I would like to have built out:
- Interface for R users.
- Plotting functions and interfacing with Logomaker.
- Faster implementations that don't rely on nested DataFrames (probably torch).
- Interfaces with JASPAR and other databases of PWMs.
- Better API so there aren't as many imports.

I use [pip-tools](https://github.com/jazzband/pip-tools) to manage package dependencies. You can add it to your environment with:
```sh
pip install pip-tools
```
To change any dependencies, update the package dependencies in `requirements.in` and run the following two commands. The first one generates a new, full list of package dependencies, and the second installs anything new.
```sh
pip-compile
pip-sync
```