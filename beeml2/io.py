import itertools

import numpy as np
import pandas as pd

def _fasta_iter(fin, sep=""):
    """A generator function to parse through one entry in a FASTA or FASTA-like file.

    Parameters
    ----------
    fin : file input stream
        Handle to the file to parse.
    sep : str
        Delimiter for adjacent bases in the file.

    Yields
    -------
    header : str
        Name of the sequence.
    sequence : str
        The sequence.
    """
    # Generator yields True if on a header line
    generator = itertools.groupby(fin, lambda x: len(x) > 0 and x[0] == ">")
    for _, header in generator:
        header = list(header)[0].strip()[1:]
        sequence = sep.join(i.strip() for i in generator.__next__()[1])
        yield header, sequence


def read_fasta(filename):
    """Parse through a FASTA file and store the sequences as a Series.

    Parameters
    ----------
    filename : str
        Name of the file.

    Returns
    -------
    sequences : pd.Series
        Index is the FASTA header, values are the sequence strings.
    """
    sequences = {}
    with open(filename) as fin:
        for header, sequence in _fasta_iter(fin):
            sequence = sequence.upper()
            sequences[header] = sequence
    
    sequences = pd.Series(sequences)
    sequences.name = "Sequence"
    return sequences


def read_meme(filename):
    """Read in a MEME file and extract motifs into a pandas Series. Each motif is a DataFrame.

    Adapted from TangerMEME https://github.com/jmschrei/tangermeme/blob/9fd7b2da3fe28b8651fbe9efb294d617c4fd4429/tangermeme/io.py#L392

    Parameters
    ----------
    filename : str
        The path to the MEME file to be read.

    Returns:
    ----------
    motifs : pd.Series
        A pandas Series where the index is the motif name and the value is a DataFrame containing 
        the position weight matrix (PWM) for each motif. The DataFrame columns are labeled 
        "A", "C", "G", and "T".
    """
    motifs = {}

    with open(filename) as fin:
        line = fin.readline().strip()
        if not line.startswith("MEME version"):
            raise ValueError("File does not appear to be a MEME file")
        
        motif, motif_len, i = None, None, 0
        for line in fin:
            if motif is None:
                if line.startswith("MOTIF"):
                    motif = line.replace("MOTIF", "").strip("\r\n")
            
            elif motif_len is None:
                if line.startswith("letter"):
                    motif_len = int(line.split()[5])
                    pwm = np.zeros((motif_len, 4))

            elif i < motif_len:
                pwm[i] = list(map(float, line.strip("\r\n").split()))
                i += 1

            else:
                motifs[motif] = pd.DataFrame(pwm, columns=["A", "C", "G", "T"])
                motif, motif_len, i = None, None, 0

        if motif is not None or motif_len is not None or i > 0:
            raise ValueError("Unexpected end of file")

    motifs = pd.Series(motifs)
    motifs.name = "PWM"
    return motifs


def save_motif_hits(motif_hits, filename):
    """Save the output of find_motifs() as a plain text file.
    
    Parameters
    ----------
    motif_hits : pd.Series
        All instances of motifs in all sequences. The values are DataFrames where the columns are 
        "motif", "strand", "start", "stop", and "rel_kd".
    filename : str
        The name of the file to save the output.
    """
    combined_df = pd.concat(motif_hits.values, keys=motif_hits.index).reset_index(level=0)
    combined_df.rename(columns={'level_0': 'sequence'}, inplace=True)
    combined_df.to_csv(filename, sep='\t', index=False)
