import numpy as np
import pandas as pd
from scipy.special import gamma

from .motifs import _ewms_to_dict
from .utils import reverse_compliment

def _scan_sequence_one_ewm(seq, ewm):
    """Internal function to scan one strand of a sequence with a single EWM."""
    motif_len = len(ewm)
    scores = np.zeros(len(seq) - motif_len + 1)
    for i in range(len(scores)):
        kmer = seq[i:i + motif_len]
        score = 0
        # This is faster than using the enumerate function
        for j in range(motif_len):
            score += ewm[j][kmer[j]]
        
        scores[i] = score

    return scores


def _scan_sequence(seq, ewms):
    """Internal function to scan both strands of a single sequence with a series of EWMs."""
    result = {}
    seq_rc = reverse_compliment(seq)
    for ewm_name, ewm in ewms.items():
        result[(ewm_name, "+")] = _scan_sequence_one_ewm(seq, ewm)
        result[(ewm_name, "-")] = _scan_sequence_one_ewm(seq_rc, ewm)[::-1]

    return result


def _scan_sequences(sequences, ewms):
    """Internal function to get the binding energy landscape for each EWM on each strand of each sequence.
    
    Parameters
    ----------
    sequences : pd.Series
        Series of sequences to scan.
    ewms : pd.Series
        Series of EWMs to use.

    Returns
    -------
    energy_landscapes : pd.Series
        Series of energy landscapes for each sequence. The values are dictionaries where the keys 
        are (motif, strand) and the values are the binding energy landscapes.
    """
    if type(sequences) is not pd.Series:
        raise ValueError("sequences must be a pandas Series")
    if sequences.name != "Sequence":
        raise ValueError("sequences does not appear to be a series of sequences. Make" +
                         "sure the name of the series is 'Sequence'.")
    if type(ewms) is not pd.Series:
        raise ValueError("ewms must be a pandas Series")
    if ewms.name != "EWM":
        raise ValueError("ewms does not appear to be a series of EWMs. Make" +
                         "sure the name of the series is 'EWM'.")
    
    ewms = _ewms_to_dict(ewms)
    energy_landscapes = sequences.apply(_scan_sequence, args=(ewms,))
    return energy_landscapes


def _get_mu(kd, rt):
    """Internal function to calculate the TF chemical potential from the minimum
    allowable relative Kd.
    
    mu corresponds to the value of ddG that results in a predicted occupancy of 0.5. This happens 
    when the two values are equal. Thus, mu = -RT ln(Kd).
    """
    return -rt * np.log(kd)


def _occupancy_landscape(energy_landscape, mu):
    """Internal function to convert an energy landscape to an occupancy landscape."""
    occupancy_landscape = {}
    for (motif, strand), landscape in energy_landscape.items():
        occupancy_landscape[(motif, strand)] = 1 / (1 + np.exp(landscape - mu))

    return occupancy_landscape


def _find_motifs_sequence(energy_landscape, kd, rt, motif_lengths):
    """Internal function to find all instances of motifs in one sequence."""
    mu = _get_mu(kd, rt)
    occupancy_landscapes = _occupancy_landscape(energy_landscape, mu)
    motif_hits = []
    colnames = ["motif", "strand", "start", "stop", "rel_kd"]
    for (motif, strand), energy in energy_landscape.items():
        energy = pd.Series(energy)
        hits = energy[occupancy_landscapes[(motif, strand)] > 0.5]
        for start, energy in hits.items():
            stop = start + motif_lengths[motif]
            motif_kd = np.exp(-energy / rt)
            motif_hits.append([motif, strand, start, stop, motif_kd])
    
    motif_hits = pd.DataFrame(motif_hits, columns=colnames)
    return motif_hits


def find_motifs(sequences, ewms, kd=0.0273, rt=2.5):
    """Find all instances of motifs in all sequences.
    
    Parameters
    ----------
    sequences : pd.Series
        Series of sequences to scan.
    ewms : pd.Series
        Series of EWMs to use.
    kd : float
        The minimum allowable relative affinity for a site to be considered bound by a TF.
    rt : float
        The value of RT to use in the formula, in kJ/mol. The default value corresponds to the ideal 
        gas constant times 300 Kelvin.
    
    Returns
    -------
    motif_hits : pd.Series
        All instances of motifs in all sequences. The values are DataFrames where the columns are 
        "motif", "strand", "start", "stop", and "rel_kd".
    """
    energy_landscapes = _scan_sequences(sequences, ewms)
    motif_lengths = ewms.apply(len)
    motif_hits = energy_landscapes.apply(_find_motifs_sequence, args=(kd, rt, motif_lengths))
    return motif_hits


def count_motifs(motif_hits):
    """Count the number of times each motif is found in each sequence.
    
    Parameters
    ----------
    motif_hits : pd.Series
        All instances of motifs in all sequences. The values are DataFrames where the columns are 
        "motif", "strand", "start", "stop", and "rel_kd".
    
    Returns
    -------
    motif_counts : pd.DataFrame
        The number of times each motif is found in each sequence. Rows are sequences with same 
        index as sequences, columns represent different motif IDs.
    """
    motif_counts = motif_hits.apply(lambda x: x["motif"].value_counts())
    motif_counts.fillna(0, inplace=True)
    return motif_counts


def _calculate_occupancy_sequence(energy_landscape, mu, motifs):
    """Internal function to calculate the total occupancy of each motif on one sequence."""
    occupancy = _occupancy_landscape(energy_landscape, mu)
    occupancy_profile = {}
    for motif in motifs:
        occupancy_profile[motif] = occupancy[(motif, "+")].sum() + occupancy[(motif, "-")].sum()

    occupancy_profile = pd.Series(occupancy_profile)
    return occupancy_profile


def calculate_occupancy(sequences, ewms, kd=0.0273, rt=2.5):
    """Calculate the total predicted occupancy of each motif on each sequence.
    
    Parameters
    ----------
    sequences : pd.Series
        Series of sequences to scan.
    ewms : pd.Series
        Series of EWMs to use.
    kd : float
        The minimum allowable relative affinity for a site to be considered bound by a TF.
    rt : float
        The value of RT to use in the formula, in kJ/mol. The default value corresponds to the ideal 
        gas constant times 300 Kelvin.
    
    Returns
    -------
    occupancies : pd.DataFrame
        Total predicted occupancy of each motif on each sequence. Rows are sequences with same 
        index as sequences, columns represent different motif IDs.
    """
    energy_landscapes = _scan_sequences(sequences, ewms)
    mu = _get_mu(kd, rt)
    motif_ids = ewms.index.values
    occupancies = energy_landscapes.apply(_calculate_occupancy_sequence, args=(mu, motif_ids))
    return occupancies


def _information_content_sequence(occupancies):
    """Internal function to calculate the information content of a single sequence."""
    total_occ = occupancies.sum()
    diversity = (occupancies > 0.5).sum()
    # Since the occupancies are continuous values, we need to use the Gamma function to compute entropy. Gamma(n+1)=n!
    # W = N! / prod(N_i!)
    microstates = gamma(total_occ + 1) / (occupancies + 1).apply(gamma).product()
    # S = log W
    info_content = np.log2(microstates)

    result = pd.Series({
        "total_occupancy": total_occ,
        "diversity": diversity,
        "information_content": info_content
    })
    return result


def information_content(occupancies):
    """Given a list of total TF occupancies for many sequences, compute the following for each sequence:
    - Total occupancy of all TFs
    - Diversity (number of TFs that are occupied)
    - Information content of the sequence
    
    Parameters
    ----------
    occupancies : pd.DataFrame
        Total predicted occupancy of each motif on each sequence. Rows are sequences with same 
        index as sequences, columns represent different motif IDs.
    
    Returns
    -------
    result : pd.DataFrame
        The total occupancy, diversity, and information content of each sequence.
    """
    result = occupancies.apply(_information_content_sequence, axis=1)
    return result
