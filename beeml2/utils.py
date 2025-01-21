def reverse_compliment(seq):
    """Take the reverse compliment of a sequence.

    Parameters
    ----------
    seq : str
        The original sequence.

    Returns
    -------
    new_seq : str
        The reverse compliment.
    """
    compliment = {"A": "T", "C": "G", "G": "C", "T": "A"}
    new_seq = seq[::-1]
    new_seq = "".join([compliment[i] for i in new_seq])
    return new_seq