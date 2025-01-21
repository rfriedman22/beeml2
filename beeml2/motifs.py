import numpy as np
import pandas as pd

def _make_ewm(pwm, pseudocount, rt):
    """Internal function to convert a single PWM to an EWM."""
    pwm += pseudocount
    pwm = pwm.apply(lambda x: x / x.max(), axis=1)
    ewm = -rt * np.log(pwm)
    ewm.columns = ["A", "C", "G", "T"]
    return ewm


def _ewms_to_dict(ewms):
    """Convert a series of EWMs to a dictionary of matrices."""
    return ewms.apply(lambda x: x.to_dict(orient="index"))


def pwms_to_ewms(pwms, pseudocount=0.0001, rt=2.5):
    """Convert a series of position weight matrices (PWMs) to energy weight matrices (EWMs).

    An EWM is a matrix of relative free energies of TF binding for each base at each position in 
    the motif. The relative free energy is calculated as ddG = -RT ln(p_b,i / p_c,i), where p_b,i is
    the probability of base b, p_c,i is the probability of the consensus (highest probability) base, 
    and ddG is relative free energy.
    
    Normalizing the PWM to the maximum letter probability at each position expresses the relative
    binding affinity (Kd) of each base at that position.

    Parameters
    ----------
    pwms : pd.Series
        Series of PWMs, where each value is a pd.DataFrame with columns A, C, G, T.
    pseudocount : float
        Pseudocount value to add to every value to account for zeros in the PWM.
    rt : float
        The value of RT to use in the formula, in kJ/mol. The default value corresponds to the ideal 
        gas constant times 300 Kelvin.

    Returns
    -------
    ewms : pd.Series
        Series of EWMs, where each value is a pd.DataFrame with columns A, C, G, T.
    """
    if type(pwms) is not pd.Series:
        raise ValueError("pwms must be a pandas Series")
    if pwms.name != "PWM":
        raise ValueError("pwms does not appear to be a series of PWMs. Make" +
                         "sure the name of the series is 'PWM'.")
    
    ewms = pwms.apply(_make_ewm, args=(pseudocount, rt))
    ewms.name = "EWM"
    return ewms
