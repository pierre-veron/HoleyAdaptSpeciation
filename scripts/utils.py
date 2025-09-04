import numpy as np
import json
from termcolor import colored
alphabet = [chr(value) for value in range(97, 123)]
    
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def round_array_by_term(x, prec):
    out = np.zeros_like(x)
    for i in range(len(x)):
        out[i] = np.around(x[i], prec[i])
    return out

def round_significant(x, sig=4):
    x = np.asarray(x)  # Assure que l'entrée est un array NumPy
    result = np.zeros_like(x)  # Crée un array du même format que x
    mask = np.isfinite(x) & (x != 0)  # Identifie les valeurs valides (non-NaN et non-nulles)
    
    # Calcul uniquement pour les valeurs valides
    result[mask] = round_array_by_term(x[mask], sig - np.floor(np.log10(np.abs(x[mask]))).astype(int) - 1)
    
    # Garde les valeurs nulles (0) et NaN telles qu'elles
    result[~mask] = x[~mask]
    return result

def savejson(fp, data):
    """Wrapper for saving dictionnaries in a JSON file. A colored message is prompted.

    Args:
        fp (str): file name, with .json extension
        data (dict): data to be saved.
    """
    with open(fp.format(dir), "w") as f:
        json.dump(data, f, indent = 4, cls = NpEncoder)
    print(colored("\nData saved:", "green"))
    print(colored("     {}".format(fp), "blue"))

def savetxt(fp, obj, full_message = False, round = False, signif = 6):
    """Wrapper for saving arrays to a text file and prompt a message.

    Args:
        fp (str): file name, with .txt extention.
        obj (np.ndarray or list): data to save
        full_message (bool, optional): if True, a longer message is prompted. 
            Defaults to False.
        round (bool, optional): if True, the data are rounded to a number of 
            significant digits. Defaults to False.
        signif (int, optional): if round is True, set the number of significant 
            digits.
    """
    if not(isinstance(obj, np.ndarray)):
        obj = np.array(obj)
    if np.issubdtype(obj.dtype, np.str_):
        with open(fp, "w") as f:
            f.writelines(obj)
    else:
        if round: 
            obj = round_significant(obj, signif)
            fmt = dict(fmt = "%." + str(signif) + "g")
        else:
            fmt = dict()
        np.savetxt(fp, obj, **fmt)
    if full_message:
        print(colored("\nData saved:", "green"))
    print(colored("     {}".format(fp), "blue"))