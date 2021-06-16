import os
import requests
# import urllib.parse

from flask import redirect, render_template, request, session
from functools import wraps


def apology(message):
    """Render message as an apology to user."""
    return render_template("apology.html", alert=message)


def mes5(seq):
    """ Call 5 prime end MaxEntScan algorithm http://genes.mit.edu/burgelab/maxent/download/"""
    # Seq example: CAGgtaagt

    if len(seq) != 9:
        raise Exception

    # Run maxentscan
    output = (os.popen('cd tools/maxentscan; perl score5_mod.pl %s' % (seq)).readline())

    # If the score is 0 it is outputted as an empty string
    if output == '':
        return 0.00
    else:
        return float(output)


def mes3(seq):
    """ Call 5 prime end MaxEntScan algorithm http://genes.mit.edu/burgelab/maxent/download/"""
    # Seq example: ttccaaacgaacttttgtagGGA

    if len(seq) != 23:
       raise Exception

    # Run maxentscan
    output = (os.popen('cd /tools/maxentscan; perl score3_mod.pl %s' % (seq)).readline())

    # If the score is 0 it is outputted as an empty string
    if output == '':
        return 0.00
    else:
        return float(output)

def fixable(sift, polyphen, wtme5, mutme5, wtme3, mutme3):
    """ Determines if this mutation is ammenable"""

    # If the amino acid isn't damaging and there is a splice prediction
    if ((sift == 'LOW' or float(sift) >= 0.05) or (polyphen < 0.75)) or ((float(wtme5) < float(mutme5)) or (float(wtme3) < float(mutme3))):
        return "There is a good chance that this mutation is ammenable!"
    # If the amino acid is damaging
    elif ((sift == 'HIGH' or float(sift) <= 0.05)) or (polyphen >= 0.8):
        return "Sorry this mutation is not ammenable by ASO therapy because the nucleotide change is damaging the protein."

    else:
        "This is a tricky one...you'll need to take a closer look."