from __future__ import division, print_function
from collections import defaultdict
import os
import pandas as pd

_ortho2human = None
_ortho2mouse = None
_human2ortho = None
_mouse2ortho = None

this_directory = os.path.dirname(os.path.abspath(__file__))


def get_human_mouse_ortho():
    """
    Returns the ortho2human translation dictionary
    Translates between orthology group and human gene symbols
    """
    global _ortho2human, _ortho2mouse, _human2ortho, _mouse2ortho
    if _ortho2human is None:
        load_mgi()

    return _ortho2human, _ortho2mouse, _human2ortho, _mouse2ortho


def load_mgi():
    """
    Loads the ortho2human and ortho2mouse dictionaries
    from the MGI exported file
    """

    global _ortho2mouse, _ortho2human, _mouse2ortho, _human2ortho

    _ortho2mouse = defaultdict(set)
    _ortho2human = defaultdict(set)
    _mouse2ortho = defaultdict(set)
    _human2ortho = defaultdict(set)

    mgi_file = os.path.join(this_directory, "MGI",
                            "HOM_MouseHumanSequence.rpt")
    data = pd.read_table(mgi_file)

    hgs = data.loc[data['Common Organism Name'] == 'human']
    mgs = data.loc[data['Common Organism Name'] == 'mouse, laboratory']

    for ortho_id, hg in zip(hgs['HomoloGene ID'], hgs['Symbol']):
        _ortho2human[ortho_id].add(hg.upper())

    for ortho_id, mg in zip(mgs['HomoloGene ID'], mgs['Symbol']):
        _ortho2mouse[ortho_id].add(mg.upper())

    # Invert the dictionaries
    for ortho_id in _ortho2human:
        for hg in _ortho2human[ortho_id]:
            _human2ortho[hg].add(ortho_id)

    for ortho_id in _ortho2mouse:
        for hg in _ortho2mouse[ortho_id]:
            _mouse2ortho[hg].add(ortho_id)
