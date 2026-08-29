"""Il tempo chimico e' un numero su cui poi si decide di NON ottimizzare:
va quindi verificato con piu' cura, non con meno."""
from __future__ import annotations

import math

import pytest

from zefiro.l0.chemistry import (
    EXTINCTION_FRACTION,
    BlowoutTable,
    psr_blowout,
)
from zefiro.l0.mixture import MixtureModel

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def model(request):
    from zefiro.schemas import FuelSpec
    return MixtureModel.from_fuel(
        FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                 thermo_source="gri30.yaml")
    )


def test_blowout_scala_come_atteso(model):
    """Due controlli di SEGNO, indipendenti dal valore assoluto.

    1. tau_blowout deve DIMINUIRE con la pressione: la velocita' di reazione
       di una combustione idrocarburo/aria e' di ordine globale ~2, quindi il
       tempo caratteristico va come p^-1 circa. Se salisse, il calcolo sarebbe
       sbagliato.
    2. tau_blowout deve avere un MINIMO vicino allo stechiometrico, perche' li'
       la temperatura di fiamma e' massima e la cinetica e' esponenziale in T.
    """
    lo = psr_blowout(model, 3.0e5, 1.0, 300.0, 290.0).tau_blowout
    hi = psr_blowout(model, 6.0e5, 1.0, 300.0, 290.0).tau_blowout
    assert hi < lo, f"tau non decresce con p: {lo:.3e} -> {hi:.3e}"

    magra = psr_blowout(model, 5.0e5, 0.70, 300.0, 290.0).tau_blowout
    stech = psr_blowout(model, 5.0e5, 1.00, 300.0, 290.0).tau_blowout
    ricca = psr_blowout(model, 5.0e5, 1.15, 300.0, 290.0).tau_blowout
    assert stech < magra and stech < ricca, (magra, stech, ricca)


def test_soglia_di_estinzione_non_conta(model):
    """Il salto fra ramo acceso e spento e' netto: cambiare la soglia da 0.5 a
    0.3 o 0.7 non deve spostare il risultato piu' della tolleranza di bisezione.
    Se lo spostasse, tau_blowout non sarebbe una bifurcazione ma un artefatto."""
    import zefiro.l0.chemistry as chem

    base = psr_blowout(model, 5.22e5, 1.0, 300.0, 290.0).tau_blowout
    for soglia in (0.3, 0.7):
        chem.EXTINCTION_FRACTION = soglia
        try:
            alt = psr_blowout(model, 5.22e5, 1.0, 300.0, 290.0).tau_blowout
        finally:
            chem.EXTINCTION_FRACTION = EXTINCTION_FRACTION
        assert abs(alt - base) / base < 0.10, (soglia, base, alt)


def test_tabella_interpola_entro_errore_misurato(model):
    """L'errore della tabella non e' assunto: e' misurato contro il calcolo
    esatto in punti che NON stanno sulla griglia."""
    import numpy as np

    p_grid = np.linspace(3.0e5, 6.0e5, 4)
    phi_grid = np.linspace(0.70, 1.15, 5)
    tab = BlowoutTable.build(model, p_grid, phi_grid, 300.0, 290.0)

    errori = []
    for p, phi in ((3.7e5, 0.81), (4.9e5, 0.97), (5.6e5, 1.08)):
        esatto = psr_blowout(model, p, phi, 300.0, 290.0).tau_blowout
        errori.append(abs(tab(p, phi) - esatto) / esatto)
    assert max(errori) < 0.15, f"errore di interpolazione {max(errori):.1%}"


def test_fuori_griglia_e_un_errore(model):
    import numpy as np

    tab = BlowoutTable.build(model, np.linspace(4.0e5, 5.0e5, 2),
                             np.linspace(0.9, 1.1, 2), 300.0, 290.0)
    with pytest.raises(ValueError, match="fuori dalla griglia"):
        tab(6.0e5, 1.0)
