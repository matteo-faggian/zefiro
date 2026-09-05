"""Che cosa costa una miscela non uniforme, e perche' le soglie sono quelle.

Questi test coprono il ponte fra le due meta' del progetto: la CFD misura una
disuniformita', il modello L0 calcola una spinta, e fino a `l0/mescolamento.py`
i due non si parlavano.
"""
from __future__ import annotations

import pytest

from zefiro.l0.mescolamento import perdita_di_mescolamento
from zefiro.l0.mixture import MixtureModel
from zefiro.schemas import FuelSpec

P_C, P_AMB = 636097.5, 101325.0
T_ARIA, T_GPL = 293.0, 283.0
MDOT_ARIA, MDOT_GPL = 0.03235383716362278, 0.0018570516373200116


@pytest.fixture(scope="module")
def modello():
    pytest.importorskip("cantera")
    return MixtureModel.from_fuel(FuelSpec(
        composition={"C3H8": 1.0}, phase_at_injection="gas",
        thermo_source="gri30.yaml"))


def perdita(modello, U):
    return perdita_di_mescolamento(modello, U, P_C, P_AMB, T_ARIA, T_GPL,
                                   MDOT_ARIA, MDOT_GPL)


def test_a_mescolamento_perfetto_non_si_perde_niente(modello):
    p = perdita(modello, 0.0)
    assert p.rendimento == pytest.approx(1.0)
    assert p.frazione_di_portata == pytest.approx(1.0)


def test_la_chiusura_a_due_parcelle_riproduce_media_e_varianza(modello):
    """La sola cosa che il modello sceglie e' la FORMA della distribuzione: i
    due momenti devono tornare esatti, altrimenti non e' la distribuzione della
    U misurata ma di un'altra."""
    Y_med = MDOT_GPL / (MDOT_ARIA + MDOT_GPL)
    for U in (0.25, 0.75, 1.50):
        p = perdita(modello, U)
        f, Y_r = p.frazione_di_portata, p.Y_ricco
        media = f * Y_r + (1.0 - f) * 0.0
        varianza = f * (Y_r - media) ** 2 + (1.0 - f) * (0.0 - media) ** 2
        assert media == pytest.approx(Y_med, rel=1e-12)
        assert varianza ** 0.5 / media == pytest.approx(U, rel=1e-12)


def test_il_rendimento_scende_e_non_risale_mai(modello):
    us = [0.0, 0.10, 0.25, 0.50, 0.75, 1.00]
    eta = [perdita(modello, u).rendimento for u in us]
    assert all(b <= a + 1e-12 for a, b in zip(eta, eta[1:])), eta


def test_le_soglie_dell_analisi_cfd_hanno_un_prezzo_accettabile(modello):
    """LE SOGLIE DI `sezioni.py` ERANO UN GIUDIZIO: qui prendono un significato.

    U < 0.10 = "la correlazione regge" e U < 0.25 = "mescolamento parziale ma
    accettabile" erano numeri scelti a mano. Se costassero il 20 % di spinta
    sarebbero soglie sbagliate, per quanto ragionevoli sembrassero."""
    assert perdita(modello, 0.10).rendimento > 0.99
    assert perdita(modello, 0.25).rendimento > 0.97


def test_una_disuniformita_di_uno_e_gia_meta_motore(modello):
    """U = 1 vuol dire scarto quadratico medio pari alla media. Non e' "un po'
    disuniforme": e' meta' della portata ad aria quasi pura."""
    p = perdita(modello, 1.0)
    assert p.rendimento < 0.70
    assert p.frazione_di_portata == pytest.approx(0.5, abs=0.02)


def test_oltre_il_limite_il_modello_si_rifiuta_invece_di_estrapolare(modello):
    """A U abbastanza grande la parcella ricca diventerebbe combustibile puro.
    Li' la chiusura a due valori non descrive piu' niente, e continuare a
    stampare un numero sarebbe peggio che fermarsi."""
    with pytest.raises(ValueError, match="due parcelle"):
        perdita(modello, 5.0)
