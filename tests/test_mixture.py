"""La cache della `ct.Solution`: e' un'ottimizzazione che tocca lo STATO
condiviso, quindi va verificata piu' di quanto costi.
"""
from __future__ import annotations

def test_la_cache_non_cambia_i_risultati(operating_point):
    """La `ct.Solution` e' condivisa e MUTABILE: se un pezzo di codice ne
    lasciasse lo stato sporco, la valutazione successiva partirebbe da li' e
    darebbe un numero diverso. Il difetto sarebbe silenzioso e dipendente
    dall'ordine, cioe' il peggior tipo di difetto possibile.

    Qui si valuta la stessa configurazione due volte, con una valutazione
    DIVERSA in mezzo, e si pretende identita' bit a bit."""
    import dataclasses

    from zefiro.geometry.parameters import default_design_vector, derive

    a = default_design_vector(p_c=4.0e5, phi_core=1.0, f_film=0.0, N_inj=12.0)
    b = default_design_vector(p_c=5.5e5, phi_core=0.8, f_film=0.2, N_inj=20.0)

    _, primo = derive(a, operating_point)
    _, _altro = derive(b, operating_point)
    _, secondo = derive(a, operating_point)

    for campo in ("T_ad", "c_star", "thrust", "Isp_s", "A_t", "q_throat",
                  "wall_volume", "gamma_c", "MW_c"):
        v1 = getattr(primo, campo)
        v2 = getattr(secondo, campo)
        assert v1 == v2, f"{campo}: {v1!r} != {v2!r} - stato sporco nella Solution"


def test_la_cache_restituisce_lo_stesso_oggetto():
    from zefiro.l0.mixture import solution

    assert solution("gri30.yaml") is solution("gri30.yaml")
