"""Il pezzo si stampa? Sbalzi, supporti impossibili, evacuazione della polvere.

Un canale interno non puo' avere supporti: nessuno puo' entrare a toglierli.
Quindi ogni superficie interna deve reggersi da sola, e la soglia non e'
un'opinione: sotto circa 45 gradi dal piano di costruzione, in SLM, la
superficie inferiore collassa nella polvere o esce con una rugosita' che in un
canale da 0.6 mm ne cambia la sezione idraulica.

Un canale con sezione QUADRATA ha il tetto orizzontale: 0 gradi. E' il caso
peggiore possibile, ed e' il motivo per cui i canali di raffreddamento veri
non sono mai quadrati.

Il secondo controllo e' la POLVERE. Un canale chiuso a valle e' una sacca:
la polvere non sinterizzata resta dentro, non esce con nessun lavaggio, e a
motore acceso viene fuori a pezzi. Serve un percorso continuo dal canale
all'esterno, e va verificato sulla topologia, non sperato.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

#: Angolo minimo fra una superficie rivolta VERSO IL BASSO e il piano di
#: costruzione perche' regga senza supporti. 45 gradi e' il valore su cui si
#: progetta di regola in SLM; il valore vero dipende da macchina, polvere e
#: parametri, ed e' parte del TODO 7. Qui e' esplicito e cambiabile.
SELF_SUPPORT_ANGLE_DEG = 45.0


def passo_elica_minimo(raggio: float, angolo_deg: float = SELF_SUPPORT_ANGLE_DEG) -> float:
    """Passo minimo di un canale elicoidale perche' il suo tetto si regga.

    Un canale elicoidale a raggio R e passo p ha tangente proporzionale a
    (1, 2 pi R / p) nel piano (asse, circonferenza). La normale della faccia
    inferiore ha componente assiale 2 pi R / p diviso la norma, quindi

        alpha = arctan( p / (2 pi R) )

    da cui **p >= 2 pi R tan(alpha)**. A 45 gradi si riduce a p >= 2 pi R: il
    canale deve avanzare lungo l'asse almeno quanto la circonferenza che
    percorre. Un canale quasi circonferenziale e' un anello, e un anello ha il
    tetto piatto.

    ATTENZIONE a una trappola: la prima stesura scriveva `2 pi R / tan(alpha)`.
    A 45 gradi tan vale 1 e le due formule COINCIDONO, quindi il caso su cui
    si verifica di solito non distingue fra la formula giusta e quella
    sbagliata. Fuori da 45 gradi divergono, e il segno dell'errore e' il
    peggiore: la formula sbagliata chiede MENO passo proprio quando si vuole
    piu' margine. Il test la verifica a 30 e a 60 gradi.

    E' il vincolo che ha bocciato il ramo della gola, che aveva p = 35 mm su
    un raggio fino a 13 mm.
    """
    return 2.0 * math.pi * raggio * math.tan(math.radians(angolo_deg))


@dataclass(frozen=True)
class RapportoStampa:
    direzione: tuple[float, float, float]
    n_triangoli: int
    area_totale: float
    area_verso_il_basso: float
    area_da_supportare: float
    area_sulla_piastra: float
    frazione_da_supportare: float
    angolo_minimo: float
    istogramma: dict[str, float]

    @property
    def ok(self) -> bool:
        return self.area_da_supportare <= 0.0


def analizza_sbalzi(vertici: np.ndarray, facce: np.ndarray,
                    direzione=(1.0, 0.0, 0.0),
                    soglia_deg: float = SELF_SUPPORT_ANGLE_DEG,
                    tolleranza_piastra: float = 2.0e-4) -> RapportoStampa:
    """Distribuzione degli angoli di sbalzo rispetto alla direzione di crescita.

    Convenzione: `alpha = arccos(|n . b|)` e' l'angolo fra la faccetta e il
    piano di costruzione. Una faccetta orizzontale ha alpha = 0 (peggiore), una
    verticale alpha = 90 (nessun problema). Serve supporto solo se la faccetta
    guarda in BASSO (`n . b < 0`) e alpha e' sotto soglia.

    L'area e' pesata: dieci triangoli minuscoli a 10 gradi non sono un problema,
    un centimetro quadrato di tetto piatto si'.
    """
    b = np.asarray(direzione, dtype=float)
    b = b / np.linalg.norm(b)
    a, c, e = vertici[facce[:, 0]], vertici[facce[:, 1]], vertici[facce[:, 2]]
    cross = np.cross(c - a, e - a)
    area = np.linalg.norm(cross, axis=1) / 2.0
    norma = np.linalg.norm(cross, axis=1, keepdims=True)
    n = cross / np.where(norma > 0, norma, 1.0)

    proiezione = n @ b
    alpha = np.degrees(np.arccos(np.clip(np.abs(proiezione), 0.0, 1.0)))
    verso_il_basso = proiezione < 0.0

    # La faccia che POGGIA sulla piastra di costruzione non e' uno sbalzo: e'
    # sostenuta dalla piastra stessa. Contarla gonfiava il verdetto di 6.5 cm2
    # su 12.5, cioe' oltre meta' del problema era immaginario.
    quota = vertici[facce].mean(axis=1) @ b
    sulla_piastra = verso_il_basso & (quota <= quota.min() + tolleranza_piastra)
    critico = verso_il_basso & (alpha < soglia_deg) & ~sulla_piastra

    bordi = [0, 10, 20, 30, 45, 60, 90.001]
    isto = {}
    for lo, hi in zip(bordi[:-1], bordi[1:]):
        m = verso_il_basso & ~sulla_piastra & (alpha >= lo) & (alpha < hi)
        isto[f"{lo:.0f}-{hi:.0f}"] = float(area[m].sum())

    return RapportoStampa(
        direzione=tuple(b.tolist()),
        n_triangoli=len(facce),
        area_totale=float(area.sum()),
        area_verso_il_basso=float(area[verso_il_basso].sum()),
        area_da_supportare=float(area[critico].sum()),
        area_sulla_piastra=float(area[sulla_piastra].sum()),
        frazione_da_supportare=float(area[critico].sum() / max(area.sum(), 1e-30)),
        angolo_minimo=(float(alpha[verso_il_basso & ~sulla_piastra].min())
                       if (verso_il_basso & ~sulla_piastra).any() else 90.0),
        istogramma=isto,
    )


def sbalzi_interni(vertici, facce, campo_gas, campo_canali,
                   direzione=(1.0, 0.0, 0.0),
                   soglia_deg: float = SELF_SUPPORT_ANGLE_DEG):
    """Come sopra, ma SOLO sulle superfici interne ai canali.

    E' la distinzione che conta: uno sbalzo sulla superficie esterna si
    supporta e poi si toglie; uno dentro un canale da 0.6 mm no, e resta li'.
    Si riconoscono le facce interne al circuito guardando dove cade il loro
    baricentro nel campo dei canali.
    """
    baricentri = vertici[facce].mean(axis=1)
    dentro = campo_canali.sample(baricentri) < 1.5 * campo_canali.grid.spacing
    if not dentro.any():
        return None
    return analizza_sbalzi(vertici, facce[dentro], direzione, soglia_deg)


def polvere_evacuabile(campo_canali, campo_solido, grid) -> tuple[bool, int, int]:
    """La polvere puo' uscire da tutta la rete di raffreddamento?

    Si etichettano le componenti connesse del VUOTO interno al pezzo e si
    guarda quante toccano il bordo del dominio, cioe' comunicano con l'esterno.
    Una componente che non lo tocca e' una sacca chiusa: quella polvere non
    esce con nessun lavaggio.

    Ritorna (tutto_evacuabile, n_componenti, n_sacche_chiuse).
    """
    from scipy import ndimage

    from zefiro.sdf.core import to_numpy

    # Si etichetta TUTTO cio' che non e' materiale - l'ambiente esterno e ogni
    # cavita' interna insieme - e si guarda in quale componente cade ciascun
    # voxel di canale. Se cade nella componente dell'esterno, comunica con
    # l'esterno e la polvere esce; se no, e' una sacca.
    #
    # La prima versione guardava se i canali toccavano il BORDO DELLA GRIGLIA,
    # che e' tutt'altra cosa: gli attacchi sbucano sulla superficie del pezzo,
    # non sul bordo del dominio di calcolo, e il test dava sempre "sacca".
    # Un test sbagliato che risponde "no" e' subdolo quanto uno che risponde
    # "si": in entrambi i casi non stai misurando quello che credi.
    non_materiale = to_numpy(campo_solido.a) >= 0
    etichette, n = ndimage.label(non_materiale)
    if n == 0:
        return True, 0, 0
    esterno = etichette[0, 0, 0]
    if esterno == 0:
        raise ValueError("l'angolo della griglia e' dentro il materiale: "
                         "griglia troppo stretta per distinguere l'esterno")
    canali = to_numpy(campo_canali.a) < 0
    componenti_canale = set(np.unique(etichette[canali]).tolist())
    componenti_canale.discard(0)
    sacche = sorted(componenti_canale - {esterno})
    return (not sacche), int(n), len(sacche)
