"""Filtri sulla densita': regolarizzazione, nitidezza, stampabilita'.

Tre filtri in cascata, ciascuno con un compito preciso:

1. **Filtro di densita'** — senza, l'ottimizzazione topologica produce la
   scacchiera: un'alternanza pieno/vuoto elemento per elemento che il FEM
   valuta artificialmente rigida. E' un artefatto numerico, non una struttura.
   Il filtro impone anche una scala di lunghezza minima, e con essa
   l'indipendenza dalla mesh: raffinando, il risultato converge invece di
   frammentarsi.

2. **Proiezione di Heaviside** — il filtro di densita' lascia bordi sfumati,
   e una densita' 0.5 non e' un materiale: non si stampa. La proiezione,
   con continuazione su beta, spinge le densita' verso 0 o 1 mantenendo la
   derivabilita'.

3. **Filtro di stampabilita' (overhang)** — un elemento pieno puo' esistere solo
   se sotto di lui c'e' materiale che lo sostenga entro l'angolo ammesso. Senza
   questo, l'ottimizzazione produce isole sospese e sbalzi che in SLM
   richiederebbero supporti dentro i condotti, dove non si tolgono. E' il filtro
   che rende il risultato un pezzo invece di un disegno.

Ogni filtro espone `apply` e `chain` (la trasposta jacobiana per la catena delle
derivate). La correttezza dell'intera catena e' verificata per differenze
finite in `tests/test_topopt.py`: se una derivata e' sbagliata, l'ottimizzatore
converge tranquillamente verso qualcosa che non e' un ottimo.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


# --------------------------------------------------------------------------- #
# 1. filtro di densita'
# --------------------------------------------------------------------------- #
def density_filter_matrix(nx: int, ny: int, rmin: float,
                          dx: float = 1.0, dy: float = 1.0) -> sp.csr_matrix:
    """Matrice H normalizzata del filtro conico di raggio `rmin`.

        x_filtrato = H x,   H_ef = max(0, rmin - dist(e,f)) / SOMMA_f (...)

    `rmin` e' in unita' FISICHE, non in elementi: cosi' raffinando la mesh la
    scala di lunghezza resta la stessa, che e' la condizione perche' il
    risultato sia indipendente dalla discretizzazione.
    """
    if rmin <= 0:
        raise ValueError("rmin deve essere positivo")
    ix = int(np.ceil(rmin / dx))
    iy = int(np.ceil(rmin / dy))
    rows, cols, vals = [], [], []
    for i in range(nx):
        for j in range(ny):
            e = i * ny + j
            for k in range(max(i - ix, 0), min(i + ix + 1, nx)):
                for m in range(max(j - iy, 0), min(j + iy + 1, ny)):
                    d = np.hypot((i - k) * dx, (j - m) * dy)
                    w = rmin - d
                    if w > 0.0:
                        rows.append(e); cols.append(k * ny + m); vals.append(w)
    H = sp.coo_matrix((vals, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    s = np.asarray(H.sum(axis=1)).ravel()
    return sp.diags(1.0 / s) @ H


# --------------------------------------------------------------------------- #
# 2. proiezione di Heaviside
# --------------------------------------------------------------------------- #
def heaviside(x: np.ndarray, beta: float, eta: float = 0.5) -> np.ndarray:
    """Proiezione regolarizzata verso 0/1.

        x_proj = [tanh(beta*eta) + tanh(beta*(x - eta))] /
                 [tanh(beta*eta) + tanh(beta*(1 - eta))]

    Per beta -> 0 e' l'identita', per beta -> inf e' un gradino in `eta`.
    Si parte da beta piccolo e lo si alza per gradi (continuazione): partire
    subito con beta alto rende il problema non convesso fin dall'inizio e
    l'ottimizzatore si blocca in un minimo locale pessimo.
    """
    if beta <= 1e-9:
        return x.copy()
    t = np.tanh(beta * eta)
    return (t + np.tanh(beta * (x - eta))) / (t + np.tanh(beta * (1.0 - eta)))


def heaviside_derivative(x: np.ndarray, beta: float, eta: float = 0.5) -> np.ndarray:
    if beta <= 1e-9:
        return np.ones_like(x)
    t = np.tanh(beta * eta)
    return beta * (1.0 - np.tanh(beta * (x - eta)) ** 2) / (t + np.tanh(beta * (1.0 - eta)))


# --------------------------------------------------------------------------- #
# 3. filtro di stampabilita'
# --------------------------------------------------------------------------- #
#: Parametri di regolarizzazione di min e max. Sono una SCELTA, non un dato:
#: piu' alti danno un filtro piu' vicino al min/max esatto ma peggio
#: condizionato. La verifica che conta non e' il valore di questi numeri, ma il
#: controllo GEOMETRICO sull'uscita (`max_overhang_violation`), che misura
#: direttamente se il risultato e' stampabile.
SMOOTH_MAX_P = 40.0
SMOOTH_MIN_Q = 60.0


def _smooth_max(v: np.ndarray) -> np.ndarray:
    """Massimo regolarizzato per norma-P su `v` (k, n), lungo l'asse 0.

        S = ( SOMMA_k v_k^P )^(1/P)

    Calcolato in forma di RAPPORTO rispetto al massimo esatto, cosi' non ci
    sono potenze di numeri piccolissimi che vanno in overflow: e' la stessa
    formula, scritta in modo che un calcolatore la sappia valutare.
    """
    M = np.max(v, axis=0)
    safe = np.where(M > 1e-12, M, 1.0)
    S = safe * np.power(np.sum(np.power(v / safe, SMOOTH_MAX_P), axis=0),
                        1.0 / SMOOTH_MAX_P)
    return np.where(M > 1e-12, S, 0.0)


def _smooth_max_grad(v: np.ndarray, S: np.ndarray) -> np.ndarray:
    """dS/dv_k = (v_k / S)^(P-1). Il rapporto sta in [0, 1]: sempre stabile."""
    safe = np.where(S > 1e-12, S, 1.0)
    g = np.power(np.clip(v / safe, 0.0, 1.0), SMOOTH_MAX_P - 1.0)
    return np.where(S > 1e-12, g, 0.0)


def _smooth_min(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Minimo regolarizzato: s = (a^-Q + b^-Q)^(-1/Q), in forma di rapporto.

        s = m * (1 + (m/M)^Q)^(-1/Q),  m = min(a,b),  M = max(a,b)

    Algebricamente identica, ma (m/M) sta in [0, 1] e non produce overflow.
    """
    m = np.minimum(a, b)
    M = np.maximum(a, b)
    safe = np.where(M > 1e-12, M, 1.0)
    ratio = np.power(np.clip(m / safe, 0.0, 1.0), SMOOTH_MIN_Q)
    return np.where(M > 1e-12, m * np.power(1.0 + ratio, -1.0 / SMOOTH_MIN_Q), 0.0)


def _smooth_min_grad(a: np.ndarray, b: np.ndarray, s: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray]:
    """ds/da = (s/a)^(Q+1), ds/db = (s/b)^(Q+1). Rapporti in [0, 1]."""
    da = np.power(np.clip(np.divide(s, np.where(a > 1e-12, a, 1.0)), 0.0, 1.0),
                  SMOOTH_MIN_Q + 1.0)
    db = np.power(np.clip(np.divide(s, np.where(b > 1e-12, b, 1.0)), 0.0, 1.0),
                  SMOOTH_MIN_Q + 1.0)
    return np.where(a > 1e-12, da, 0.0), np.where(b > 1e-12, db, 0.0)


def printable(x: np.ndarray, nx: int, ny: int) -> np.ndarray:
    """Densita' effettivamente stampabile, strato per strato dal basso.

    Direzione di costruzione: +y (j crescente). Un elemento puo' essere pieno
    solo fin dove lo sostiene il materiale sottostante, preso fra i tre elementi
    che lo toccano nello strato precedente. Tre elementi corrispondono a un
    angolo di sbalzo di 45 gradi su una griglia quadrata, che e' il limite
    tipico dell'SLM (quello vero della tua macchina e' il TODO K).

        xi_1 = x_1
        xi_i = min( x_i , max(xi_(i-1, j-1), xi_(i-1, j), xi_(i-1, j+1)) )
    """
    return _printable_forward(x, nx, ny)[0]


def _printable_forward(x: np.ndarray, nx: int, ny: int):
    """Passata in avanti, conservando cio' che serve all'adjoint."""
    X = x.reshape(nx, ny)
    Xi = np.empty_like(X)
    Xi[:, 0] = X[:, 0]
    cache = []
    for j in range(1, ny):
        below = Xi[:, j - 1]
        stack = np.vstack([np.concatenate(([0.0], below[:-1])),
                           below,
                           np.concatenate((below[1:], [0.0]))])
        S = _smooth_max(stack)
        xi = _smooth_min(X[:, j], S)
        cache.append((stack, S, X[:, j].copy(), xi.copy()))
        Xi[:, j] = xi
    return Xi.reshape(-1), cache


def printable_chain(x: np.ndarray, nx: int, ny: int,
                    dfdxi: np.ndarray) -> np.ndarray:
    """Trasposta jacobiana ESATTA del filtro, per passata all'indietro.

    Il filtro e' ricorsivo verso l'alto, quindi la sua adjoint scorre verso il
    basso accumulando: la sensibilita' di uno strato ricade sui tre elementi
    che lo sostenevano. Non e' un'approssimazione — e' la regola della catena
    applicata alla ricorsione, e i test la confrontano con differenze finite.
    """
    _, cache = _printable_forward(x, nx, ny)
    lam = dfdxi.reshape(nx, ny).astype(float).copy()
    g = np.zeros((nx, ny))
    for j in range(ny - 1, 0, -1):
        stack, S, xj, xi = cache[j - 1]
        da, db = _smooth_min_grad(xj, S, xi)
        g[:, j] = lam[:, j] * da
        up = lam[:, j] * db                       # verso il supporto S
        w = _smooth_max_grad(stack, S)            # (3, nx) pesi sui vicini
        contrib = np.zeros(nx)
        contrib[:-1] += (up * w[0])[1:]           # il vicino "sinistro" e' j-1 a i-1
        contrib += up * w[1]
        contrib[1:] += (up * w[2])[:-1]
        lam[:, j - 1] += contrib
    g[:, 0] = lam[:, 0]
    return g.reshape(-1)


def max_overhang_violation(xi: np.ndarray, nx: int, ny: int,
                           threshold: float = 0.5) -> float:
    """Controllo GEOMETRICO sull'uscita: quanta materia sta senza sostegno.

    Ritorna la frazione di elementi pieni (sopra `threshold`) che non hanno
    almeno un elemento pieno fra i tre sottostanti. Deve essere zero, o quasi,
    su un progetto stampabile. E' la verifica che conta davvero: misura la
    proprieta' voluta invece di fidarsi dei parametri del filtro.
    """
    X = (xi.reshape(nx, ny) > threshold)
    if not X.any():
        return 0.0
    bad = 0
    for j in range(1, ny):
        below = X[:, j - 1]
        left = np.concatenate(([False], below[:-1]))
        right = np.concatenate((below[1:], [False]))
        supported = left | below | right
        bad += int(np.sum(X[:, j] & ~supported))
    return bad / max(int(X.sum()), 1)
