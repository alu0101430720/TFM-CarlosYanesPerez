# =============================================================================
# S-AES - Simplified AES (16-bit block, 16-bit key)
# =============================================================================

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

SBOX_INT = [0x9, 0x4, 0xA, 0xB, 0xD, 0x1, 0x8, 0x5,
            0x6, 0x2, 0x0, 0x3, 0xC, 0xE, 0xF, 0x7]
SBOX     = {i: SBOX_INT[i] for i in range(16)}
SBOX_INV = {v: k for k, v in SBOX.items()}

RCON = {1: 0x80, 2: 0x30}

# ---------------------------------------------------------------------------
# Aritmética en GF(2^4)
# ---------------------------------------------------------------------------

def gf_mult(a, b, mod=0x13):
    """Multiplicación en GF(2^4)."""
    result = 0
    for _ in range(4):
        if b & 1:
            result ^= a
        hi = a & 0x8
        a  = (a << 1) & 0xF
        if hi:
            a ^= (mod & 0xF)
        b >>= 1
    return result

# ---------------------------------------------------------------------------
# Expansión de clave  ->  3 subclaves de 16 bits
# ---------------------------------------------------------------------------

def key_expansion(key):
    """Genera las tres subclaves k0, k1, k2 a partir de la clave maestra."""
    def sub_word(w): return (SBOX[(w >> 4) & 0xF] << 4) | SBOX[w & 0xF]
    def rot_word(w): return ((w & 0xF) << 4) | ((w >> 4) & 0xF)

    w    = [0] * 6
    w[0] = (key >> 8) & 0xFF
    w[1] =  key       & 0xFF
    w[2] = w[0] ^ RCON[1] ^ sub_word(rot_word(w[1]))
    w[3] = w[2] ^ w[1]
    w[4] = w[2] ^ RCON[2] ^ sub_word(rot_word(w[3]))
    w[5] = w[4] ^ w[3]

    return (w[0] << 8) | w[1], (w[2] << 8) | w[3], (w[4] << 8) | w[5]

# ---------------------------------------------------------------------------
# Cifrado S-AES
# ---------------------------------------------------------------------------

def saes_encrypt(P, K):
    """Cifra P (16 bits) con clave K (16 bits). Devuelve el criptograma C."""
    k0, k1, k2 = key_expansion(K)

    def nibs(x):       return [(x >> 12) & 0xF, (x >> 8) & 0xF,
                                (x >>  4) & 0xF,  x       & 0xF]
    def from_nibs(n):  return (n[0] << 12) | (n[1] << 8) | (n[2] << 4) | n[3]
    def ark(s, k):     kn = nibs(k); return [s[i] ^ kn[i] for i in range(4)]
    def sub(s):        return [SBOX[x] for x in s]
    def shr(s):        return [s[0], s[3], s[2], s[1]]   # intercambia n1 y n3
    def mc(s):
        return [
            gf_mult(1, s[0]) ^ gf_mult(4, s[1]),   # n0' = n0 + 4·n1
            gf_mult(4, s[0]) ^ gf_mult(1, s[1]),   # n1' = 4·n0 + n1
            gf_mult(1, s[2]) ^ gf_mult(4, s[3]),   # n2' = n2 + 4·n3
            gf_mult(4, s[2]) ^ gf_mult(1, s[3]),   # n3' = 4·n2 + n3
        ]

    s = nibs(P)
    s = ark(s, k0)                        # Ronda inicial
    s = sub(s); s = shr(s); s = mc(s)    # Ronda 1
    s = ark(s, k1)
    s = sub(s); s = shr(s)               # Ronda 2  (sin MixColumns)
    s = ark(s, k2)
    return from_nibs(s)


def saes_decrypt(C, K):
    """Descifra C (16 bits) con clave K (16 bits). Devuelve el texto P."""
    k0, k1, k2 = key_expansion(K)

    def nibs(x):       return [(x >> 12) & 0xF, (x >> 8) & 0xF,
                                (x >>  4) & 0xF,  x       & 0xF]
    def from_nibs(n):  return (n[0] << 12) | (n[1] << 8) | (n[2] << 4) | n[3]
    def ark(s, k):     kn = nibs(k); return [s[i] ^ kn[i] for i in range(4)]
    def sub_inv(s):    return [SBOX_INV[x] for x in s]
    def shr(s):        return [s[0], s[3], s[2], s[1]]   # misma operación (swap)
    def mc_inv(s):
        # Matriz inversa en GF(2^4): [[9,2],[2,9]]
        return [
            gf_mult(9, s[0]) ^ gf_mult(2, s[1]),
            gf_mult(2, s[0]) ^ gf_mult(9, s[1]),
            gf_mult(9, s[2]) ^ gf_mult(2, s[3]),
            gf_mult(2, s[2]) ^ gf_mult(9, s[3]),
        ]

    s = nibs(C)
    s = ark(s, k2)                        # Deshace ronda 2
    s = shr(s); s = sub_inv(s)
    s = ark(s, k1)                        # Deshace ronda 1
    s = mc_inv(s); s = shr(s); s = sub_inv(s)
    s = ark(s, k0)                        # Deshace ronda inicial
    return from_nibs(s)

# ---------------------------------------------------------------------------
# Vectorización NumPy  -  precomputa encrypt(P, k) para las 2^16 claves
#
# Equivalente exacto a saes_encrypt escalar, pero opera sobre todos los
# valores de k simultáneamente mediante arrays NumPy.
# ---------------------------------------------------------------------------

import numpy as np

# Lookup tables derivadas de las constantes escalares (fuente única de verdad)
_SBOX_V = np.array(SBOX_INT, dtype=np.uint16)

def _gf_mult_const(c, mod=0x13):
    """Precomputa la tabla de multiplicación GF(2^4) por la constante c."""
    table = np.zeros(16, dtype=np.uint16)
    for i in range(16):
        table[i] = gf_mult(c, i, mod)
    return table

_GF_MUL4 = _gf_mult_const(4)   # tabla para ×4
_GF_MUL9 = _gf_mult_const(9)   # tabla para ×9 (descifrado)
_GF_MUL2 = _gf_mult_const(2)   # tabla para ×2 (descifrado)


def _key_expansion_v(K: np.ndarray):
    """
    Expansión de clave vectorizada.
    K : array uint32 con shape (N,), valores en [0, 0xFFFF].
    Devuelve (k0, k1, k2), cada uno array uint32 shape (N,).
    """
    W0 = (K >> 8) & 0xFF
    W1 =  K       & 0xFF

    def sub_word(w):
        return (_SBOX_V[(w >> 4) & 0xF].astype(np.uint32) << 4) | _SBOX_V[w & 0xF]

    def rot_word(w):
        return ((w & 0xF) << 4) | ((w >> 4) & 0xF)

    W2 = (W0 ^ RCON[1] ^ sub_word(rot_word(W1))) & 0xFF
    W3 = (W2 ^ W1) & 0xFF
    W4 = (W2 ^ RCON[2] ^ sub_word(rot_word(W3))) & 0xFF
    W5 = (W4 ^ W3) & 0xFF

    return (W0 << 8) | W1, (W2 << 8) | W3, (W4 << 8) | W5


def build_ct_table(plaintext: int) -> np.ndarray:
    """
    Precomputa saes_encrypt(plaintext, k) para las 2^16 claves posibles.

    Parámetros
    ----------
    plaintext : int  - texto en claro de 16 bits.

    Devuelve
    --------
    table : np.ndarray uint16, shape (65536,)
        table[k] == saes_encrypt(plaintext, k)  para k en [0, 65535].
    """
    K = np.arange(65536, dtype=np.uint32)
    k0, k1, k2 = _key_expansion_v(K)

    # Descomponer subclaves en nibbles
    def nibs(x):
        return (x >> 12) & 0xF, (x >> 8) & 0xF, (x >> 4) & 0xF, x & 0xF

    n0_k0, n1_k0, n2_k0, n3_k0 = nibs(k0)
    n0_k1, n1_k1, n2_k1, n3_k1 = nibs(k1)
    n0_k2, n1_k2, n2_k2, n3_k2 = nibs(k2)

    # Estado inicial
    s0 = (plaintext >> 12) & 0xF
    s1 = (plaintext >>  8) & 0xF
    s2 = (plaintext >>  4) & 0xF
    s3 =  plaintext        & 0xF

    # Ronda inicial: AddRoundKey con k0
    s0 = s0 ^ n0_k0;  s1 = s1 ^ n1_k0
    s2 = s2 ^ n2_k0;  s3 = s3 ^ n3_k0

    # Ronda 1: SubNibbles -> ShiftRows -> MixColumns -> AddRoundKey(k1)
    s0 = _SBOX_V[s0]; s1 = _SBOX_V[s1]
    s2 = _SBOX_V[s2]; s3 = _SBOX_V[s3]
    s1, s3 = s3, s1                                          # ShiftRows
    ns0 = s0 ^ _GF_MUL4[s1];  ns1 = _GF_MUL4[s0] ^ s1     # MixColumns
    ns2 = s2 ^ _GF_MUL4[s3];  ns3 = _GF_MUL4[s2] ^ s3
    s0, s1, s2, s3 = ns0, ns1, ns2, ns3
    s0 = s0 ^ n0_k1;  s1 = s1 ^ n1_k1
    s2 = s2 ^ n2_k1;  s3 = s3 ^ n3_k1

    # Ronda 2: SubNibbles -> ShiftRows -> AddRoundKey(k2)
    s0 = _SBOX_V[s0]; s1 = _SBOX_V[s1]
    s2 = _SBOX_V[s2]; s3 = _SBOX_V[s3]
    s1, s3 = s3, s1                                          # ShiftRows
    s0 = s0 ^ n0_k2;  s1 = s1 ^ n1_k2
    s2 = s2 ^ n2_k2;  s3 = s3 ^ n3_k2

    table = (
        (s0.astype(np.uint16) << 12) |
        (s1.astype(np.uint16) <<  8) |
        (s2.astype(np.uint16) <<  4) |
         s3.astype(np.uint16)
    )
    return table


def verify_ct_table(plaintext: int, n_samples: int = 64) -> bool:
    """
    Comprueba que build_ct_table coincide con saes_encrypt escalar
    en n_samples claves aleatorias. Lanza AssertionError si divergen.
    """
    import random
    table = build_ct_table(plaintext)
    for _ in range(n_samples):
        k = random.randint(0, 0xFFFF)
        ref = saes_encrypt(plaintext, k)
        assert table[k] == ref, (
            f"Divergencia en P={plaintext:#06x}, K={k:#06x}: "
            f"escalar={ref:#06x}, vectorizado={table[k]:#06x}"
        )
    return True


# ---------------------------------------------------------------------------
# Demo / tests
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("S-Box S-AES:")
    for row in range(4):
        print("  " + "  ".join(
            f"{row*4+col:04b}->{SBOX_INT[row*4+col]:04b}"
            for col in range(4)
        ))

    print(f"\nVerificación GF(2^4): 4x4 = {gf_mult(4,4):#x}  (esperado: 0x3)")

    K = 0x544A
    k0, k1, k2 = key_expansion(K)
    print(f"\nExpansión de clave para K={K:#06x}:")
    print(f"  k0={k0:#06x}  k1={k1:#06x}  k2={k2:#06x}")

    test_vectors = [
        (0xDF64, 0x544A),
        (0x0000, 0x0000),
        (0xFFFF, 0xFFFF),
        (0x1234, 0xABCD),
    ]
    print("\nTest cifrado/descifrado:")
    for P, K in test_vectors:
        C = saes_encrypt(P, K)
        P2 = saes_decrypt(C, K)
        ok = "(Ok)" if P2 == P else "(Fallo)"
        print(f"  {ok}  P={P:#06x}  K={K:#06x}  ->  C={C:#06x}  ->  P'={P2:#06x}")

    print("\nEjemplo build_ct_table:")
    P = 0xDF64
    table = build_ct_table(P)
    print(f"  table[0x544A] = {table[0x544A]:#06x}  "
          f"(esperado: {saes_encrypt(P, 0x544A):#06x})")