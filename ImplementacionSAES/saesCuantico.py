
# Depende de algunas funciones del SAES

from saes import saes_encrypt, key_expansion, SBOX_INT, SBOX_INV
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
import time

def build_sbox_nibble():
    """
    S-Box cuántica para un nibble de 4 bits.
    """
    qin  = QuantumRegister(4, 'in')
    qout = QuantumRegister(4, 'out')
    qc   = QuantumCircuit(qin, qout, name='SBOX')
    for inp, out in enumerate(SBOX_INT):
        if out == 0: continue
        ctrl_state = format(inp, '04b')
        for bit in range(4):
            if (out >> bit) & 1:
                qc.mcx(list(qin), qout[bit], ctrl_state=ctrl_state)
    return qc

def build_sbox_inv_nibble():
    """S-Box inversa cuántica para un nibble."""
    qin  = QuantumRegister(4, 'in')
    qout = QuantumRegister(4, 'out')
    qc   = QuantumCircuit(qin, qout, name='SBOX_INV')
    for inp, out in SBOX_INV.items():
        if out == 0: continue
        ctrl_state = format(inp, '04b')
        for bit in range(4):
            if (out >> bit) & 1:
                qc.mcx(list(qin), qout[bit], ctrl_state=ctrl_state)
    return qc

def build_shift_rows():
    """ShiftRows: SWAP entre n1 (state[8..11]) y n3 (state[0..3])."""
    state = QuantumRegister(16, 'state')
    qc    = QuantumCircuit(state)
    for i in range(4):
        qc.swap(state[8+i], state[0+i])
    return qc

def apply_mult4(qc, src, tgt):
    """
    tgt ^= 4*src en GF(2^4) con x^4+x+1.
    Mapeo verificado: [a2, a3^a2, a0^a3, a1] (little-endian)
    """
    qc.cx(src[2], tgt[0])
    qc.cx(src[3], tgt[1]); qc.cx(src[2], tgt[1])
    qc.cx(src[0], tgt[2]); qc.cx(src[3], tgt[2])
    qc.cx(src[1], tgt[3])

def build_mix_columns():
    """
    MixColumns: opera sobre pares de nibbles en GF(2^4).
    Escribe el resultado en aux_mc y hace SWAP con state.
    Mapeo de columnas post-ShiftRows:
      col0: [state[12..15], state[8..11]] -> n0' en aux_mc[12..15], n1' en aux_mc[8..11]
      col1: [state[4..7],   state[0..3]]  -> n2' en aux_mc[4..7],   n3' en aux_mc[0..3]
    """
    state  = QuantumRegister(16, 'state')
    aux_mc = QuantumRegister(16, 'aux_mc')
    qc     = QuantumCircuit(state, aux_mc)

    # col0: n0'=na0+4*nb0, n1'=4*na0+nb0
    na0 = [state[12+i]  for i in range(4)]
    nb0 = [state[8+i]   for i in range(4)]   # post-SR: state[8..11] tiene n1
    oa0 = [aux_mc[12+i] for i in range(4)]   # n0'
    ob0 = [aux_mc[8+i]  for i in range(4)]   # n1'
    for i in range(4): qc.cx(na0[i], oa0[i])
    apply_mult4(qc, nb0, oa0)
    for i in range(4): qc.cx(nb0[i], ob0[i])
    apply_mult4(qc, na0, ob0)

    # col1: n2'=na1+4*nb1, n3'=4*na1+nb1
    na1 = [state[4+i]   for i in range(4)]
    nb1 = [state[0+i]   for i in range(4)]   # post-SR: state[0..3] tiene n3
    oa1 = [aux_mc[4+i]  for i in range(4)]   # n2'
    ob1 = [aux_mc[0+i]  for i in range(4)]   # n3'
    for i in range(4): qc.cx(na1[i], oa1[i])
    apply_mult4(qc, nb1, oa1)
    for i in range(4): qc.cx(nb1[i], ob1[i])
    apply_mult4(qc, na1, ob1)

    for i in range(16): qc.swap(aux_mc[i], state[i])
    return qc

def build_saes_circuit(P, K):
    """
    Circuito S-AES cuántico completo (48 qubits).
    Las subclaves se pasan clásicamente para reducir el número de qubits.
    """
    k0, k1, k2 = key_expansion(K)

    state   = QuantumRegister(16, 'state')
    aux_sub = QuantumRegister(16, 'aux_sub')
    aux_mc  = QuantumRegister(16, 'aux_mc')
    cr      = ClassicalRegister(16, 'c')
    qc      = QuantumCircuit(state, aux_sub, aux_mc, cr)

    sbox     = build_sbox_nibble()
    sbox_inv = build_sbox_inv_nibble()

    def apply_sub(qc):
        """SubNibbles reversible via truco de Bennett."""
        # 1. Computar SBOX en aux_sub
        for n in range(4):
            qin  = [state[n*4+i]    for i in range(4)]
            qout = [aux_sub[n*4+i]  for i in range(4)]
            qc.append(sbox, qin+qout)
        # 2. SWAP: state <-> aux_sub  (state = SBOX(input), aux_sub = input)
        for i in range(16): qc.swap(state[i], aux_sub[i])
        # 3. Descomputar: SBOX_INV limpia aux_sub a |0⟩
        for n in range(4):
            qin  = [state[n*4+i]    for i in range(4)]
            qout = [aux_sub[n*4+i]  for i in range(4)]
            qc.append(sbox_inv, qin+qout)

    def apply_sr(qc):
        """ShiftRows: SWAP(state[8..11], state[0..3])."""
        for i in range(4): qc.swap(state[8+i], state[0+i])

    def apply_mc(qc):
        """
        MixColumns con mapeo de columnas correcto post-ShiftRows.
        Post-SR: state[12..15]=n0, state[8..11]=n1, state[4..7]=n2, state[0..3]=n3
        """
        na0 = [state[12+i]  for i in range(4)]
        nb0 = [state[8+i]   for i in range(4)]
        oa0 = [aux_mc[12+i] for i in range(4)]
        ob0 = [aux_mc[8+i]  for i in range(4)]
        for i in range(4): qc.cx(na0[i], oa0[i])
        apply_mult4(qc, nb0, oa0)
        for i in range(4): qc.cx(nb0[i], ob0[i])
        apply_mult4(qc, na0, ob0)

        na1 = [state[4+i]   for i in range(4)]
        nb1 = [state[0+i]   for i in range(4)]
        oa1 = [aux_mc[4+i]  for i in range(4)]
        ob1 = [aux_mc[0+i]  for i in range(4)]
        for i in range(4): qc.cx(na1[i], oa1[i])
        apply_mult4(qc, nb1, oa1)
        for i in range(4): qc.cx(nb1[i], ob1[i])
        apply_mult4(qc, na1, ob1)

        for i in range(16): qc.swap(aux_mc[i], state[i])

    # ── Inicializar plaintext ─────────────────────────────────────────────────
    for i in range(16):
        if (P >> i) & 1: qc.x(state[i])
    qc.barrier()

    # ── AddRoundKey(k0) ───────────────────────────────────────────────────────
    for i in range(16):
        if (k0 >> i) & 1: qc.x(state[i])
    qc.barrier(label='ARK0')

    # ── Ronda 1 ───────────────────────────────────────────────────────────────
    apply_sub(qc); qc.barrier(label='Sub1')
    apply_sr(qc);  qc.barrier(label='SR1')
    apply_mc(qc);  qc.barrier(label='MC1')
    for i in range(16):
        if (k1 >> i) & 1: qc.x(state[i])
    qc.barrier(label='ARK1')

    # ── Ronda 2 (sin MixColumns) ──────────────────────────────────────────────
    apply_sub(qc); qc.barrier(label='Sub2')
    apply_sr(qc);  qc.barrier(label='SR2')
    for i in range(16):
        if (k2 >> i) & 1: qc.x(state[i])
    qc.barrier(label='ARK2')

    qc.measure(state, cr)
    return qc

if __name__ == "__main__":
    qc_info = build_saes_circuit(0xDF64, 0x544A)
    print(f"Circuito S-AES:")
    print(f"  Qubits     : {qc_info.num_qubits}")
    print(f"  Profundidad: {qc_info.decompose().depth()}")
    print(f"  Puertas    : {qc_info.size()}")

    print("=== Verificación ===")
    sim = AerSimulator(method='matrix_product_state')

    P_test = 0xDF64
    K_test = 0x544A

    print("=" * 52)
    print(" COMPARACIÓN S-AES: CLÁSICO vs CUÁNTICO")
    print("=" * 52)
    print(f"Texto Plano (P) : {P_test:#06x} = {P_test:016b}")
    print(f"Clave (K)       : {K_test:#06x} = {K_test:016b}\n")

    # Clásico
    t0 = time.time()
    C_c = saes_encrypt(P_test, K_test)
    t_c = time.time() - t0
    print(f"--- RESULTADO CLÁSICO ---")
    print(f"Cifrado (C)     : {C_c:#06x} = {C_c:016b}")
    print(f"Tiempo          : {t_c:.6f} s\n")

    # Cuántico
    print("--- RESULTADO CUÁNTICO ---")
    qc = build_saes_circuit(P_test, K_test)
    print(f"Qubits totales  : {qc.num_qubits}")
    qct = transpile(qc, sim)
    t0 = time.time()
    result = sim.run(qct, shots=1).result()
    t_q = time.time() - t0
    C_q = int(list(result.get_counts().keys())[0].split(' ')[-1], 2)
    print(f"Cifrado (C)     : {C_q:#06x} = {C_q:016b}")
    print(f"Tiempo          : {t_q:.6f} s\n")

    print("=" * 52)
    if C_c == C_q:
        print("(Ok) ÉXITO - Cuántico coincide con clásico")
    else:
        print("(Fallo) DISCREPANCIA")
    print("=" * 52)

    print("=== Verificación multi-vector ===")
    test_vectors = [
        (0xDF64, 0x544A),
        (0x0000, 0x0000),
        (0xFFFF, 0xFFFF),
        (0x1234, 0xABCD),
        (0x2D55, 0x4A61),
        (0xABCD, 0x1234),
        (0x6F6B, 0xCAFE),
    ]
    all_ok = True
    for P, K in test_vectors:
        C_c = saes_encrypt(P, K)
        qc  = build_saes_circuit(P, K)
        qct = transpile(qc, sim)
        C_q = int(list(sim.run(qct, shots=1).result().get_counts()
                    .keys())[0].split(' ')[-1], 2)
        ok  = C_c == C_q
        if not ok: all_ok = False
        print(f"  P={P:#06x}, K={K:#06x} -> C={C_c:#06x}  {'(Ok)' if ok else '(Fallo)'}")

    print(f"\n{'(Ok) Todos los test vectors correctos' if all_ok else '(Fallo) Hay discrepancias'}")