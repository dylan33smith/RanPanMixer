import numpy as np, itertools, math
np.set_printoptions(precision=4, suppress=True)
NE,R = 10_000, 1.26
NAMES=['D1','D2','D3','D4']; K=4
H = np.array([[0,1, 1,0,0, 1,1],
              [0,1, 0,1,1, 0,0],
              [1,0, 1,0,0, 1,0],
              [1,0, 0,1,1, 0,1]])
target = np.array([0,1, 1,0,0, 0,1])
BLK = {0:[0,1], 1:[2,3,4], 2:[5,6]}
cMv  = np.array([0.0, 0.00002, 0.00012, 0.00014, 0.00016, 0.00066, 0.00068])
def ls(dx, n):
    d = dx*NE*R; e = math.exp(-d/n)
    return e + (1-e)/n, (1-e)/n
def Amat(dx, n=K):
    s,q = ls(dx,n); M=np.full((n,n),q); np.fill_diagonal(M,s); return M
beta_v = np.full(7, 1/7)                       # equal per-variant weight, sums to 1
phi_v  = (H == target).astype(float)           # K x 7 : per-variant agreement
TAU=0.5; ETA=math.atanh(TAU)

print("=== the toy, model B: ONE chain, units = BLOCKS ===")
print("        v1 v2 | v3 v4 v5 | v6 v7")
for j in range(K): print(f"  {NAMES[j]}    {H[j,0]}  {H[j,1]}  |  {H[j,2]}  {H[j,3]}  {H[j,4]}  |  {H[j,5]}  {H[j,6]}")
print(f"  target {target[0]}  {target[1]}  |  {target[2]}  {target[3]}  {target[4]}  |  {target[5]}  {target[6]}")

gapB = [cMv[BLK[t+1][0]] - cMv[BLK[t][-1]] for t in range(2)]
print(f"\ninter-block genetic distances: block1->2 {gapB[0]:.5f} cM, block2->3 {gapB[1]:.5f} cM")
for i,g in enumerate(gapB):
    s,q = ls(g,K); print(f"  gap {i+1}: d = {g*NE*R:7.4f}   p_stay {s:.4f}   q(each other) {q:.4f}   total switch {(K-1)*q:.4f}")

beta_b = np.array([beta_v[BLK[t]].sum() for t in range(3)])
phi_b  = np.array([[phi_v[j,BLK[t]].mean() for j in range(K)] for t in range(3)])
psi_b  = np.exp(ETA*beta_b[:,None]*phi_b)
print(f"\nbeta (block share of utility weight): {beta_b}   sum {beta_b.sum():.4f}")
print("phi (fraction of block's variants the donor matches):"); print(phi_b)
print(f"psi = exp(eta*beta*phi), tau={TAU}, eta_tau={ETA:.4f}:"); print(psi_b)

# ---- forward filter, model B
T=3; alpha=np.zeros((T,K)); c=np.zeros(T)
rho=np.full(K,1/K)
r=rho*psi_b[0]; c[0]=r.sum(); alpha[0]=r/c[0]
for t in range(1,T):
    pred = alpha[t-1] @ Amat(gapB[t-1]); r = psi_b[t]*pred
    c[t]=r.sum(); alpha[t]=r/c[t]
print("\n=== forward pass ===")
for t in range(T): print(f"  alpha_{t+1} = {alpha[t]}   c_{t+1} = {c[t]:.6f}")
logZ=np.log(c).sum(); print(f"  log Z_p = {logZ:.6f}   in [0, {ETA:.4f}] -> {0<=logZ<=ETA}")

# ---- backward sample
rng=np.random.default_rng(11); z=[0]*T
z[T-1]=rng.choice(K,p=alpha[T-1]); print(f"\n=== backward sampling ===\n  z_3 ~ {alpha[T-1]} -> {NAMES[z[2]]}")
for t in range(T-2,-1,-1):
    w=alpha[t]*Amat(gapB[t])[:,z[t+1]]; w=w/w.sum(); z[t]=rng.choice(K,p=w)
    print(f"  z_{t+1} | z_{t+2}={NAMES[z[t+1]]}:  {w} -> {NAMES[z[t]]}")
rel=np.concatenate([H[z[t],BLK[t]] for t in range(T)])
print(f"  path {' -> '.join(NAMES[j] for j in z)}   released {rel}   target {target}")
print(f"  u = {sum(beta_b[t]*phi_b[t,z[t]] for t in range(T)):.4f}")

# ---- exact enumeration
paths=list(itertools.product(range(K),repeat=T)); w=[]
for zz in paths:
    p_=rho[zz[0]]*psi_b[0,zz[0]]
    for t in range(1,T): p_*=Amat(gapB[t-1])[zz[t-1],zz[t]]*psi_b[t,zz[t]]
    w.append(p_)
w=np.array(w); Q=w/w.sum()
print(f"\n=== enumeration over {len(paths)} paths ===")
print(f"  Z_p direct {w.sum():.6f} vs forward {math.exp(logZ):.6f}  agree {np.isclose(w.sum(),math.exp(logZ))}")
U=np.array([sum(beta_b[t]*phi_b[t,zz[t]] for t in range(T)) for zz in paths])
base=[]
for zz in paths:
    p_=rho[zz[0]]
    for t in range(1,T): p_*=Amat(gapB[t-1])[zz[t-1],zz[t]]
    base.append(p_)
base=np.array(base); base/=base.sum()
print(f"  E[u] baseline {(base*U).sum():.4f}  ->  tilted {(Q*U).sum():.4f}   gain {(Q*U).sum()-(base*U).sum():+.4f}")
for i in np.argsort(-Q)[:4]:
    print(f"    {' -> '.join(NAMES[j] for j in paths[i]):<14} R_D={base[i]:.4f}  Q_p={Q[i]:.4f}  u={U[i]:.4f}")

# ================= the tilt vanishing as T grows (model B, tiled) =================
def Etilt(Tblocks):
    pat=[t%3 for t in range(Tblocks)]
    sizes=np.array([len(BLK[p]) for p in pat]); bw=sizes/sizes.sum()
    ph=np.array([phi_b[p] for p in pat])
    ps=np.exp(ETA*bw[:,None]*ph)
    gaps=[gapB[t%2] for t in range(Tblocks-1)]
    a=np.zeros((Tblocks,K)); r=rho*ps[0]; a[0]=r/r.sum()
    for t in range(1,Tblocks):
        r=ps[t]*(a[t-1]@Amat(gaps[t-1])); a[t]=r/r.sum()
    b=np.ones((Tblocks,K))
    for t in range(Tblocks-2,-1,-1):
        v=Amat(gaps[t])@(ps[t+1]*b[t+1]); b[t]=v/v.sum()
    g=a*b; g/=g.sum(axis=1,keepdims=True)
    EQ=(bw*(g*ph).sum(axis=1)).sum(); ER=(bw*ph.mean(axis=1)).sum()
    return ER,EQ,ps.max()
print("\n\n=== HOW THE TILT VANISHES AS THE CHAIN GROWS (model B, pattern tiled) ===")
print(f"{'T blocks':>10}{'max psi':>12}{'E[u] baseline':>16}{'E[u] tilted':>14}{'gain':>12}")
for Tb in [3,30,300,3000,30000,100757]:
    Tb=Tb-(Tb%3) if Tb%3 else Tb
    ER,EQ,mx=Etilt(Tb)
    print(f"{Tb:>10,}{mx:>12.7f}{ER:>16.4f}{EQ:>14.6f}{EQ-ER:>+12.2e}")

# ================= same toy under models A and C =================
def chain_over_variants(inter_block_uniform):
    gaps=[cMv[i+1]-cMv[i] for i in range(6)]
    bnd={1,4}                                    # steps 1->2 and 4->5 cross a block boundary
    ps=np.exp(ETA*beta_v[None,:]*phi_v).T        # 7 x K
    P=list(itertools.product(range(K),repeat=7))
    wq=[];wr=[]
    for zz in P:
        q_=rho[zz[0]]*ps[0,zz[0]]; r_=rho[zz[0]]
        for t in range(1,7):
            M=np.full((K,K),1/K) if (t-1 in bnd and inter_block_uniform) else Amat(gaps[t-1])
            q_*=M[zz[t-1],zz[t]]*ps[t,zz[t]]; r_*=M[zz[t-1],zz[t]]
        wq.append(q_); wr.append(r_)
    wq=np.array(wq); wq/=wq.sum(); wr=np.array(wr); wr/=wr.sum()
    U=np.array([ (np.array([phi_v[zz[t],t] for t in range(7)])*beta_v).sum() for zz in P])
    sws=np.array([sum(zz[t]!=zz[t+1] for t in range(6)) for zz in P])
    return (wr*U).sum(),(wq*U).sum(),(wq*sws).sum(),(wr*sws).sum()
print("\n=== the SAME toy under all three models (tau=0.5, identical u) ===")
print(f"{'model':<46}{'E[u] base':>11}{'E[u] tilt':>11}{'gain':>9}{'E[switch]':>11}")
rA=chain_over_variants(True); rC=chain_over_variants(False)
UB=(base*U).sum(); QB=(Q*U).sum()
swB=(Q*np.array([sum(zz[t]!=zz[t+1] for t in range(2)) for zz in paths])).sum()
print(f"{'A. per-block chains, uniform reset at boundaries':<46}{rA[0]:>11.4f}{rA[1]:>11.4f}{rA[1]-rA[0]:>+9.4f}{rA[2]:>11.3f}")
print(f"{'B. one chain, units = blocks (blocks atomic)':<46}{UB:>11.4f}{QB:>11.4f}{QB-UB:>+9.4f}{swB:>11.3f}")
print(f"{'C. one chain, units = variants':<46}{rC[0]:>11.4f}{rC[1]:>11.4f}{rC[1]-rC[0]:>+9.4f}{rC[2]:>11.3f}")
print(f"\n  (E[switch] for A and C counts donor changes across 6 variant steps;")
print(f"   for B across 2 block steps -- B cannot switch inside a block by construction)")
