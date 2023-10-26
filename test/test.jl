using PyCall
using Plots: plot, plot!
using LinearAlgebra

pd = pyimport("pyscf")
scf = pyimport("pyscf.scf")
gto = pyimport("pyscf.gto")
dft = pyimport("pyscf.dft")
lib = pyimport("pyscf.lib")

lib.num_threads(2)

#use_functional="camb3lyp"
use_functional="wB97XD"
#use_functional="wB97M_V"
#use_basis="def2tzvppd"
use_basis="def2tzvppd"
#use_basis="631+g*"

eth_mon_A = gto.M(atom =
"C  -0.667121   0.000000   3.000000;
 C   0.667121  -0.000000   3.000000;
 H  -1.239234   0.924947   3.000000;
 H  -1.239234  -0.924947   3.000000;
 H   1.239234  -0.924947   3.000000;
 H   1.239234   0.924947   3.000000",#;
 #X-C  -0.667121   0.000000   0.000000;
 #X-C   0.667121   0.000000   0.000000;
 #X-H  -1.239234   0.924947   0.000000;
 #X-H  -1.239234  -0.924947   0.000000;
 #X-H   1.239234  -0.924947   0.000000;
 #X-H   1.239234   0.924947   0.000000",
symmetry = true, basis = use_basis, spin = 0, verbose = 5)

rks_mon_A = dft.RKS(eth_mon_A)
rks_mon_A.xc = use_functional
@time rks_mon_A.kernel()

dmA = rks_mon_A.make_rdm1()
#n = length(dmA[:,1]) ÷ 2;
#dmb = [dmA[n+1:end,n+1:end] transpose(dmA[1:n,n+1:end]); transpose(dmA[n+1:end,1:n]) dmA[1:n,1:n]];

eth_mon_B = gto.M(atom =
#"X-C  -0.667121   0.000000   3.000000;
# X-C   0.667121  -0.000000   3.000000;
# X-H  -1.239234   0.924947   3.000000;
# X-H  -1.239234  -0.924947   3.000000;
# X-H   1.239234  -0.924947   3.000000;
# X-H   1.239234   0.924947   3.000000;
" C  -0.667121   0.000000   0.000000;
 C   0.667121   0.000000   0.000000;
 H  -1.239234   0.924947   0.000000;
 H  -1.239234  -0.924947   0.000000;
 H   1.239234  -0.924947   0.000000;
 H   1.239234   0.924947   0.000000",
symmetry = true, basis = use_basis, spin = 0, verbose = 5)

rks_mon_B = dft.RKS(eth_mon_B)
rks_mon_B.xc = use_functional
@time rks_mon_B.kernel(dm0=dmA)

dmB = rks_mon_B.make_rdm1();
n = length(dmA[:,1])
dmDim = [dmA zeros(n,n); zeros(n,n) dmB];
#dmDim = dmA + dmB;

eth_dim = gto.M(atom =
"C  -0.667121   0.000000   3.000000;
 C   0.667121  -0.000000   3.000000;
 H  -1.239234   0.924947   3.000000;
 H  -1.239234  -0.924947   3.000000;
 H   1.239234  -0.924947   3.000000;
 H   1.239234   0.924947   3.000000;
 C  -0.667121   0.000000   0.000000;
 C   0.667121   0.000000   0.000000;
 H  -1.239234   0.924947   0.000000;
 H  -1.239234  -0.924947   0.000000;
 H   1.239234  -0.924947   0.000000;
 H   1.239234   0.924947   0.000000",
symmetry = true, basis = use_basis, spin = 0, verbose = 5)

# Q1. Altering Linear Dependence Tolerances? ** DONE ** --> Comparable speed, Canonical Orthogonalization
# Q2. Enforcing Partial Cholesky? ** DONE ** --> Comparable speed, Thresholds aggression equivalent to QChem default
# Q3. Counterpoise Correction = Better Initial Guess? --> Still 7 cycles for some reason.

# Alternative Initial Guesses:
# minao -> 9 cycles
# huckel -> 9 cycles (didn't obviously switch to something else/still stated minao)
# mod huckel -> 9 cycles
# vsap -> 9 cycles, wtf.

#m = length(dmA[:,1])
#dmDim = [dmA zeros(m,m); zeros(m,m) dmA]

rks_dim = dft.RKS(eth_dim)
scf.addons.remove_linear_dep_(rks_dim,threshold=1e-06,lindep=1e-06,cholesky_threshold=1e-06,force_pivoted_cholesky=true)
rks_dim.xc = use_functional
#rks_dim.init_guess = "vsap"
@time rks_dim.kernel()#dm0=dmDim)
