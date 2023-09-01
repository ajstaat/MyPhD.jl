using PyCall
using Plots: plot, plot!
using LinearAlgebra

pd = pyimport("pyscf")
gto = pyimport("pyscf.gto")
dft = pyimport("pyscf.dft")

use_functional="wB97M_V"
use_basis="def2tzvppd"

eth_mon_A = gto.M(atom =
"C  -0.667121   0.000000   3.000000;
 C   0.667121  -0.000000   3.000000;
 H  -1.239234   0.924947   3.000000;
 H  -1.239234  -0.924947   3.000000;
 H   1.239234  -0.924947   3.000000;
 H   1.239234   0.924947   3.000000;
 X-C  -0.667121   0.000000   0.000000;
 X-C   0.667121   0.000000   0.000000;
 X-H  -1.239234   0.924947   0.000000;
 X-H  -1.239234  -0.924947   0.000000;
 X-H   1.239234  -0.924947   0.000000;
 X-H   1.239234   0.924947   0.000000",
symmetry = true, basis = use_basis, spin = 0, verbose = 5)

rks_mon_A = dft.RKS(eth_mon_A)
rks_mon_A.xc = use_functional
rks_mon_A.kernel()
