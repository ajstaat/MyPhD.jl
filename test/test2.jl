using PyCall
using Plots: plot, plot!
using LinearAlgebra
using DelimitedFiles

pd = pyimport("pyscf")
scf = pyimport("pyscf.scf")
gto = pyimport("pyscf.gto")
dft = pyimport("pyscf.dft")
molden = pyimport("pyscf.tools.molden")
ao2mo = pyimport("pyscf.ao2mo")
lib = pyimport("pyscf.lib")

lib.num_threads(2)

cf = 27.2;

use_functional="wB97XD"
use_basis="321g"

A = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/mon_A.xyz",  skipstart = 2)
dma = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/density_matrix_A.txt")

eth_mon_A = gto.M(atom = A, symmetry = true, basis = use_basis, spin = 0, verbose = 5)

rks_mon_A = dft.RKS(eth_mon_A)
rks_mon_A.xc = use_functional
@time rks_mon_A.kernel(dm0=dma)

#dmA = rks_mon_A.make_rdm1()

#writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/density_matrix_A.txt",dmA)
