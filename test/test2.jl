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
B = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/mon_B.xyz",  skipstart = 2)
dma = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/density_matrix_A.txt")
dmb = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/density_matrix_B.txt")

#eth_mon_A = gto.M(atom = A, symmetry = true, basis = use_basis, spin = 0, verbose = 5)
#rks_mon_A = dft.RKS(eth_mon_A)
#rks_mon_A.xc = use_functional
#@time rks_mon_A.kernel(dm0=dma)

#MO_mon_a = rks_mon_A.mo_coeff[:,8:9]
#writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/frontier_orbitals_a.txt", MO_mon_a)

#eth_mon_B = gto.M(atom = B, symmetry = true, basis = use_basis, spin = 0, verbose = 5)
#rks_mon_B = dft.RKS(eth_mon_B)
#rks_mon_B.xc = use_functional
#@time rks_mon_B.kernel(dm0=dmb)

#MO_mon_b = rks_mon_B.mo_coeff[:,8:9]
#writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/frontier_orbitals_b.txt", MO_mon_b)

MO_mon_a = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/frontier_orbitals_a.txt")
MO_mon_b = readdlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/frontier_orbitals_b.txt")

eth_dim = gto.M(atom = [A ; B], symmetry = true, basis = use_basis, spin = 0, verbose = 5)
n = size(dma)[1]
dmdim = [dma zeros(n,n); zeros(n,n) dmb];

rks_dim = dft.RKS(eth_dim)
rks_dim.xc = use_functional
@time rks_dim.kernel(dm0=dmdim)

#molden.from_mo(eth_dim,"eth_dim.molden", rks_dim.mo_coeff[:,16:17], ene=rks_dim.mo_energy[16:17])

AO_Overlap = eth_dim.intor_symmetric("int1e_ovlp");

println("Overlap Matrix determined...")

MOs_energy = Diagonal(rks_dim.mo_energy);
MOs_dim = rks_dim.mo_coeff;
MOs_mon = [MO_mon_a zeros(n,2); zeros(n,2) MO_mon_b];

println("Monomer direct sum constructed...")

cf = 27.2;

E = AO_Overlap * MOs_dim * MOs_energy * inv(MOs_dim)
H = transpose(MOs_mon) * E * MOs_mon;
S = transpose(MOs_mon) * AO_Overlap * MOs_mon;
Slowdin = inv(sqrt(S));

Heff = cf * Slowdin * H * Slowdin

g = transpose(MOs_mon) * AO_Overlap * MOs_dim

s = g * transpose(g)
slowdin = inv(sqrt(s))

h = g * MOs_energy * transpose(g)
heff = cf * slowdin * h * slowdin

v2e = ao2mo.kernel(eth_dim,MOs_mon,intor="int2e",aosym = 1)
V2e = cf*reshape(v2e,(4,4,4,4))

writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/t1_1storder_couplings.txt", heff)
writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/t1_2ndorder_couplings.txt", V2e)

#dmA = rks_mon_A.make_rdm1()
#writedlm("/Users/alex/Documents/VSCode_Julia/MyPhD/test/density_matrix_B.txt",dmA)