
module OrderEquivalence

using DelimitedFiles

#export generatecluster

# To Do --> Project: Enforce Intermonomer Atom Order Equivalence
#           (1) Read input coordinates ** DONE **
#           (2) Compute orientational details for all other monomers based on target
#           (3) Identify equivalent atoms
#           (4) Reorder coordinate lists for other monomers

# function (1): Reader [Coordinates]

struct Atom                                  # no need to mutate as monomer is array of structs!
    label::String                            # atom type
    r::Vector{Float64}                       # position vector
end                                          # choose atom as struct to include charges, momenta, fixes unclear molecule size

# Constructor? Seems to do fine with default. Could potentially rewrite code with struct cluster and inner constructor.

function generatecluster(path::String)           # path to folder containing coordinate files for each monomer in cluster

    cluster = Vector{Atom}[]                     # cluster = [ monomera monomerb ... ] for monomera = [atom1; atom2; ... ]

    for i in readdir(path; join=true)            # looping over files in folder
    
        data = readdlm(i,  skipstart = 3)        # read the coordinates into an array called data, skipping header lines
    
        monomer = [Atom(data[j,1],data[j,2:4]) for j in axes(data,1)]
        push!(cluster,monomer)
    end

    return cluster

end

# function (2): Orientor [Monomers|Target]

function localcoordinate(Cluster::Vector{Vector{Atom}},Technique::String="SVD")

    # Isolate Perylene Diimide Core

    

    # Construct local coordinate basis with three carbons

    

    # Construct local coordinate basis from moment of inertia

end

function orientor(Target::Vector{Atom},Monomer::Vector{Atom})

end

# function (3): Equivocator
# function (4): Shuffler

end