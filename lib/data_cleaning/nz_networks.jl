# Download and prep New-Zealand food webs for analysis


NZ_foodwebs = DataFrame.(CSV.File.(glob("*.csv", 
                                        joinpath("data", "raw", "new_zealand", "adjacency_matrices")),
                                        drop=[1])) # the first column is row names
N_NZ = Binary.(NZ_foodwebs)
N_NZ = UnipartiteNetwork.(N_NZ, names.(NZ_foodwebs))