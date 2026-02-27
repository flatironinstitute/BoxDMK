using Test

include("../src/BoxDMK.jl")
using .BoxDMK

@test isdefined(BoxDMK, :vol_tree_mem!)
@test isdefined(BoxDMK, :vol_tree_build!)
@test isdefined(BoxDMK, :bdmk!)
