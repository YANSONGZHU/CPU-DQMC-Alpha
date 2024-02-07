include("ACFlow-main/src/ACFlow.jl")
using DelimitedFiles

file = open("tmp.toml")
txt=readlines(file, keep=true)
close(file)
L = 6
beta = 1/1
txt[2] = "finput = \"myGk L6 U-4 -6 -8/Gk L6 U-8 T1.0/kx=1 ky=1 kz=1.dat\"\n"
txt[8] = "ngrid  = 80\n"
txt[12] = "beta   = " * string(beta) * "\n"
io = open("ac.toml", "w")
for i in eachindex(txt)
    write(io, txt[i])
end
close(io)
ACFlow.setup_args("ac.toml")
ACFlow.read_param()
mesh, Aout, Gout = ACFlow.solve(ACFlow.read_data())
# Ak = readdlm("Aout.data")
# Akfile = "ak.dat"
# io = open(Akfile, "a")
# writedlm(io, Ak, '\t')
# close(io)