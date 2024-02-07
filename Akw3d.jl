include("ACFlow-main/src/ACFlow.jl")
using DelimitedFiles

file = open("tmp.toml")
txt=readlines(file, keep=true)
close(file)
L = 6
# kx = [1 2 3 4 4 4 4 4 4 4 3 2]
# ky = [1 1 1 1 2 3 4 4 4 4 3 2]
# kz = [1 1 1 1 1 1 1 2 3 4 3 2]
T = [2.0 1.6 1.2 1.0 0.8 0.6]
for i = 1:6
    beta = 1/T[i]
    for kx = 1:4
        for ky = 1:4
            for kz = 1:4
                txt[2] = "finput = \"myGk L6 U-4 -6 -8/Gk L6 U-8 T" * string(T[i]) * "/kx=" * string(kx) * " ky=" * string(ky) * " kz=" * string(kz) * ".dat\"\n"
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
                Ak = readdlm("Aout.data")
                Akfile = "myAk L6 U-4 -6 -8/Ak L6 U-8 T" * string(T[i]) * "/kx=" * string(kx) * " ky=" * string(ky) * " kz=" * string(kz) * ".dat"
                io = open(Akfile, "a")
                writedlm(io, Ak, '\t')
                close(io)
                rm("Aout.data")
                rm("chi2.data")
                rm("Gout.data")
                rm("model.data")
                rm("repr.data")
            end
        end
    end
end