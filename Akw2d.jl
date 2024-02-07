include("ACFlow-main/src/ACFlow.jl")
using DelimitedFiles

file = open("tmp.toml")
txt=readlines(file, keep=true)
close(file)
for kx = 1:7
    for ky = 1:7
        txt[2] = "finput = \"Gk L12 U0 T1.0/kx=" * string(kx) * " ky=" * string(ky) * ".dat\"\n"
        io = open("ac.toml", "w")
        for i in eachindex(txt)
            write(io, txt[i])
        end
        close(io)
        ACFlow.setup_args("ac.toml")
        ACFlow.read_param()
        mesh, Aout, Gout = ACFlow.solve(ACFlow.read_data())
        Ak = readdlm("Aout.data")
        Akfile = "Ak L12 U0 T1.0/kx=" * string(kx) * " ky=" * string(ky) * ".dat"
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