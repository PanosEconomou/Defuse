module DefuseMain

import Defuse: foldedMinimalModularData

function julia_main()::Cint
    try
        args = ARGS

        p  = length(args) >= 1 ? parse(Int, args[1]) : 5
        pp = length(args) >= 2 ? parse(Int, args[2]) : 4

        s, t, z, sv = foldedMinimalModularData(p,pp)

        display(s)
        return 0
    catch e
        println(stderr,"Error: ",e)
        return 1
    end
end

end

