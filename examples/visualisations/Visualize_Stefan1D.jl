using ExactFieldSolutions, Plots, StaticArrays, Printf

function main()

    # Define domain
    Ny  = 51
    y   = LinRange(0, 4, Ny)

    # Allocate arrays
    t   = [0.1 1e3 1e4 1e5]
    T   = zeros(length(t),Ny)

    # Evaluate solution
    for it in eachindex(t)
        for j = 1:Ny
            X = @SVector[y[j]; t[it]]
            sol = Diffusion1D_StefanProblem(X)
            T[it,j] = sol.T
        end
    end

    # Visualise
    p1 = plot(xlabel = "T", ylabel= "y", yflip=true)
    p1 = plot!(T[1,:], y, title=@sprintf("T @ t = %1.2e", t[1]), label="T₀")
    p1 = plot!(T[1,:], y, title=@sprintf("T @ t = %1.2e", t[1]), label="T")
    p2 = plot(xlabel = "T", ylabel= "y", yflip=true)
    p2 = plot!(T[1,:], y, title=@sprintf("T @ t = %1.2e", t[1]), label="T₀")
    p2 = plot!(T[2,:], y, title=@sprintf("T @ t = %1.2e", t[2]), label="T")
    p3 = plot(xlabel = "T", ylabel= "y", yflip=true)
    p3 = plot!(T[1,:], y, title=@sprintf("T @ t = %1.2e", t[1]), label="T₀")
    p3 = plot!(T[3,:], y, title=@sprintf("T @ t = %1.2e", t[3]), label="T")
    p4 = plot(xlabel = "T", ylabel= "y", yflip=true)
    p4 = plot!(T[1,:], y, title=@sprintf("T @ t = %1.2e", t[1]), label="T₀")
    p4 = plot!(T[4,:], y, title=@sprintf("T @ t = %1.2e", t[4]), label="T")

    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 

end

main()
