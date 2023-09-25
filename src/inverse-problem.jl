function optimise_particulate(target_S::DiscreteStructureFactor{Dim}, specie::Specie{Dim};
        correlation_length::Float64 = 4.0, 
        numberofparticles::Int = 200,
        optimoptions::Optim.Options{Float64} = Optim.Options(f_tol = 1e-5, iterations = 10)
    ) where Dim

    # Define two regions R1 and R2 to place the particles
    cell_number = (numberofparticles)^(1/Dim) |> round

    cell_volume = 1 / number_density(specie);
    cell_length = cell_volume ^ (1/Dim)
    box_length = cell_length * cell_number

    if Dim != 2 
        error("the box regions only written for 2D, where as Dim = $Dim")
    end    

    reg1 = Box([
        [-box_length/2, -box_length/2],
        [box_length/2, box_length/2]
    ]);

    cell_number_boundary = ceil(correlation_length/(2 * cell_length))
    boundary_length = cell_number_boundary * cell_length

    reg2 = Box([
        [-box_length/2 - boundary_length, -box_length/2 - boundary_length],
        [box_length/2 + boundary_length, box_length/2 + boundary_length]
    ]);

    particles = periodic_particles(reg2,specie; random_perturbation = true)
    particles1 = filter(p -> p ⊆ reg1, particles)  

    S = structure_factor(particles, target_S.k; inner_box = reg1)

    points = origin.(particles)
    x0 = Vector{Float64}(vcat(points...))

    function objective(x)
        points = Iterators.partition(x,Dim) |> collect
        S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)

        return sum(abs2.(S.S - target_S.S))
    end

    ks = abs.(target_S.k)

    function objective_g!(G, x)
        points = Iterators.partition(x,Dim) |> collect
        points1 = filter(p -> p ∈ reg1, points)
        J1 = length(points1)  

        S = DiscreteStructureFactor(points, target_S.k; inner_box = reg1)
        dS = S.S - target_S.S # add quadradture weight w here if wanted

        G[:] = vcat(map(points) do p
            sum(
                if p1 == p
                    zeros(Float64, Dim)
                else     
                    Rji = p - p1
                    nRji = norm(Rji)
                    
                    Gji = (p ∈ reg1 ? 2 : 1) * sum(diffbesselj.(0, ks .* nRji) .* dS) / J1
                    2 * Gji * Rji / nRji
                end    
            for p1 in points1)
        end...)
    end    

    res = optimize(objective, objective_g!, x0, LBFGS(), 
        optimoptions
    )

    println("The result of the global step was:") 
    show(res)

    x0 = res.minimizer
    minf = res.minimum
    A = 10 * minf + 0.1

    function penalise(x)
        points = Iterators.partition(x,Dim) |> collect
        α = 1 / (2*outer_radius(specie))^2

        return A * sum(
            sum(
                begin
                    R12 = sum(abs2.(p1 - p2)) * α
                    if 0 < R12 < 1.0 
                        exp(- 4 * R12) 
                    else 0.0
                    end
                end        
            for p2 in points)
        for p1 in points)
    end
    
    function penalise_g!(G, x)
        points = Iterators.partition(x,Dim) |> collect
        α = 1 / (2*outer_radius(specie))^2

        G[:] =  - (A * 16 * α) .* vcat(
            map(points) do p2
                sum(
                    begin
                        R21 = sum(abs2.(p2 - p1)) .* α
                        if 0 < R21 < 1.0 
                            exp(- 4 * R21) .* (p2 - p1)
                        else zeros(Float64, Dim)
                        end
                    end        
                for p1 in points)
            end...)
    end    

    # refine the optimisation to not allow particles to overlap
    f(x) = objective(x) + penalise(x)
    g(G, x) = objective_g!(G, x) + penalise_g!(G, x)

    res = optimize(f, g, x0, LBFGS(), 
        optimoptions
    )

    # x0 = res.minimizer
    res
end