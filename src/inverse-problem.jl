function optimise_particulate(target_S::DiscreteStructureFactor{Dim}, specie::Specie{Dim};
        correlation_length::Float64 = 4.0, 
        numberofparticles::Int = 200,
        method = LBFGS(),
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
                    
                    Gji = (p ∈ reg1 ? 2 : 1) * sum(ks .* diffbesselj.(0, ks .* nRji) .* dS) / J1
                    2 * Gji * Rji / nRji
                end    
            for p1 in points1)
        end...)
    end    

    function penalise(x)
        points = Iterators.partition(x,Dim) |> collect
        α = 1 / (2*outer_radius(specie))^2

        return sum(
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
    
    function penalise_g(x)
        points = Iterators.partition(x,Dim) |> collect
        α = 1 / (2*outer_radius(specie))^2

        return - (16 * α) .* vcat(
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

    # add a soft constraint
    A0 = objective(x0) / numberofparticles
    f0(x) = objective(x) + A0 * penalise(x)
    # f(x) = objective(x)

    function g0!(G, x) 
        objective_g!(G, x) 
        G[:] = G + A0 .* penalise_g(x)
    end
    
    res1 = optimize(f0, g0!, x0, method, optimoptions)
    # res1 = optimize(objective, x0, optimoptions)

    println("The result of the global step was:") 
    show(res1)

    x0 = res1.minimizer
    minf = res1.minimum
    A = 10minf + 4A0

    # refine the optimisation to not allow particles to overlap
    f(x) = objective(x) + A * penalise(x)
    function g!(G, x) 
        objective_g!(G, x) 
        G[:] = G + A .* penalise_g(x)
    end

    res = optimize(f, g!, x0, method, optimoptions)

    # x0 = res.minimizer
    return res1, res
end