using Gadfly
using Cairo
using Fontconfig

type Coordinate
    x :: Float64
    y :: Float64
end

function input_coordinates(filename::String)
    coords = Array(Coordinate, 0)
    open(filename, "r") do file
        for line in eachline(file)
            coord_ary = split(strip(line, ['\n', ',']), ' ')
            coord = Coordinate(float(coord_ary[1]), float(coord_ary[2]))
            push!(coords, coord)
        end
    end

    return coords
end

function calc_distance{N <: Integer}(path::Array{N}, coords::Array{Coordinate})
    distance = 0
    for i in 1:(length(coords) - 1)
        present_idx, next_idx = path[i], path[i + 1]
        distance += sqrt((coords[next_idx].x - coords[present_idx].x)^2
                         + (coords[next_idx].y - coords[present_idx].y)^2)
    end

    return distance
end

function random_path{N <: Integer}(n::N)
    path = collect(1:(n - 1))
    path = shuffle(path)
    push!(path, path[1])
    return path
end

function random_path(coords::Array{Coordinate})
    n = length(coords)
    random_path(n)
end

function metropolis{N <: Integer, R <: Real}(current_path::Array{N},
                                             coords::Array{Coordinate},
                                             temp::R)
    current_distance = calc_distance(current_path, coords)
    candidate_path = copy(current_path)
    n = length(current_path)

    swap_cities_indices = Array(Any, 2)
    while true
        swap_cities_indices = rand(2:(length(current_path) - 1), 2)
        if swap_cities_indices[1] != swap_cities_indices[2]
            break
        end
    end
    candidate_path[swap_cities_indices[1]], candidate_path[swap_cities_indices[2]] =
        candidate_path[swap_cities_indices[2]], candidate_path[swap_cities_indices[1]]
    candidate_distance = calc_distance(candidate_path, coords)

    Δd = candidate_distance - current_distance
    if Δd < 0
        return candidate_path, candidate_distance
    end

    prob = exp(-Δd / temp)
    if rand() <= prob
        return candidate_path, candidate_distance
    else
        return current_path, current_distance
    end
end

function anneal{N <: Integer}(path::Array{N}, coords::Array{Coordinate}; n_iter::N=100000,
                              pmelt=0.7, tgt=0.1, stagfactor=0.05, procplt::Bool=false)
    n_cities = length(path)
    initial_distance = calc_distance(path, coords)
    min_distance, max_distance = initial_distance, initial_distance

    optimized_distances = Array(Float64, 0)
    distances = Array(Float64, 0)
    push!(optimized_distances, initial_distance)
    push!(distances, initial_distance)
    for i in 1:(max(0.01 * n_cities, 2))
        path_, disntance = metropolis(path, coords, typemax(Int64))
        if disntance < min_distance
            min_distance = disntance
        end
        if max_distance < disntance
            max_distance = disntance
        end
    end

    range = (max_distance - min_distance) * pmelt
    temp = tgt ^ (1 / n_iter)
    optimized_distance = initial_distance
    optimized_step = 1
    optimized_path = copy(path)
    path_ = copy(path)

    for i in 1:n_iter
        print("$(i) / $(n_iter) processing...\r")
        flush(STDOUT)

        dtemp = range * (temp ^ i)
        path_, distance = metropolis(path_, coords, dtemp)

        if distance < optimized_distance
            optimized_distance = distance
            optimized_path = copy(path_)
            optimized_step = i
        end
        push!(optimized_distances, optimized_distance)
        push!(distances, distance)

        # Reheat
        if i - optimized_step == stagfactor * n_iter
            temp = temp ^ (0.05 * i / n_iter)
        end
    end
    print("\n")

    if procplt
        white_panel = Theme(background_color="white")
        p = plot(layer(x=1:length(distances), y=distances, Geom.line,
                       Theme(default_color="gray")),
                 layer(x=1:length(distances), y=optimized_distances, Geom.line,
                       Theme(default_color="orange", line_width=2px)),
                 white_panel,
                 Coord.Cartesian(xmin=0, xmax=length(distances)))
        img = PNG("result.png", 1280px, 720px)
        draw(img, p)
    end

    return optimized_path
end

coords = input_coordinates("prefs.in")
init_path = random_path(coords)
path = anneal(init_path, coords, n_iter=100000, procplt=true)
