struct CyclicContainer{T} <: AbstractVector{T}
    c::Vector{T}
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.length(c::CyclicContainer) = length(c.c)
Base.size(c::CyclicContainer) = size(c.c)
Base.getindex(c::CyclicContainer, i::Int) = c.c[mod1(i, length(c.c))]
function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end
Base.iterate(c::CyclicContainer, i = 1) = iterate(c.c, i)
Base.getindex(c::CyclicContainer, i) = [c[j] for j in i]


# COLORS = [
# "#377e7f",
# "#8719cb",
# "#f09235",
# "crimson",
# "nayv" 
# ]

COLORS = [
"darkgreen",
"#f95738",
"#272727",
"#2cacc9",
"#885a89",
"#e09200",
"#dcddde"
]


# COLORS = [
#     "#1B1B1B",
#     "#6D44D0",
#     "#2CB3BF",
#     "#DA5210",
#     "#03502A",
#     "#866373",
#     "white",
#     "blue",
# ]

LINESTYLES = [
    '-', ':', "--", "-."
]


MARKERS = [
    :rtriangle,:ltriangle, :circle, :xcross, '□','∘', :xcross, :rtriangle, :ltriangle
]

MARKERSIZES = [4, 4, 2]

CCOLORS = CyclicContainer(COLORS)
CLINESTYLES = CyclicContainer(LINESTYLES)
CMARKERS = CyclicContainer(MARKERS)
CMARKERSIZES = CyclicContainer(MARKERSIZES)


export COLORS, CCOLORS, LINESTYLES, CLINESTYLES, MARKERS, CMARKERS, MARKERSIZES, CMARKERSIZES