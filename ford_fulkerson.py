import copy

def find_path(graph: list[list[int]], source: int, target: int):
    """Return a path -- any path -- from source to target in the graph"""

    # Initialize return item
    path: list[int] = None

    # Make sure inputs are ok
    if graph is not None:
        n = len(graph)
        if n > 0 and (0 <= source < n) and (0 <= target < n):

            # Initialize DFS tools
            no_edge: int = graph[0][0]  # absence of edge
            marked: list[int] = [source]  # vertices already processed
            found: bool = False  # Flags detection of path

            # What vertex to explore next and what is the path
            # to it. The information is stored as a tuple in
            # the form:
            #  (vertex, path_to_this_vertex)
            # with path_to_this_vertex being a list of the
            # vertices alonγ the path.
            stack: list[(int, list[int])] = [(source, [source])]

            while len(stack) > 0 and not found:
                # Explore the next vertex from the stack
                (u, path_from_source_to_u) = stack.pop()
                found = (u == target)
                if found:
                    # u is the end of the path, so we got what we are 
                    # looking for
                    path = path_from_source_to_u
                else:
                    # Explore the neighbors of u, hopefully one of them
                    # will get us a stop closer to the target vertex.
                    v: int = n - 1
                    while v >= 0:
                        if graph[u][v] != no_edge and v not in marked:
                            marked.append(v)
                            stack.append((v, path_from_source_to_u + [v]))
                        v -= 1
    return path


def reachability_of(s: int, G: list[list[int]]) -> list[int]:
    """
    DFS returns a list of vertices in the same component as/reachable from s.

    Args:
    s: source vertex
    G: graph

    Returns:
    A list of all vertices that are reachable from s (same component of graph).
    """
    no_edge = G[0][0] # local infinity/non-edge
    reach = [] # which vertices are reachable from s
    visit_next = [s] # which vertices to visit next

    while visit_next:  
        v = visit_next.pop()
        if v not in reach:
            reach.append(v)
            for u in range(len(G)):
                if G[v][u] != no_edge:
                    visit_next.append(u)
    return reach


def ford_fulkerson(graph, source, target):
    """
    Finds the maximum flow from source to target and the minimum
    cut (smallest set of edges which disrupts source->target flow)

    Returns
    _______
    
    max_flow : int
        Maximum flow of path from source to target
    min_cut : [(int, int)...]
        Edges which cut off flow from source to target
    """
    residual = copy.deepcopy(graph)
    max_flow = 0
    min_cut = []

    path = find_path(residual, source, target)

    # while any path between source and target exists, find that paths maximum flow
    while path is not None:
        path_flow = float('inf')
        # find maximum flow of current path
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i+1]
            if residual[u][v] < path_flow:
                path_flow = residual[u][v]

        # maximum flow of every path is accumulated onto max_flow
        max_flow += path_flow

        # update residual graph's flows
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]
            residual[u][v] -= path_flow
            residual[v][u] += path_flow
        
        # find the next path
        path = find_path(residual, source, target)

    # we got max_flow; now for min_cut
    # no path exists between source and target in residual graph anymore, so
    # they're in different components in it. min_cut edges
    reachable = reachability_of(source, residual)

    # if an edge between the reachable and unreachable component exists,
    # it is a min_cut edge
    for u in reachable:
        for v in range (len(graph)):
            if v not in reachable and graph[u][v] > 0:
                min_cut.append((u,v))
    
    return max_flow, min_cut





# Test graph

G = [  #  A   B   C   D   E
    [0, 20, 0, 0, 0],  # A
    [0, 0, 5, 6, 0],  # B
    [0, 0, 0, 3, 7],  # C
    [0, 0, 0, 0, 8],  # D
    [0, 0, 0, 0, 0],  # E
]



src, tgt = 0, 3

flow, cuts = ford_fulkerson(G, src, tgt)

print(f"[{src} --> {tgt}] Flow: {flow}  Edges to cut: {cuts}")
