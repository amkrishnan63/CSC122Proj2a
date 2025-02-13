from collections import defaultdict, deque


def build_de_bruijn_graph(spectrum):
    graph = defaultdict(deque)
    for entry in spectrum:
        read_id, kmer = entry.split()
        graph[kmer[:-1]].append((kmer[1:], read_id))
    return graph

def eulerian_path(graph):
    in_degree = defaultdict(int)
    out_degree =  defaultdict(int)
    
    for node, neighbors in graph.items():
        out_degree[node] += len(neighbors)
        neighbors_list = (pair[0] for pair in neighbors)
        for neighbor in neighbors_list:
            in_degree[neighbor] += 1
    
    start_node = None
    for node in out_degree:
        if out_degree[node] - in_degree[node] == 1:
            start_node = node
            break
    if start_node == None:
        start_node = next(iter(graph))
    
    stack = [start_node]
    path = deque()
    read_order = deque()
    
    while stack:
        node = stack[-1]
        if len(graph[node]) != 0:
            next_kmer, read_id = graph[node].popleft()
            stack.append(next_kmer)
            read_order.append(read_id) 
        else:
            path.appendleft(stack.pop())
    
    return list(read_order)

def read_spectrum_from_fasta(file_path):
    spectrum = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                read_id = line.strip()
                kmer = next(f).strip()
                spectrum.append(f"{read_id} {kmer}")
    return spectrum

def main():
    input_file = "project2a_spectrum.fasta" 
    output_file = "predictions.csv"
    spectrum = read_spectrum_from_fasta(input_file)
    read_order = eulerian_path(build_de_bruijn_graph(spectrum))
    
    with open(output_file, 'w') as f:
        for read_id in read_order:
            f.write(read_id + "\n")  
    

if __name__ == "__main__":
    main()
