import torch

def node_ptr(batch):
    """Pointer vector for node-level tensors (replaces __slices__['x'] etc.)."""
    return batch.ptr            # shape: [num_graphs+1]

def edge_ptr(batch):
    """Pointer vector for edge-level tensors (bonds, edge_index, …)."""
    #   graph-id of each edge = graph-id of its source node
    edge_batch = batch.batch[batch.edge_index[0]]
    num_graphs = batch.ptr.numel() - 1
    counts = torch.bincount(edge_batch, minlength=num_graphs)
    return torch.cat([batch.ptr.new_zeros(1), counts.cumsum(0)])  # same semantics
