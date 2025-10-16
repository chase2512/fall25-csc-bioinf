# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.sequence.phylo"
__author__ = "Patrick Kunzmann"
__all__ = ["neighbor_joining"]

from tree_codon import Tree, TreeNode
from typing import List
import numpy as np


def neighbor_joining(distances: np.ndarray) -> Tree:
    """
    neighbor_join(distances)
    
    Perform hierarchical clustering using the
    *neighbor joining* algorithm.

    In contrast to UPGMA this algorithm does not assume a constant
    evolution rate. The resulting tree is considered to be unrooted.

    Parameters
    ----------
    distances : ndarray, shape=(n,n)
        Pairwise distance matrix.

    Returns
    -------
    tree : Tree
        A rooted tree. The `index` attribute in the leaf
        :class:`TreeNode` objects refer to the indices of `distances`.

    Raises
    ------
    ValueError
        If the distance matrix is not symmetric
        or if any matrix entry is below 0.
    
    Notes
    -----
    The created tree is binary except for the root node, that has three
    child notes

    Examples
    --------
    
    >>> distances = np.array([
    ...     [0, 1, 7, 7, 9],
    ...     [1, 0, 7, 6, 8],
    ...     [7, 7, 0, 2, 4],
    ...     [7, 6, 2, 0, 3],
    ...     [9, 8, 4, 3, 0],
    ... ])
    >>> tree = neighbor_joining(distances)
    >>> print(tree.to_newick(include_distance=False))
    (3,(2,(1,0)),4);
    """
    MAX_FLOAT = float(np.finfo(np.float32).max)

    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if distances.shape[0] < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    # Keep track on clustered indices
    nodes: List[TreeNode] = []
    for i in range(distances.shape[0]):
        nodes.append(TreeNode(index=i))
    
    # Indicates whether an index in the distance matrix has already been
    # clustered and the respective rows and columns can be ignored
    is_clustered_v: List[bool] = []
    for i in range(distances.shape[0]):
        is_clustered_v.append(False)
    
    n_rem_nodes = len(distances) - sum(is_clustered_v)
    
    # The divergence of a 'taxon'
    # describes the relative evolution rate
    divergence_v: List[float] = []
    for i in range(distances.shape[0]):
        divergence_v.append(0.0)
    
    # Triangular matrix for storing the divergence corrected distances
    corr_distances_v = np.zeros((distances.shape[0], distances.shape[0]), dtype=np.float32)
    
    distances_v = distances.astype(np.float32, copy=True)

    # Cluster indices
    # Exit loop via 'return'
    while True:
        # Calculate divergence
        for i in range(distances_v.shape[0]):
            if is_clustered_v[i]:
                continue
            dist_sum = 0.0
            for k in range(distances_v.shape[0]):
                if is_clustered_v[k]:
                    continue
                dist_sum += float(distances_v[i, k])
            divergence_v[i] = dist_sum
        
        # Calculate corrected distance matrix
        for i in range(distances_v.shape[0]):
            if is_clustered_v[i]:
                continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                corr_distances_v[i, j] = (
                    float(n_rem_nodes - 2) * float(distances_v[i, j])
                    - divergence_v[i] - divergence_v[j]
                )

        # Find minimum corrected distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        
        for i in range(corr_distances_v.shape[0]):
            if is_clustered_v[i]:
                continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                dist = float(corr_distances_v[i, j])
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        # Check if all nodes have been clustered
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            # -> exit loop
            break
        
        # Cluster the nodes with minimum distance
        # replacing the node at position i_min
        # leaving the node at position j_min empty
        # (is_clustered_v -> True)
        node_dist_i = 0.5 * (
            float(distances_v[i_min, j_min])
            + 1.0 / float(n_rem_nodes - 2) * (divergence_v[i_min] - divergence_v[j_min])
        )
        node_dist_j = 0.5 * (
            float(distances_v[i_min, j_min])
            + 1.0 / float(n_rem_nodes - 2) * (divergence_v[j_min] - divergence_v[i_min])
        )
        
        if n_rem_nodes > 3:
            # Clustering is not finished
            # -> Create a node with two children
            nodes[i_min] = TreeNode(
                [nodes[i_min], nodes[j_min]],
                [node_dist_i, node_dist_j]
            )
            # Mark position j_min as clustered
            nodes[j_min] = None
            is_clustered_v[j_min] = True
        else:
            # Clustering is finished
            # Combine last three nodes into root node
            # Find the index of the remaining one of the three nodes
            # (other than i_min and j_min)
            is_clustered_v[i_min] = True
            is_clustered_v[j_min] = True
            
            # The index of the remaining one
            k = -1
            for idx in range(len(is_clustered_v)):
                if not is_clustered_v[idx]:
                    k = idx
                    break
            
            node_dist_k = 0.5 * (
                float(distances_v[i_min, k]) + float(distances_v[j_min, k])
                - float(distances_v[i_min, j_min])
            )
            
            root = TreeNode(
                [nodes[i_min], nodes[j_min], nodes[k]],
                [node_dist_i, node_dist_j, node_dist_k]
            )
            # Clustering is finished -> put into tree and return
            return Tree(root)
        
        # Update distance matrix
        # Calculate distances of new node to all other nodes
        for k in range(distances_v.shape[0]):
            if not is_clustered_v[k] and k != i_min:
                dist = 0.5 * (
                    float(distances_v[i_min, k]) + float(distances_v[j_min, k])
                    - float(distances_v[i_min, j_min])
                )
                distances_v[i_min, k] = dist
                distances_v[k, i_min] = dist

        # Update the amount of remaining nodes
        n_rem_nodes = len(distances) - sum(is_clustered_v)
