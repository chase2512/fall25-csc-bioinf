import time
import numpy as np
from biotite.sequence.phylo import Tree, TreeNode, upgma, neighbor_joining

def test_distances():
    """Test distance calculations between nodes"""
    leaf1 = TreeNode(index=0)
    leaf2 = TreeNode(index=1)
    leaf3 = TreeNode(index=2)
    inter = TreeNode([leaf1, leaf2], [5.0, 7.0])
    root = TreeNode([inter, leaf3], [3.0, 10.0])
    tree = Tree(root)
    
    assert abs(tree.get_distance(0, 1) - 12.0) < 0.001
    assert abs(tree.get_distance(0, 2) - 18.0) < 0.001
    assert abs(tree.get_distance(1, 2) - 20.0) < 0.001

def test_upgma():
    """Test UPGMA algorithm"""
    distances = np.array([
        [0, 1, 7, 7, 9],
        [1, 0, 7, 6, 8],
        [7, 7, 0, 2, 4],
        [7, 6, 2, 0, 3],
        [9, 8, 4, 3, 0],
    ], dtype=float)
    
    tree = upgma(distances)
    assert len(tree) == 5
    assert tree.root is not None

def test_neighbor_joining():
    """Test Neighbor Joining algorithm"""
    dist = np.array([
        [0, 5, 4, 7, 6, 8],
        [5, 0, 7, 10, 9, 11],
        [4, 7, 0, 7, 6, 8],
        [7, 10, 7, 0, 5, 9],
        [6, 9, 6, 5, 0, 8],
        [8, 11, 8, 9, 8, 0],
    ], dtype=float)

    test_tree = neighbor_joining(dist)
    
    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )
    
    assert test_tree == ref_tree

if __name__ == "__main__":
    start = time.time()
    
    test_distances()
    test_upgma()
    test_neighbor_joining()
    
    end = time.time()
    elapsed_ms = int((end - start) * 1000)
    
    # Ensure at least 1ms is reported
    if elapsed_ms == 0:
        elapsed_ms = 1
    
    print(f"python      {elapsed_ms}ms")