import time
import numpy as np
from tree_codon import Tree, TreeNode
from upgma_codon import upgma
from nj_codon import neighbor_joining

def test_distances():
    """Test distance calculations between nodes"""
    leaf1 = TreeNode(index=0)
    leaf2 = TreeNode(index=1)
    leaf3 = TreeNode(index=2)
    inter = TreeNode(children=[leaf1, leaf2], distances=[5.0, 7.0])
    root = TreeNode(children=[inter, leaf3], distances=[3.0, 10.0])
    tree = Tree(root)
    
    assert abs(tree.get_distance(0, 1) - 12.0) < 0.001
    assert abs(tree.get_distance(0, 2) - 18.0) < 0.001
    assert abs(tree.get_distance(1, 2) - 20.0) < 0.001

def test_upgma():
    """Test UPGMA algorithm"""
    distances = np.array([
        [0.0, 1.0, 7.0, 7.0, 9.0],
        [1.0, 0.0, 7.0, 6.0, 8.0],
        [7.0, 7.0, 0.0, 2.0, 4.0],
        [7.0, 6.0, 2.0, 0.0, 3.0],
        [9.0, 8.0, 4.0, 3.0, 0.0],
    ], dtype=np.float64)
    
    tree = upgma(distances)
    newick = tree.to_newick(include_distance=False)
    # The exact structure might vary but should be a valid tree
    assert len(tree) == 5  # 5 leaves
    assert tree.root is not None

def test_neighbor_joining():
    """Test Neighbor Joining algorithm"""
    distances = np.array([
        [0.0, 1.0, 7.0, 7.0, 9.0],
        [1.0, 0.0, 7.0, 6.0, 8.0],
        [7.0, 7.0, 0.0, 2.0, 4.0],
        [7.0, 6.0, 2.0, 0.0, 3.0],
        [9.0, 8.0, 4.0, 3.0, 0.0],
    ], dtype=np.float64)
    
    tree = neighbor_joining(distances)
    newick = tree.to_newick(include_distance=False)
    # The exact structure might vary but should be a valid tree
    assert len(tree) == 5  # 5 leaves
    assert tree.root is not None

def main():
    start = time.time()
    
    test_distances()
    test_upgma()
    test_neighbor_joining()
    
    end = time.time()
    elapsed_ms = int((end - start) * 1000)
    
    print(f"codon       {elapsed_ms}ms")

main()