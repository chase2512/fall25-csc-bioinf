# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.sequence.phylo"
__author__ = "Patrick Kunzmann, Tom David Müller"
__all__ = ["Tree", "TreeNode", "as_binary", "TreeError"]

import copy
from typing import List, Tuple, Optional


class TreeError(Static[Exception]):
    """
    An exception that occurs in context of tree topology.
    """
    pass


# Add hash support for set (Codon doesn't have frozenset)
@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h


class TreeNode:
    _index: int
    _distance: float
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: List[TreeNode]
    
    def __init__(self, children=None, distances=None, index: Optional[int] = None):
        if children is None:
            self._children = []
        else:
            self._children = [i for i in children]
        
        if index is not None:
            self._index = index
        else:
            self._index = -1
        
        self._distance = 0.0
        self._parent = None
        self._is_root = False
        
        if distances is not None:
            if len(self._children) != len(distances):
                raise ValueError(
                    f"Expected {len(self._children)} distances, "
                    f"got {len(distances)}"
                )
            for child, distance in zip(self._children, distances):
                child._parent = self
                child._distance = distance
    
    def _set_parent(self, parent: TreeNode, distance: float):
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance
    
    def copy(self) -> TreeNode:
        """
        Create a deep copy of this TreeNode.
        """
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            children_copy: List[TreeNode] = []
            distances_copy: List[float] = []
            for child in self._children:
                children_copy.append(child.copy())
                distances_copy.append(child._distance)
            return TreeNode(children_copy, distances_copy)

    @property
    def index(self) -> Optional[int]:
        return None if self._index == -1 else self._index
    
    @property
    def children(self) -> Optional[List[TreeNode]]:
        return None if len(self._children) == 0 else self._children
    
    @property
    def parent(self) -> Optional[TreeNode]:
        return self._parent
    
    @property
    def distance(self) -> Optional[float]:
        return None if self._parent is None else self._distance

    def is_leaf(self) -> bool:
        """
        Check if the node is a leaf node.
        """
        return len(self._children) == 0
    
    def is_root(self) -> bool:
        """
        Check if the node is a root node.
        """
        return self._is_root
    
    def as_root(self):
        """
        Finalize this node as a root node.
        """
        if self._parent is not None:
            raise TreeError("Node with parent cannot be root")
        self._is_root = True
        self._distance = 0.0
    
    def distance_to(self, node: TreeNode, topological: bool = False) -> float:
        """
        Get the distance from this node to another node.
        """
        lca = self.lowest_common_ancestor(node)
        
        if topological:
            # Count edges
            dist = 0.0
            current = self
            while current is not lca:
                dist += 1.0
                current = current._parent
            current = node
            while current is not lca:
                dist += 1.0
                current = current._parent
            return dist
        else:
            # Sum distances
            dist = 0.0
            current = self
            while current is not lca:
                dist += current._distance
                current = current._parent
            current = node
            while current is not lca:
                dist += current._distance
                current = current._parent
            return dist
    
    def lowest_common_ancestor(self, node: TreeNode) -> TreeNode:
        """
        Find the lowest common ancestor of this node and another node.
        """
        path1 = _create_path_to_root(self)
        path2 = _create_path_to_root(node)
        
        # Find first common node
        for node1 in path1:
            for node2 in path2:
                if node1 is node2:
                    return node1
        
        raise TreeError("Nodes do not share a common ancestor")
    
    def get_indices(self) -> List[int]:
        """
        Get all leaf indices under this node.
        """
        indices: List[int] = []
        _get_indices(self, indices)
        return indices

    def get_leaves(self) -> List[TreeNode]:
        """
        Get all leaf nodes under this node.
        """
        leaf_list: List[TreeNode] = []
        _get_leaves(self, leaf_list)
        return leaf_list
    
    def get_leaf_count(self) -> int:
        """
        Get the number of leaf nodes under this node.
        """
        return _get_leaf_count(self)
    
    def to_newick(self, labels=None, include_distance: bool = True, 
                  round_distance: Optional[int] = None) -> str:
        """
        Convert this node to Newick notation.
        """
        if self.is_leaf():
            if labels is not None:
                label = str(labels[self._index])
            else:
                label = str(self._index)
        else:
            # Intermediate node
            child_strs: List[str] = []
            for child in self._children:
                child_strs.append(child.to_newick(labels, include_distance, round_distance))
            label = f"({','.join(child_strs)})"
        
        if include_distance and self._parent is not None:
            if round_distance is not None:
                dist_str = str(round(self._distance, round_distance))
            else:
                dist_str = str(self._distance)
            return f"{label}:{dist_str}"
        else:
            return label
    
    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None) -> Tuple[TreeNode, float]:
        """
        Create a TreeNode from Newick notation.
        Returns (node, distance_to_parent).
        """
        newick = newick.strip()
        if len(newick) == 0:
            raise ValueError("Empty Newick string")
        
        # Parse distance if present
        distance = 0.0
        if ':' in newick:
            parts = newick.rsplit(':', 1)
            newick = parts[0]
            distance = float(parts[1])
        
        if newick[0] == '(':
            # Intermediate node
            if newick[-1] != ')':
                raise ValueError("Mismatched parentheses")
            
            # Find matching parentheses and split children
            children_str = newick[1:-1]
            children: List[TreeNode] = []
            distances: List[float] = []
            
            # Simple parser for comma-separated children
            depth = 0
            start = 0
            for i in range(len(children_str)):
                if children_str[i] == '(':
                    depth += 1
                elif children_str[i] == ')':
                    depth -= 1
                elif children_str[i] == ',' and depth == 0:
                    child, child_dist = TreeNode.from_newick(children_str[start:i], labels)
                    children.append(child)
                    distances.append(child_dist)
                    start = i + 1
            
            # Last child
            if start < len(children_str):
                child, child_dist = TreeNode.from_newick(children_str[start:], labels)
                children.append(child)
                distances.append(child_dist)
            
            node = TreeNode(children, distances)
            return node, distance
        else:
            # Leaf node
            if labels is not None:
                # Find index in labels
                idx = -1
                for i in range(len(labels)):
                    if labels[i] == newick:
                        idx = i
                        break
                if idx == -1:
                    raise ValueError(f"Label '{newick}' not found in labels list")
                return TreeNode(index=idx), distance
            else:
                # Parse as integer index
                return TreeNode(index=int(newick)), distance
    
    def __str__(self) -> str:
        return self.to_newick()
    
    def __eq__(self, item) -> bool:
        if not isinstance(item, TreeNode):
            return False
        
        if self._index != item._index:
            return False
        
        if len(self._children) == 0 and len(item._children) == 0:
            return True
        
        if len(self._children) == 0 or len(item._children) == 0:
            return False
        
        if len(self._children) != len(item._children):
            return False
        
        for i in range(len(self._children)):
            if self._children[i] != item._children[i]:
                return False
        
        return True
    
    def __hash__(self) -> int:
        if self.is_leaf():
            return hash(self._index)
        else:
            # Use set instead of tuple for hashing
            return hash(set(self._children))


def _get_leaves(node: TreeNode, leaf_list: List[TreeNode]):
    """Helper function to collect all leaf nodes."""
    if node._index == -1:
        # Intermediate node
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        # Leaf node
        leaf_list.append(node)


def _get_indices(node: TreeNode, indices: List[int]):
    """Helper function to collect all leaf indices."""
    if node._index == -1:
        # Intermediate node
        for child in node._children:
            _get_indices(child, indices)
    else:
        # Leaf node
        indices.append(node._index)


def _get_leaf_count(node: TreeNode) -> int:
    """Helper function to count leaf nodes."""
    if node._index == -1:
        # Intermediate node
        count = 0
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        # Leaf node
        return 1


def _create_path_to_root(node: TreeNode) -> List[TreeNode]:
    """
    Create a list of nodes representing the path from this node to the root.
    """
    path: List[TreeNode] = []
    current_node = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path


class Tree:
    """
    A Tree represents a rooted tree (e.g. alignment guide tree or phylogenetic tree).
    """
    
    _root: TreeNode
    _leaves: List[TreeNode]
    
    def __init__(self, root: TreeNode):
        root.as_root()
        self._root = root
        
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        
        # Create a temporary list that allows None
        temp_leaves: List[Optional[TreeNode]] = [None] * leaf_count
        
        for leaf in leaves_unsorted:
            idx = leaf.index
            if idx is None:
                raise ValueError("Leaf node has no index")
            if idx >= leaf_count:
                raise ValueError(f"Index {idx} is out of range for {leaf_count} leaves")
            if temp_leaves[idx] is not None:
                raise ValueError(f"Duplicate leaf index {idx}")
            temp_leaves[idx] = leaf
        
        # Now convert to non-optional list (all should be filled)
        self._leaves = []
        for leaf in temp_leaves:
            if leaf is None:
                raise ValueError("Missing leaf in tree")
            self._leaves.append(leaf)
    
    def __copy_create__(self):
        return Tree(self._root.copy())
    
    @property
    def root(self) -> TreeNode:
        return self._root
    
    @property
    def leaves(self) -> List[TreeNode]:
        return copy.copy(self._leaves)
    
    def get_distance(self, index1: int, index2: int, topological: bool = False) -> float:
        """
        Get the distance between two leaf nodes.
        """
        return self._leaves[index1].distance_to(
            self._leaves[index2], topological
        )
    
    def to_newick(self, labels=None, include_distance: bool = True, 
                  round_distance: Optional[int] = None) -> str:
        """
        Obtain the Newick notation of the tree.
        """
        return self._root.to_newick(
            labels, include_distance, round_distance
        ) + ";"
    
    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None) -> Tree:
        """
        Create a tree from a Newick notation.
        """
        newick = newick.strip()
        if len(newick) == 0:
            raise ValueError("Newick string is empty")
        # Remove terminal semicolon
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __str__(self) -> str:
        return self.to_newick()
    
    def __len__(self) -> int:
        return len(self._leaves)
    
    def __eq__(self, item) -> bool:
        if not isinstance(item, Tree):
            return False
        return self._root == item._root
    
    def __hash__(self) -> int:
        return hash(self._root)


def as_binary(tree_or_node):
    """
    Convert a tree into a binary tree.
    """
    if isinstance(tree_or_node, Tree):
        node = _as_binary(tree_or_node._root)
        return Tree(node)
    elif isinstance(tree_or_node, TreeNode):
        return _as_binary(tree_or_node)
    else:
        raise TypeError("Argument must be Tree or TreeNode")


def _as_binary(node: TreeNode) -> TreeNode:
    """
    The actual logic for converting to binary tree.
    Returns the converted node.
    """
    children = node.children
    if children is None:
        # Leaf node
        return TreeNode(index=node.index)
    elif len(children) == 1:
        # Single child - skip this node
        child_binary = _as_binary(children[0])
        # Add this node's distance to child's distance
        new_dist = node.distance if node.distance is not None else 0.0
        new_dist += children[0].distance
        new_node = child_binary.copy()
        if new_node._parent is not None:
            new_node._distance = new_dist
        return new_node
    elif len(children) > 2:
        # Multiple children - create binary divisions
        # Convert first two children
        child0_binary = _as_binary(children[0])
        child1_binary = _as_binary(children[1])
        
        current_div_node = TreeNode(
            [child0_binary, child1_binary],
            [children[0].distance, children[1].distance]
        )
        
        # Add remaining children
        for i in range(2, len(children)):
            child_binary = _as_binary(children[i])
            current_div_node = TreeNode(
                [current_div_node, child_binary],
                [0.0, children[i].distance]
            )
        
        return current_div_node
    else:
        # Exactly 2 children
        child0_binary = _as_binary(children[0])
        child1_binary = _as_binary(children[1])
        return TreeNode(
            [child0_binary, child1_binary],
            [children[0].distance, children[1].distance]
        )
