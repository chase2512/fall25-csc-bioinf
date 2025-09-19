import copy
from python import matplotlib.pyplot as plt
from typing import Optional


def reverse_complement(key: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complement[key[i]]
    return ''.join(key)


class Node:
    _children: set[int]
    _count: int
    kmer: str
    visited: bool
    depth: int
    max_depth_child: Optional[int]

    def __init__(self, kmer: str):
        self._children = set[int]()
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def add_child(self, kmer: int):
        self._children.add(kmer)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> list[int]:
        return list(self._children)

    def remove_children(self, target: set[int]):
        self._children = self._children - target


class DBG:
    k: int
    nodes: dict[int, Node]
    kmer2idx: dict[str, int]
    kmer_count: int

    def __init__(self, k: int, data_list: list[list[str]]):
        self.k = k
        self.nodes = {}
        self.kmer2idx = {}
        self.kmer_count = 0
        # build
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: list[list[str]]):
        # check data list
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list: list[list[str]]):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                    self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def show_count_distribution(self):
        count = [0] * 30
        for idx in self.nodes:
            count[self.nodes[idx].get_count()] += 1
        print(count[0:10])
        # Uncomment below to use matplotlib plotting
        # plt.plot(count)
        # plt.show()

    def _add_node(self, kmer: str) -> int:
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child: int) -> int:
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx: int) -> list[int]:
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth, max_child = 0, None
            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth, max_child = depth, child
            self.nodes[idx].depth, self.nodes[idx].max_depth_child = max_depth + 1, max_child
        return self.nodes[idx].depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self) -> list[int]:
        max_depth, max_idx = 0, None
        for idx in self.nodes.keys():
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth, max_idx = depth, idx

        path: list[int] = []
        while max_idx is not None:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: list[int]):
        for idx in path:
            del self.nodes[idx]
        path_set: set[int] = set(path)
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path: list[int]) -> Optional[str]:
        if len(path) < 1:
            return None
        concat = copy.copy(self.nodes[path[0]].kmer)
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self) -> Optional[str]:
        # reset params in nodes for getting longest path
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig