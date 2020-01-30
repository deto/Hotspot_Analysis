import random
import ete3
from ete3 import Tree
import os

random.seed(1028)

tree_file = snakemake.input['tree']
tree_out = snakemake.output['tree']

os.makedirs(os.path.dirname(tree_out), exist_ok=True)

t = Tree(tree_file, format=1)

leaf_order = []
for tn in t.traverse('postorder'):
    if tn.is_leaf():
        leaf_order.append(tn.name)


leaf_order_shuff = leaf_order.copy()
random.shuffle(leaf_order_shuff)

leaf_map = {a: b for a, b in zip(leaf_order, leaf_order_shuff)}

for tn in t.traverse('postorder'):
    if tn.is_leaf():
        old_name = tn.name
        new_name = leaf_map[old_name]
        tn.name = new_name

t.write(format=1, outfile=tree_out)
