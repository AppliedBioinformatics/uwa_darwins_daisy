from ete3 import Tree, TreeStyle, faces, AttrFace
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image
from scripts.sgsgeneloss.popcolors import pop_colors

# Load trees.
tree_folder = Path("../../data/sgsgeneloss")
mash_tree = Tree(str(tree_folder / "mashtree_pav.nwk"))
nj_tree = Tree(str(tree_folder / "nj_pav_tree.nwk"))

# Create tanglegram layout.
def layout(node):
    # Add node names to leaves for visualisation.
    if node.is_leaf():
        name = node.name
        color = pop_colors.get(name, "black")
        face = AttrFace("name", fgcolor=color)
        faces.add_face_to_node(face, node, column=0)

# Create treestyle for both trees.
ts1 = TreeStyle()
ts1.show_leaf_name = False
ts1.layout_fn = layout
ts1.title.add_face(faces.TextFace("Mash Tree", fsize=20), column=0)

ts2 = TreeStyle()
ts2.show_leaf_name = False
ts2.layout_fn = layout
ts2.title.add_face(faces.TextFace("Neighbor Joining Tree", fsize=20), column=0)

# Render both trees side by side as images
mash_png = "../../data/sgsgeneloss_/mash_tree.png"
nj_png = "../../data/sgsgeneloss_/nj_tree.png"
mash_tree.render(mash_png, tree_style=ts1, w=500)
nj_tree.render(nj_png, tree_style=ts2, w=500)

# Load PNGs
img1 = Image.open(mash_png)
img2 = Image.open(nj_png)

# Plot side-by-side with matplotlib
fig, axs = plt.subplots(1, 2, figsize=(14, 8))
axs[0].imshow(img1)
axs[0].axis('off')
axs[1].imshow(img2)
axs[1].axis('off')
plt.tight_layout()
plt.savefig("../../data/sgsgeneloss_/mash_v_nj_pav.png")
