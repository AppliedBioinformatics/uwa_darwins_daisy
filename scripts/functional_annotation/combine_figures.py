from PIL import Image
import math
import os

# === CONFIG ===
input_folder = "../../plots/functional_annotation"       # folder with your PNGs
output_path = "../../plots/functional_annotation/figure_panel.png"  # output file
grid_cols = 4                   # e.g. 4 columns → 2 rows (7 figs)
padding = 20                    # space between images

# === LOAD IMAGES ===
image_files = sorted([
    os.path.join(input_folder, f)
    for f in os.listdir(input_folder)
    if f.endswith(".png")
])

images = [Image.open(f) for f in image_files]
img_width, img_height = images[0].size

# === CALCULATE GRID SIZE ===
n_images = len(images)
grid_rows = math.ceil(n_images / grid_cols)

panel_width = grid_cols * img_width + (grid_cols - 1) * padding
panel_height = grid_rows * img_height + (grid_rows - 1) * padding

panel = Image.new("RGB", (panel_width, panel_height), color=(255, 255, 255))

# === PASTE IMAGES ===
for idx, img in enumerate(images):
    row = idx // grid_cols
    col = idx % grid_cols
    x = col * (img_width + padding)
    y = row * (img_height + padding)
    panel.paste(img, (x, y))

# === SAVE COMBINED PANEL ===
panel.save(output_path)
print(f"Saved panel to: {output_path}")