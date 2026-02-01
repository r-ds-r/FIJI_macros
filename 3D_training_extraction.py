import numpy as np
from pathlib import Path
from tifffile import imread, imwrite

train_dir = Path("C:/Users/rrodri42/OneDrive - McGill University/Rose/Microscopy/Analysis/Cellpose_Training/3D_cpsam-s100b/3D/3D_Ch3")
output_dir = Path("C:/Users/rrodri42/OneDrive - McGill University/Rose/Microscopy/Analysis/Cellpose_Training/3D_cpsam-s100b/3D/train")
output_dir.mkdir(exist_ok=True)

anisotropy = 3.85
min_masks = 1  # Minimum number of masks required to save a slice

# Find all image files (not masks or flows)
image_files = [f for f in train_dir.glob("*.tif") 
               if "_cp_masks" not in f.name 
               and "_flows" not in f.name
               and "_seg" not in f.name]

print(f"Found {len(image_files)} image files\n")

total_saved = 0

for img_file in image_files:
    print(f"Processing {img_file.name}")
    
    # Find corresponding seg file
    seg_file = img_file.with_name(img_file.stem + "_seg.npy")
    
    if not seg_file.exists():
        print(f"  WARNING: No _seg.npy file found, skipping")
        continue
    
    # Load image and masks
    img = imread(img_file)
    dat = np.load(seg_file, allow_pickle=True).item()
    
    if 'masks' not in dat:
        print(f"  WARNING: No masks in _seg.npy file, skipping")
        continue
    
    masks = dat['masks']
    
    print(f"  Image shape: {img.shape}, Masks shape: {masks.shape}")
    
    # Handle different dimensional formats
    if len(img.shape) == 4:
        if img.shape[-1] == 1:
            img = img[..., 0]
    
    nz, ny, nx = img.shape
    saved_this_file = 0
    
    # ===== XY slices (along Z axis) - ALL SLICES =====
    xy_saved = 0
    
    for z in range(nz):  # ALL Z slices
        mask_slice = masks[z, :, :]
        num_masks = len(np.unique(mask_slice)) - 1  # -1 to exclude background (0)
        
        if num_masks >= min_masks:
            img_slice = img[z, :, :]
            
            img_out_name = f"{img_file.stem}_XY_Z{z:03d}.tif"
            mask_out_name = f"{img_file.stem}_XY_Z{z:03d}_cp_masks.tif"
            
            imwrite(output_dir / img_out_name, img_slice)
            imwrite(output_dir / mask_out_name, mask_slice.astype(np.uint16))
            xy_saved += 1
            saved_this_file += 1
    
    print(f"  XY: saved {xy_saved}/{nz} slices")
    
    # ===== XZ slices (along Y axis) - ALL SLICES =====
    xz_saved = 0
    
    for y in range(ny):  # ALL Y slices
        mask_slice = masks[:, y, :]
        num_masks = len(np.unique(mask_slice)) - 1
        
        if num_masks >= min_masks:
            img_slice = img[:, y, :]
            
            img_out_name = f"{img_file.stem}_XZ_Y{y:03d}.tif"
            mask_out_name = f"{img_file.stem}_XZ_Y{y:03d}_cp_masks.tif"
            
            imwrite(output_dir / img_out_name, img_slice)
            imwrite(output_dir / mask_out_name, mask_slice.astype(np.uint16))
            xz_saved += 1
            saved_this_file += 1
    
    print(f"  XZ: saved {xz_saved}/{ny} slices")
    
    # ===== YZ slices (along X axis) - ALL SLICES =====
    yz_saved = 0
    
    for x in range(nx):  # ALL X slices
        mask_slice = masks[:, :, x]
        num_masks = len(np.unique(mask_slice)) - 1
        
        if num_masks >= min_masks:
            img_slice = img[:, :, x]
            
            img_out_name = f"{img_file.stem}_YZ_X{x:03d}.tif"
            mask_out_name = f"{img_file.stem}_YZ_X{x:03d}_cp_masks.tif"
            
            imwrite(output_dir / img_out_name, img_slice)
            imwrite(output_dir / mask_out_name, mask_slice.astype(np.uint16))
            yz_saved += 1
            saved_this_file += 1
    
    print(f"  YZ: saved {yz_saved}/{nx} slices")
    print(f"  Total for this file: {saved_this_file} slices\n")
    total_saved += saved_this_file

print(f"{'='*50}")
print(f"Done! Check {output_dir}")

# Count all slices
xy_count = len(list(output_dir.glob("*_XY_*_cp_masks.tif")))
xz_count = len(list(output_dir.glob("*_XZ_*_cp_masks.tif")))
yz_count = len(list(output_dir.glob("*_YZ_*_cp_masks.tif")))
total = xy_count + xz_count + yz_count

print(f"Created:")
print(f"  {xy_count} XY slice pairs")
print(f"  {xz_count} XZ slice pairs") 
print(f"  {yz_count} YZ slice pairs")
print(f"  {total} total training pairs")
print(f"{'='*50}")