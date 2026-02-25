# Spatial Proteomics Portfolio Website

## Context

This repo (`skelly001/3D_Islet_Mapping`) contains a 6-stage spatial proteomics pipeline that maps 2,215 proteins at cellular resolution in human pancreatic islet tissue. The user wants a visually striking portfolio website on GitHub Pages to showcase this work during their job search. The current `plan.md` in the repo is a starting sketch; this plan supersedes it with a much richer design.

**Key visual assets already in the repo:**
- 2,215 protein map PNGs (~140KB each, 318MB total) — navy-to-yellow heatmaps on cell polygons
- Cell type assignment map (acinar/alpha/beta), segmented cells, pixel assignment, DAPI mask
- H&E stain JPEG (526KB, web-ready), islet + acinar marker 3-panels
- 4 TIFF tissue images (~48MB each) — need conversion to web format

**Critical spatial alignment fact:** The cell type assignment map and ALL protein maps share identical cell polygon geometries (same `sf` objects from the R pipeline), so they overlay perfectly for slider comparisons.

**Dimension mismatch note:** The H&E image (`Islet23_Slide60_10X_H&E.jpg`) is from a different tissue section (Slide 60) than the protein maps (Slide 61), and likely has different image dimensions/field of view. The H&E should be displayed as standalone tissue context (not in a slider with protein maps). Only images sharing the same `sf` polygon geometry should be paired in sliders.

---

## Files to Create/Modify

| File | Action | Purpose |
|------|--------|---------|
| `index.html` | Create | Single-page website (all CSS/JS inline, ~40KB + embedded 38KB JSON manifest) |
| `convert_tiffs.py` | Create | Python/Pillow script to convert 48MB TIFFs → web JPEGs |
| `website/converted/*.jpg` | Generate | 4 converted tissue images (~200-500KB each) |
| `.github/workflows/pages.yml` | Create | GitHub Actions Pages deployment |
| `.gitignore` | Modify | Comment out `docs/` line, add clarifying comments |
| `README.md` | Create | Repo description + live website link |
| `plan.md` | Overwrite | Replace with this plan (keep in repo for future reference) |

---

## Website Design

### Aesthetic Direction: "Dark Scientific Editorial"
- **Background:** Near-black with blue undertone (`#0d0d1a`) — makes fluorescence images and protein maps POP
- **Accents:** Cyan (`#00d4ff`) from DAPI, yellow (`#ffd700`) from protein maps
- **Typography:** System font stack — geometric sans headings (Avenir/Segoe UI), refined serif body (Charter/Georgia)
- **Animations:** Scroll-triggered fade-in reveals, smooth slider interactions
- **Self-contained:** Zero external dependencies, all inline

### Section Layout

**1. Hero (full viewport)**
- Dark background with subtle CSS radial gradient (faint cyan glow, mimicking DAPI scatter)
- Title: "Pseudo-Single-Cell Spatial Proteomics of Human Pancreatic Islets"
- Subtitle: "2,215 protein abundance maps at cellular resolution"
- Scroll-down chevron with CSS bounce animation

**2. Project Overview (two-column)**
- Left (60%): Brief research description — tissue type, technique (nanopots MS + IF), outcome
- Right (40%): H&E stain image (`data/images/Islet23_Slide60_10X_H&E.jpg`) with caption — shown as standalone tissue context (different slide/dimensions than protein maps, NOT used in slider overlays)
- Below: 3 stat callout boxes — "2,215 Proteins", "3 Cell Types", "6-Stage Pipeline"

**3. Interactive Slider Comparisons (full-width)**
- Header: "From Tissue to Protein Maps"
- **Slider 1:** Cell Type Assignment ↔ Insulin map (`INS_P01308.png`) — perfect polygon alignment
- **Slider 2:** Cell Type Assignment ↔ Glucagon map (`GCG_P01275.png`) — perfect polygon alignment
- **Slider 3:** DAPI fluorescence ↔ Cell boundaries overlay (converted TIFFs)
- Each slider: draggable divider, CSS `clip-path` approach, touch/pointer support, left/right labels
- Images pre-cropped (via `convert_tiffs.py`) to remove ggplot legends for clean overlay

**4. Analysis Pipeline Visualization**
- Horizontal CSS flexbox flowchart with 6 stage cards connected by arrows
- Each card: stage number badge, title, tool/method, representative thumbnail
  1. 3D Spatial Proteomics — Giotto (thumbnail: `RD0/spatPlot2D.png`)
  2. Cell Segmentation — CellPose 2.0 (thumbnail: `DAPI_mask_for_fig_5.png`)
  3. Cell Type Assignment — IF RGB classification (thumbnail: `cell_type_assignment.png`)
  4. Pixel-to-Cell Mapping — ImageJ macros (thumbnail: `Pixel Assignment.png`)
  5. MS Integration + RNA-Seq — Azimuth/Seurat (thumbnail: `islet_protein_markers.png`)
  6. Protein Map Generation — cell-type adjusted (thumbnail: `INS_P01308.png`)
- Stacks vertically on mobile

**5. Marker Protein Showcase (full-bleed)**
- Islet markers 3-panel (Secretogranin-2, Insulin, Glucagon)
- Acinar markers 3-panel (REG1A, REG1B, CPA1)
- Dark framed with biology captions

**6. Searchable Protein Map Gallery (main interactive feature)**
- Sticky search bar: live filtering by gene name or UniProt ID (debounced 150ms)
- Result count badge + A-Z sort
- CSS Grid (`auto-fill, minmax(220px, 1fr)`) of thumbnail cards
- Each card: lazy-loaded image, gene name, UniProt ID, hover glow effect
- **Load strategy:** Show 48 at a time + "Load More" button
- **Lightbox modal:** Full-size image, gene name header, UniProt link, arrow/keyboard navigation, ESC close
- Navy placeholder backgrounds on cards before images load

**7. Methods Summary**
- Brief technical description, link to GitHub repo, manuscript reference

**8. Footer**
- Researcher info, GitHub/contact links, "Built with vanilla HTML/CSS/JS"

---

## TIFF Conversion Script (`convert_tiffs.py`)

Converts 4 TIFF files (48MB each) to web-optimized JPEGs. Also crops slider comparison images to remove ggplot legends.

**Dependencies:** `Pillow` (+ `numpy` for 16-bit TIFFs)

**Conversions:**

| Source (in `data/1-ImageJ_cell_ROI/`) | Output (in `website/converted/`) | Notes |
|---|---|---|
| `Slide61_Max_Projection_DAPI.tif` | `DAPI_fluorescence.jpg` | JPEG q85, max 2000px edge |
| `Result of Slide61_Max_Projection_INS_GCG_DAPI.tif` | `INS_GCG_DAPI_composite.jpg` | JPEG q85, max 2000px edge |
| `cropped.tif` | `tissue_cropped.jpg` | JPEG q85, max 2000px edge |
| `Slide61-islet23_with_cell_boundaries_with_edges.tif` | `cell_boundaries_overlay.jpg` | JPEG q85, max 2000px edge |

**Additional crops for slider:**
- Crop `cell_type_assignment.png` → `cell_type_cropped.jpg` (remove right legend, ~72% width)
- Crop `INS_P01308.png` → `INS_cropped.jpg` (same crop for alignment)
- Crop `GCG_P01275.png` → `GCG_cropped.jpg` (same crop for alignment)

---

## Performance Strategy (318MB of protein maps)

1. Gallery shows 48 images at a time + "Load More" (6.7MB per batch)
2. `loading="lazy"` + Intersection Observer with 200px margin
3. Explicit `width`/`height` on all `<img>` tags to prevent layout shift
4. Navy CSS placeholder backgrounds on cards
5. Search reduces image load count (search "INS" → only matching images load)
6. `decoding="async"` on gallery images
7. Hero/pipeline sections load first; gallery images only load on scroll
8. GitHub Pages Fastly CDN provides caching after first load
9. Start without thumbnails (140KB per image is acceptable); add as Phase 2 if needed

---

## GitHub Pages Deployment

**`.github/workflows/pages.yml`:**
- Trigger: push to `main`
- Uses: `actions/checkout@v4` → `actions/configure-pages@v5` → `actions/upload-pages-artifact@v3` (path: `.`) → `actions/deploy-pages@v4`
- Deploys entire repo root — `index.html` at root becomes homepage

**User must:** Enable Pages in repo Settings → Source: GitHub Actions

**Live URL:** `https://skelly001.github.io/3D_Islet_Mapping/`

---

## `.gitignore` Changes

- Comment out `docs/` line (allow future use)
- Add comment block documenting website files
- No other changes needed — `*.jpg` and `*.png` are not globally ignored, `website/` is not ignored

---

## Execution Order

### Phase 0: Replace plan.md
1. Overwrite `plan.md` in the repo root with this plan content (replaces the old sketch)

### Phase 1: Image Preparation
2. Create `convert_tiffs.py`
3. Run it to generate `website/converted/*.jpg` (4 tissue images + 3 slider crops)
4. Verify outputs look correct

### Phase 2: Build Website
5. Generate 2,215-entry protein JSON array by listing `output/RD5-final_protein_maps/final_protein_maps/`
6. Create `index.html` with all 8 sections, inline CSS/JS, embedded JSON, slider component, gallery, lightbox

### Phase 3: Deployment Config
7. Create `.github/workflows/pages.yml`
8. Update `.gitignore`
9. Create `README.md`

### Phase 4: Test + Stage
10. Test locally with `python -m http.server` from repo root
11. Verify: slider interaction, search filtering, lazy loading, lightbox, mobile responsive
12. Stage all new/modified files (do NOT commit — leave for user review)

---

## Verification Plan

1. Open `index.html` via local HTTP server (`python -m http.server`)
2. Check hero section renders with dark theme and correct typography
3. Verify H&E image loads in overview section
4. Test all 3 slider overlays — drag divider, check alignment, test on touch device
5. Verify pipeline flowchart displays correctly on desktop and mobile
6. Test gallery search: type "INS" → should show insulin-related proteins
7. Test "Load More" button: click → 48 more images appear
8. Test lightbox: click image → modal opens, arrow keys navigate, ESC closes
9. Check lazy loading: open Network tab, scroll slowly, verify images load on demand
10. Test mobile responsive: resize to 375px width, verify all sections stack correctly
11. After push: verify `https://skelly001.github.io/3D_Islet_Mapping/` is live
12. Confirm `plan.md` is in the repo for future reference
