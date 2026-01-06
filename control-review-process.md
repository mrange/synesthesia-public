# Scene Control Review Process

This document outlines the process for reviewing and improving control descriptions and names in synesthesia scenes.

## Overview

When reviewing a scene, we want to ensure that:
1. Control descriptions accurately match what they actually do in the shader
2. Control descriptions are succinct and clear
3. Control names are descriptive and user-friendly
4. Control names match between scene.json and the shader code

## Step-by-Step Process

### 1. Read All Scene Files

Read the following files for the scene being reviewed:
- `scene.json` - Contains control definitions with NAME and DESCRIPTION fields
- `main.glsl` (or shader file) - Contains the actual implementation
- `script.js` - May contain additional logic (if relevant)

### 2. Analyze Control Descriptions vs Implementation

For each control in scene.json:
- Find where it's used in the shader code (search for the control NAME)
- Understand what the control actually does in the code
- Compare the DESCRIPTION to the actual implementation
- Check if the description is accurate and succinct

**Common issues to look for:**
- Descriptions that are too technical (e.g., "FFT distinct" vs "audio threshold")
- Descriptions that don't match the actual behavior (e.g., "density" when it's actually a "threshold")
- Descriptions that are too verbose
- Descriptions missing important details (e.g., negative values invert behavior)

### 3. Update Control Descriptions

Edit scene.json to update DESCRIPTION fields where needed:
- Keep descriptions succinct (typically 5-10 words)
- Use clear, non-technical language when possible
- Format: "Control purpose – specific behavior details"
- Examples:
  - "Tower threshold – higher values reveal more neon buildings"
  - "Audio reactivity – building heights respond to spectrum"
  - "Media intensity – amplify blend (negative inverts)"

### 4. Review Control Names

Analyze each control NAME for clarity:
- Is the name self-descriptive?
- Could it be confused with other controls?
- Does it use clear terminology?

**Common improvements:**
- `height` → `camera_height` (more specific)
- `height_mul` → `audio_reactivity` (more descriptive)
- `fft_distinct` → `audio_threshold` (less technical)
- `media_transparency` → `media_opacity` (standard terminology)
- `path_control` → `path_wave` (describes what it controls)

### 5. Update Control Names in Both Files

When renaming controls, update BOTH files:

**In scene.json:**
- Update the NAME field for each control

**In the shader (main.glsl):**
- Update the constant declarations (typically in #ifdef KODELIFE block)
- Find and replace all usages of the old name with the new name
- Common places to check:
  - Function parameters and calculations
  - Conditional statements
  - Blending/mixing operations

**Important:** Use search/find to ensure all occurrences are updated.

### 6. Verify Consistency

After updates:
- Verify all control names match between scene.json and shader
- Ensure no old variable names remain in the shader
- Check that descriptions still make sense with new names

### 7. Review Scene Thumbnail (Optional)

Check the scene.png thumbnail:
- Does it represent the scene described in the DESCRIPTION field?
- Are key visual features visible?
- Is the author credit appropriately sized (not obscuring the scene)?

Note: Some scenes may intentionally use a title card style thumbnail, which is acceptable.

## Example Changes

### Description Improvements

**Before:**
```json
"DESCRIPTION": "Spectrum threshold – filter out low frequency audio"
```

**After:**
```json
"DESCRIPTION": "Audio threshold – filter out quiet audio values"
```

### Name Improvements

**Before (scene.json):**
```json
"NAME": "neon_towers"
```

**Before (main.glsl):**
```glsl
const float neon_towers = .2;
// Later in code:
x1 = HH.w < neon_towers && isr && ...
```

**After (scene.json):**
```json
"NAME": "tower_threshold"
```

**After (main.glsl):**
```glsl
const float tower_threshold = .2;
// Later in code:
x1 = HH.w < tower_threshold && isr && ...
```

## Checklist

- [ ] Read scene.json, main.glsl, and script.js
- [ ] Verify each control description matches implementation
- [ ] Update descriptions for accuracy and clarity
- [ ] Review control names for clarity
- [ ] Update control names in scene.json
- [ ] Update variable names in shader to match
- [ ] Verify all occurrences updated in shader
- [ ] Test that no compilation errors exist
- [ ] Review scene.png thumbnail (optional)

## Notes

- Keep descriptions succinct but informative
- Prioritize user-friendly language over technical accuracy
- Always update both scene.json AND shader code when renaming
- Use consistent naming conventions across similar controls
- Document any special cases (e.g., negative values, ranges)
