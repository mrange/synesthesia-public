# Scene Review Process

This document outlines the complete process for reviewing and improving synesthesia scenes, including controls, descriptions, tags, and visual presentation.

## Overview

When reviewing a scene, we want to ensure that:
1. Control descriptions accurately match what they actually do in the shader
2. Control descriptions are succinct and clear
3. Control names are descriptive and user-friendly
4. Control names match between scene.json and the shader code
5. Scene description is engaging and captures the visual/interactive essence
6. Tags are relevant and improve discoverability
7. Scene thumbnail accurately represents the scene

## Step-by-Step Process

### 1. Read All Scene Files

Read the following files for the scene being reviewed:
- `scene.json` - Contains control definitions, scene description, and tags
- `main.glsl` (or shader file) - Contains the actual implementation
- `script.js` - May contain additional logic (if relevant)
- `scene.png` - Visual thumbnail of the scene

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
- Keep descriptions succinct (typically 7-10 words)
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

### 7. Review and Update Scene Description

Look at the scene.png thumbnail and understand the visual aesthetic:
- What is the dominant visual style? (retro, cyberpunk, abstract, etc.)
- What are the main interactive elements? (audio-reactive, reflective, animated, etc.)
- What emotional response does it evoke?

Write or update the DESCRIPTION field in scene.json:
- Keep it engaging and evocative (2-4 sentences)
- Mention key visual elements visible in the thumbnail
- Describe the interactive/audio-reactive aspects
- Convey the overall aesthetic and mood
- Format: Narrative description, then attributions

**Example:**
"A retro-futuristic satellite orbits a neon dreamscape, its antenna responding to the audio beat. A glowing pink triangle floats above an infinite grid of electric blue, reflecting the cosmic sunset beyond. The entire scene pulses with 80s synthwave energy. Created for SceneSat."

### 8. Update Scene Tags

Add relevant tags to improve discoverability. Consider:
- **Aesthetic:** retro, synthwave, cyberpunk, abstract, geometric, organic, etc.
- **Interactive:** audio-reactive, animated, responsive, interactive, etc.
- **Visual effects:** crt, neon, reflection, grid, distortion, glow, etc.
- **Subjects:** satellite, antenna, tunnel, landscape, abstract, etc.
- **Creator name(s)**

**Guidelines:**
- Include 8-12 relevant tags
- Use lowercase with hyphens for multi-word tags (e.g., "audio-reactive")
- Include creator name as first tag (e.g., "mrange")
- Avoid redundant tags

**Example tags:** `["mrange", "neon", "satellite", "antenna", "audio-reactive", "crt", "retro", "synthwave", "reflection", "atmospheric"]`

### 9. Verify Scene Thumbnail

Check the scene.png thumbnail:
- Does it represent the scene described in the DESCRIPTION field?
- Are key visual features visible?
- Is the author credit appropriately sized (not obscuring the scene)?
- Does it accurately show the default control state?

Note: Some scenes may intentionally use a title card style thumbnail, which is acceptable.

### 10. Attribution and Credits

If using external assets (fonts, techniques, etc.), add proper attribution:
- **Font credits:** Format as "Font: [Name] by [Artist] ([URL]), licensed as [License]"
- Include links to artist portfolios or licensing information
- Maintain consistency with existing attributions in the project

## Example Changes

### Control Description Improvements

**Before:**
```json
"DESCRIPTION": "Spectrum threshold – filter out low frequency audio"
```

**After:**
```json
"DESCRIPTION": "Audio threshold – filter out quiet audio values"
```

### Control Name Improvements

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

### Scene Description Improvement

**Before:**
```json
"DESCRIPTION": "Neon Scenesat"
```

**After:**
```json
"DESCRIPTION": "A retro-futuristic satellite orbits a neon dreamscape, its antenna responding to the audio beat. A glowing pink triangle floats above an infinite grid of electric blue, reflecting the cosmic sunset beyond. Media content bleeds through as the satellite rotates, while CRT scan lines flicker across the screen. The entire scene pulses with 80s synthwave energy. Created for SceneSat (https://www.scenesat.com/). Font: Cyber Brush by burhanafif (https://hanscostudio.com/), licensed as Non-Commercial."
```

## Checklist

- [ ] Read scene.json, main.glsl, script.js, and scene.png
- [ ] Verify each control description matches implementation
- [ ] Check for missing shader constant declarations
- [ ] Update control descriptions for accuracy and clarity
- [ ] Review control names for clarity
- [ ] Update control names in scene.json
- [ ] Update variable names in shader to match
- [ ] Verify all occurrences updated in shader
- [ ] Review scene thumbnail and visual elements
- [ ] Write engaging scene description (2-4 sentences)
- [ ] Add relevant tags (8-12 tags)
- [ ] Add/verify proper attribution and credits
- [ ] Check grammar and language in descriptions
- [ ] Test that no shader compilation errors exist

## Notes

- Keep descriptions succinct but informative
- Prioritize user-friendly language over technical accuracy
- Always update both scene.json AND shader code when renaming
- Use consistent naming conventions across similar controls
- Document any special cases (e.g., negative values, ranges)
- Scene descriptions should be evocative and capture the essence of what you see in the thumbnail
- Tags should improve discoverability—include aesthetic, interactive, and visual effect tags
- Proper attribution respects artist contributions and follows community best practices
