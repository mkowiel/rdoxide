//! Updating `MolDrawOptions` from JSON (`MolDrawOptionsJSONParsers.cpp` of
//! RDKit 2026.09.1). Keys that are absent, or whose values cannot be read,
//! leave the option as it was.
//!
//! Copyright (C) 2020-2022 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::helpers::{
    ColourPalette, DrawColour, LegendPosition, MolDrawOptions, MultiColourHighlightStyle,
    avalon_palette, bw_palette, cdk_palette, dark_mode_palette, default_palette, draw_element,
};
use crate::Result;
use crate::rdgeneral::ptree::PTree;

fn get_f64(pt: &PTree, key: &str, default: f64) -> f64 {
    pt.get_str(key)
        .and_then(|v| v.trim().parse().ok())
        .unwrap_or(default)
}

fn get_string(pt: &PTree, key: &str, default: &str) -> String {
    pt.get_str(key).unwrap_or(default).to_string()
}

/// Reads `[r, g, b]` or `[r, g, b, a]` (through `float` like upstream).
fn get_rgba(node: &PTree, colour: &mut DrawColour) {
    let vals: Vec<f64> = node
        .children()
        .filter_map(|c| c.value().trim().parse::<f32>().ok())
        .map(f64::from)
        .collect();
    if vals.len() >= 3 {
        colour.r = vals[0];
        colour.g = vals[1];
        colour.b = vals[2];
        if let Some(a) = vals.get(3) {
            colour.a = *a;
        }
    }
}

fn get_colour_option(pt: &PTree, key: &str, colour: &mut DrawColour) {
    if let Some(node) = pt.child(key) {
        get_rgba(node, colour);
    }
}

fn get_colour_palette_option(pt: &PTree, key: &str, palette: &mut ColourPalette) {
    let Some(node) = pt.child(key) else { return };
    if node.children().next().is_none() {
        // a terminal value names a preset palette
        let name = node.value().to_lowercase();
        let preset = match name.as_str() {
            "default" => default_palette(),
            "avalon" => avalon_palette(),
            "cdk" => cdk_palette(),
            "darkmode" => dark_mode_palette(),
            "bw" => bw_palette(),
            _ => return,
        };
        *palette = preset;
    } else {
        for (atomic_num, colour_node) in node.members() {
            if let Ok(num) = atomic_num.trim().parse::<i32>() {
                let mut colour = DrawColour::default();
                get_rgba(colour_node, &mut colour);
                palette.insert(num, colour);
            }
        }
    }
}

/// `flagsFromJson<DrawElement>`: the flags set to true and not to false, and
/// whether any of the names was found.
fn draw_element_flags(node: &PTree) -> Option<u32> {
    const NAMES: [(&str, u32); 9] = [
        ("NONE", draw_element::NONE),
        ("PRESHAPES", draw_element::PRESHAPES),
        ("BONDS", draw_element::BONDS),
        ("ATOMLABELS", draw_element::ATOMLABELS),
        ("HIGHLIGHTS", draw_element::HIGHLIGHTS),
        ("ANNOTATIONS", draw_element::ANNOTATIONS),
        ("RADICALS", draw_element::RADICALS),
        ("POSTSHAPES", draw_element::POSTSHAPES),
        ("ALL", draw_element::ALL),
    ];
    let mut by_value = [0u32; 2];
    let mut found = false;
    for (name, bits) in NAMES {
        if let Some(v) = node.child(name) {
            found = true;
            let on = matches!(v.value().trim(), "true" | "1");
            by_value[usize::from(on)] |= bits;
        }
    }
    found.then_some(by_value[1] & !by_value[0])
}

/// `updateMolDrawOptionsFromJSON`
pub fn update_mol_draw_options_from_json(opts: &mut MolDrawOptions, json: &str) -> Result<()> {
    if json.is_empty() {
        return Ok(());
    }
    let pt = PTree::from_json(json)?;

    macro_rules! bools {
        ($($key:literal => $field:ident),* $(,)?) => {
            $( opts.$field = pt.get_bool($key, opts.$field); )*
        };
    }
    macro_rules! ints {
        ($($key:literal => $field:ident),* $(,)?) => {
            $( opts.$field = pt.get_i32($key, opts.$field); )*
        };
    }
    macro_rules! doubles {
        ($($key:literal => $field:ident),* $(,)?) => {
            $( opts.$field = get_f64(&pt, $key, opts.$field); )*
        };
    }
    macro_rules! strings {
        ($($key:literal => $field:ident),* $(,)?) => {
            $( opts.$field = get_string(&pt, $key, &opts.$field); )*
        };
    }

    bools! {
        "atomLabelDeuteriumTritium" => atom_label_deuterium_tritium,
        "dummiesAreAttachments" => dummies_are_attachments,
        "circleAtoms" => circle_atoms,
        "splitBonds" => split_bonds,
        "continuousHighlight" => continuous_highlight,
        "fillHighlights" => fill_highlights,
        "includeAtomTags" => include_atom_tags,
        "clearBackground" => clear_background,
        "noAtomLabels" => no_atom_labels,
        "scaleBondWidth" => scale_bond_width,
        "scaleHighlightBondWidth" => scale_highlight_bond_width,
        "prepareMolsBeforeDrawing" => prepare_mols_before_drawing,
        "addAtomIndices" => add_atom_indices,
        "addBondIndices" => add_bond_indices,
        "isotopeLabels" => isotope_labels,
        "dummyIsotopeLabels" => dummy_isotope_labels,
        "addStereoAnnotation" => add_stereo_annotation,
        "showAllCIPCodes" => show_all_cip_codes,
        "atomHighlightsAreCircles" => atom_highlights_are_circles,
        "centreMoleculesBeforeDrawing" => centre_molecules_before_drawing,
        "explicitMethyl" => explicit_methyl,
        "includeMetadata" => include_metadata,
        "includeRadicals" => include_radicals,
        "comicMode" => comic_mode,
        "includeChiralFlagLabel" => include_chiral_flag_label,
        "simplifiedStereoGroupLabel" => simplified_stereo_group_label,
        "addStereoGroupAnnotation" => add_stereo_group_annotation,
        "unspecifiedStereoIsUnknown" => unspecified_stereo_is_unknown,
        "singleColourWedgeBonds" => single_colour_wedge_bonds,
        "singleColourBonds" => single_colour_bonds,
        "useMolBlockWedging" => use_mol_block_wedging,
        "drawMolsSameScale" => draw_mols_same_scale,
        "useComplexQueryAtomSymbols" => use_complex_query_atom_symbols,
        "bracketsAroundAtomLists" => brackets_around_atom_lists,
        "standardColoursForHighlightedAtoms" => standard_colours_for_highlighted_atoms,
    }
    doubles! {
        "highlightRadius" => highlight_radius,
        "legendFraction" => legend_fraction,
        "baseFontSize" => base_font_size,
        "annotationFontScale" => annotation_font_scale,
        "multipleBondOffset" => multiple_bond_offset,
        "padding" => padding,
        "componentPadding" => component_padding,
        "additionalAtomLabelPadding" => additional_atom_label_padding,
        "bondLineWidth" => bond_line_width,
        "fixedScale" => fixed_scale,
        "fixedBondLength" => fixed_bond_length,
        "rotate" => rotate,
        "variableAtomRadius" => variable_atom_radius,
        "scalingFactor" => scaling_factor,
    }
    ints! {
        "flagCloseContactsDist" => flag_close_contacts_dist,
        "legendFontSize" => legend_font_size,
        "maxFontSize" => max_font_size,
        "minFontSize" => min_font_size,
        "fixedFontSize" => fixed_font_size,
        "highlightBondWidthMultiplier" => highlight_bond_width_multiplier,
        "variableBondWidthMultiplier" => variable_bond_width_multiplier,
    }
    strings! {
        "fontFile" => font_file,
        "stereoGroupAndLabel" => stereo_group_and_label,
        "stereoGroupOrLabel" => stereo_group_or_label,
        "stereoGroupAbsLabel" => stereo_group_abs_label,
    }

    // this sets a whole family of colours at once, so it has to be applied
    // before the individual colour options, which are allowed to override it
    if let Some(scheme) = pt.get_str("colourScheme") {
        match scheme.to_lowercase().as_str() {
            "dark" | "darkmode" => opts.set_dark_mode(),
            "monochrome" => opts.set_monochrome_mode(
                DrawColour::new(0.0, 0.0, 0.0),
                DrawColour::new(1.0, 1.0, 1.0),
            ),
            _ => {}
        }
    }
    get_colour_option(&pt, "highlightColour", &mut opts.highlight_colour);
    get_colour_option(&pt, "backgroundColour", &mut opts.background_colour);
    get_colour_option(&pt, "queryColour", &mut opts.query_colour);
    get_colour_option(&pt, "legendColour", &mut opts.legend_colour);
    get_colour_option(&pt, "symbolColour", &mut opts.symbol_colour);
    get_colour_option(&pt, "annotationColour", &mut opts.annotation_colour);
    get_colour_option(&pt, "atomNoteColour", &mut opts.atom_note_colour);
    get_colour_option(&pt, "bondNoteColour", &mut opts.bond_note_colour);
    get_colour_option(
        &pt,
        "variableAttachmentColour",
        &mut opts.variable_attachment_colour,
    );
    get_colour_palette_option(&pt, "atomColourPalette", &mut opts.atom_colour_palette);
    if let Some(labels) = pt.child("atomLabels") {
        for (idx, label) in labels.members() {
            if let Ok(idx) = idx.trim().parse::<i32>() {
                opts.atom_labels.insert(idx, label.value().to_string());
            }
        }
    }
    if let Some(style) = pt.get_str("multiColourHighlightStyle") {
        match style {
            "Lasso" => opts.multi_colour_highlight_style = MultiColourHighlightStyle::Lasso,
            "CircleAndLine" => {
                opts.multi_colour_highlight_style = MultiColourHighlightStyle::CircleAndLine
            }
            _ => {}
        }
    }
    if let Some(pos) = pt.get_str("legendPosition") {
        opts.legend_position = match pos {
            "Top" => LegendPosition::Top,
            "Left" => LegendPosition::Left,
            "Right" => LegendPosition::Right,
            _ => LegendPosition::Bottom,
        };
    }
    opts.legend_vertical_text = pt.get_bool("legendVerticalText", opts.legend_vertical_text);
    if let Some(flags) = pt
        .child("drawingExtentsInclude")
        .and_then(draw_element_flags)
    {
        opts.drawing_extents_include = flags;
    }
    Ok(())
}
