//! Differential tests of the Molfile reader (`mol_fmt1.c`, `mol_fmt2.c`)
//! against fixtures generated from the upstream C code
//! (`tools/inchi_molfmt1_reference.py`).

use rdoxide::external::inchi::ichi_io::{INCHI_IOS_TYPE_STRING, InchiIostream};
use rdoxide::external::inchi::molfmt::*;

fn unhex(h: &str) -> Vec<u8> {
    if h == "-" {
        return Vec::new();
    }
    (0..h.len() / 2)
        .map(|i| u8::from_str_radix(&h[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn hex(s: &[u8]) -> String {
    if s.is_empty() {
        "-".into()
    } else {
        s.iter().map(|c| format!("{c:02x}")).collect()
    }
}

fn pd(d: f64) -> String {
    if d.is_nan() {
        "nan".into()
    } else {
        format!("{:016x}", d.to_bits())
    }
}

fn stream(text: &[u8]) -> InchiIostream {
    let mut ios = InchiIostream {
        type_: INCHI_IOS_TYPE_STRING,
        ..Default::default()
    };
    ios.s.p_str = Some(text.to_vec());
    ios.s.n_used_length = text.len() as i32;
    ios.s.n_allocated_length = text.len() as i32 + 2;
    ios
}

fn field_case(t: char, len: i32, text: &[u8]) -> String {
    let mut pos = 0;
    let (ret, val) = match t {
        'S' => {
            let mut d = Vec::new();
            let r = molfile_read_field_str(&mut d, len, text, &mut pos);
            (r, hex(&d))
        }
        'C' | 'N' | 'L' => {
            let kind = match t {
                'C' => IntKind::Char,
                'N' => IntKind::Short,
                _ => IntKind::Long,
            };
            let (v, r) = molfile_read_field_int(kind, len, text, &mut pos);
            (r, v.to_string())
        }
        'D' => {
            let (v, r) = molfile_read_field_double(len, text, &mut pos);
            (r, pd(v))
        }
        'F' => {
            let (v, r) = molfile_read_field_float(len, text, &mut pos);
            (r, pd(f64::from(v)))
        }
        _ => (molfile_read_field_jump(len, text, &mut pos), "-".into()),
    };
    format!("{ret} {pos} {val}")
}

fn mol_case(bcoord: bool, text: &[u8]) -> Vec<String> {
    let mut out = Vec::new();
    let mut ios = stream(text);
    let mut err = String::new();
    let mut m = MolFmtData::default();
    let mut retcode = molfile_read_header_lines(&mut m.hdr, &mut ios, &mut err);
    let mut data_ended = false;
    let h = &m.hdr;
    out.push(format!(
        "H {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        hex(&h.molname),
        hex(&h.line2),
        hex(&h.user_initls),
        hex(&h.prog_name),
        h.month,
        h.day,
        h.year,
        h.hour,
        h.minute,
        hex(&h.dim_code),
        h.scaling_factor1,
        pd(h.scaling_factor2),
        pd(h.energy),
        h.internal_regno,
        hex(&h.comment),
        retcode
    ));
    out.push(format!("N {}", molfile_extract_struc_num(h)));
    'fin: {
        if retcode != 0 {
            retcode += 10;
            break 'fin;
        }
        retcode = molfile_read_counts_line(&mut m.ctab, &mut ios, &mut err);
        let c = &m.ctab;
        out.push(format!(
            "C {} {} {} {} {} {} {}",
            c.n_atoms,
            c.n_bonds,
            c.chiral_flag,
            c.n_stext_entries,
            c.n_property_lines,
            hex(&c.version_string),
            i32::from(c.v3000.is_some())
        ));
        if retcode != 0 {
            retcode += 20;
            break 'fin;
        }
        if m.ctab.v3000.is_some() {
            break 'fin;
        }
        let na = m.ctab.n_atoms.max(1) as usize;
        m.ctab.atoms = Some(vec![MolFmtAtom::default(); na]);
        if bcoord {
            m.ctab.coords = Some(vec![[0u8; 32]; na]);
        }
        retcode = molfile_read_atoms_block(&mut m.ctab, &mut ios, retcode, &mut err);
        if retcode != 0 {
            if retcode < 0 {
                retcode = -retcode;
                data_ended = true;
            }
            retcode += 30;
        }
        if retcode < 30 && !data_ended {
            m.ctab.bonds = Some(vec![MolFmtBond::default(); m.ctab.n_bonds.max(1) as usize]);
        }
        let mut prevcode = retcode;
        if !data_ended {
            retcode = molfile_read_bonds_block(&mut m.ctab, &mut ios, retcode, &mut err);
            if retcode != 0 {
                if retcode < 0 {
                    retcode = -retcode;
                    data_ended = true;
                }
                retcode = if prevcode != 0 {
                    prevcode
                } else {
                    retcode + 40
                };
            }
        }
        prevcode = retcode;
        if !data_ended {
            retcode = molfile_read_stext_block(&m.ctab, &mut ios, retcode, &mut err);
            if retcode != 0 {
                retcode = if prevcode != 0 {
                    prevcode
                } else {
                    retcode + 50
                };
            }
        }
        let c = &m.ctab;
        for i in 0..(c.n_atoms.max(0) as usize).min(300) {
            let a = &c.atoms.as_ref().unwrap()[i];
            out.push(format!(
                "A {i} {} {} {} {} {} {} {} {} {}",
                pd(a.fx),
                pd(a.fy),
                pd(a.fz),
                hex(&a.symbol),
                a.mass_difference,
                a.charge,
                a.radical,
                a.stereo_parity,
                a.valence
            ));
            if let Some(co) = &c.coords {
                let n = co[i].iter().position(|&x| x == 0).unwrap_or(32);
                out.push(format!("O {i} {}", hex(&co[i][..n])));
            }
        }
        if let Some(b) = &c.bonds {
            for i in 0..(c.n_bonds.max(0) as usize).min(300) {
                out.push(format!(
                    "B {i} {} {} {} {}",
                    b[i].atnum1, b[i].atnum2, b[i].bond_type, b[i].bond_stereo
                ));
            }
        }
        out.push(format!(
            "HN {}",
            i32::from(molfile_has_no_chem_struc(Some(&m)))
        ));
        for f in [false, true] {
            let mut e2 = 0;
            let mut e2s = String::new();
            let x = molfile_get_xyz_dim_and_norm_factors(&m, f, &mut e2, &mut e2s);
            out.push(format!(
                "X {} {} {} {} {} {} {} {} {} {}",
                x.num_dimensions,
                e2,
                pd(x.x0),
                pd(x.y0),
                pd(x.z0),
                pd(x.xmin),
                pd(x.ymin),
                pd(x.zmin),
                pd(x.scaler),
                hex(e2s.as_bytes())
            ));
        }
    }
    out.push(format!(
        "R {} {} {}",
        retcode,
        i32::from(data_ended),
        hex(err.as_bytes())
    ));
    out
}

#[test]
fn molfile_reader_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_molfmt1.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let got: Vec<String> = if t[0] == "F" {
            vec![field_case(
                t[1].chars().next().unwrap(),
                t[2].parse().unwrap(),
                &unhex(t[3]),
            )]
        } else {
            mol_case(t[1] == "1", &unhex(t[2]))
        };
        assert_eq!(got, want, "case {ci}: {inp}");
    }
}

fn dump_hdr(h: &MolFmtHeaderBlock, out: &mut Vec<String>) {
    out.push(format!(
        "H {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        hex(&h.molname),
        hex(&h.line2),
        hex(&h.user_initls),
        hex(&h.prog_name),
        h.month,
        h.day,
        h.year,
        h.hour,
        h.minute,
        hex(&h.dim_code),
        h.scaling_factor1,
        pd(h.scaling_factor2),
        pd(h.energy),
        h.internal_regno,
        hex(&h.comment)
    ));
}

fn dump_ctab(c: &MolFmtCtab, out: &mut Vec<String>) {
    out.push(format!(
        "C {} {} {} {} {} {} {}",
        c.n_atoms,
        c.n_bonds,
        c.chiral_flag,
        c.n_stext_entries,
        c.n_property_lines,
        hex(&c.version_string),
        i32::from(c.v3000.is_some())
    ));
}

fn dump_all(c: &MolFmtCtab, out: &mut Vec<String>) {
    for i in 0..(c.n_atoms.max(0) as usize).min(300) {
        let a = &c.atoms.as_ref().unwrap()[i];
        out.push(format!(
            "A {i} {} {} {} {} {} {} {} {} {} {}",
            pd(a.fx),
            pd(a.fy),
            pd(a.fz),
            hex(&a.symbol),
            a.mass_difference,
            a.charge,
            a.radical,
            a.stereo_parity,
            a.valence,
            a.atom_aliased_flag
        ));
        if let Some(co) = &c.coords {
            let n = co[i].iter().position(|&x| x == 0).unwrap_or(32);
            out.push(format!("O {i} {}", hex(&co[i][..n])));
        }
    }
    if let Some(b) = &c.bonds {
        for i in 0..(c.n_bonds.max(0) as usize).min(300) {
            out.push(format!(
                "B {i} {} {} {} {}",
                b[i].atnum1, b[i].atnum2, b[i].bond_type, b[i].bond_stereo
            ));
        }
    }
    for s in &c.sgroups {
        let mut l = format!(
            "G {} {} {} {} {} ",
            s.id, s.type_, s.subtype, s.conn, s.label
        );
        for x in s.xbr1.iter().chain(s.xbr2.iter()) {
            l += &format!("{} ", pd(*x));
        }
        l += &hex(&s.smt);
        l += &format!(" {}", s.alist.len());
        for x in &s.alist {
            l += &format!(" {x}");
        }
        l += &format!(" {}", s.blist.len());
        for x in &s.blist {
            l += &format!(" {x}");
        }
        out.push(l);
    }
}

#[test]
fn molfile_data_lines_match_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_molfmt2.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let mode: i32 = t[0].parse().unwrap();
        let (bcoord, polymers, nowarn, pseudo) =
            (t[1] == "1", t[2] == "1", t[3] == "1", t[4] == "1");
        let mut ios = stream(&unhex(t[5]));
        let mut perr = String::new();
        let mut out = Vec::new();
        let err;
        if mode == 0 {
            let (m, e) =
                molfile_read_data_lines(&mut ios, None, bcoord, polymers, &mut perr, nowarn);
            err = e;
            if let Some(mut m) = m {
                dump_hdr(&m.hdr, &mut out);
                dump_ctab(&m.ctab, &mut out);
                let mut e2 = 0;
                let mut s2 = String::new();
                let nzz = molfile_treat_pseudo_element_atoms(&mut m.ctab, pseudo, &mut e2, &mut s2);
                dump_all(&m.ctab, &mut out);
                out.push(format!("Z {nzz} {e2} {}", hex(s2.as_bytes())));
            }
        } else {
            let mut only = MolFmtData::default();
            err = molfile_read_data_lines(
                &mut ios,
                Some(&mut only),
                bcoord,
                polymers,
                &mut perr,
                nowarn,
            )
            .1;
            dump_hdr(&only.hdr, &mut out);
            dump_ctab(&only.ctab, &mut out);
        }
        out.push(format!("R {err} {}", hex(perr.as_bytes())));
        assert_eq!(out, want, "case {ci}: {inp}");
    }
}

fn dump_v3(c: &MolFmtCtab, out: &mut Vec<String>) {
    let Some(v) = c.v3000.as_ref() else {
        return;
    };
    let mut l = format!(
        "V {} {} {} {} {} {} {} {} {} {} {}",
        v.n_non_star_atoms,
        v.n_star_atoms,
        v.n_sgroups,
        v.n_3d_constraints,
        v.n_collections,
        v.n_non_haptic_bonds,
        v.n_haptic_bonds,
        v.n_steabs,
        v.n_sterel,
        v.n_sterac,
        i32::from(!v.atom_index_orig.is_empty())
    );
    if !v.atom_index_orig.is_empty() {
        for i in 0..(v.n_non_star_atoms + v.n_star_atoms) as usize {
            l += &format!(" {}:{}", v.atom_index_orig[i], v.atom_index_fin[i]);
        }
    }
    out.push(l);
    fn lists(ls: &[Vec<i32>], haptic: bool) -> String {
        let mut s = format!(" {}", ls.len());
        for x in ls {
            let n = if haptic { x[2] } else { x[1] };
            let cnt = (n + if haptic { 3 } else { 2 }).max(0) as usize;
            s += " [";
            for k in 0..cnt.min(40) {
                s += &format!("{},", x[k]);
            }
            s += "]";
        }
        s
    }
    out.push(format!("W{}", lists(&v.haptic_bonds, true)));
    out.push(format!(
        "S{}{}{}",
        lists(&v.steabs, false),
        lists(&v.sterel, false),
        lists(&v.sterac, false)
    ));
}

#[test]
fn molfile_v3000_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_molfmt3.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let (bcoord, polymers, nowarn, pseudo) =
            (t[1] == "1", t[2] == "1", t[3] == "1", t[4] == "1");
        let mut ios = stream(&unhex(t[5]));
        let mut perr = String::new();
        let mut out = Vec::new();
        let (m, err) = molfile_read_data_lines(&mut ios, None, bcoord, polymers, &mut perr, nowarn);
        if let Some(mut m) = m {
            dump_hdr(&m.hdr, &mut out);
            dump_ctab(&m.ctab, &mut out);
            let mut e2 = 0;
            let mut s2 = String::new();
            let nzz = molfile_treat_pseudo_element_atoms(&mut m.ctab, pseudo, &mut e2, &mut s2);
            dump_v3(&m.ctab, &mut out);
            dump_all(&m.ctab, &mut out);
            out.push(format!("Z {nzz} {e2} {}", hex(s2.as_bytes())));
        }
        out.push(format!("R {err} {}", hex(perr.as_bytes())));
        assert_eq!(out, want, "case {ci}: {inp}");
    }
}

#[test]
fn read_molfile_with_sdf_data_matches_c() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_molfmt4.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let bcoord = t[0] == "1";
        let polymers: i32 = t[1].parse().unwrap();
        let npzz: i32 = t[2].parse().unwrap();
        let nowarn = t[3] == "1";
        let (use_name, use_id, use_val) = (t[4] == "1", t[5] == "1", t[6] == "1");
        let labels: [Option<&[u8]>; 12] = [
            None,
            Some(b""),
            Some(b"CAS"),
            Some(b"CAS_NO"),
            Some(b"MyLabel"),
            Some(b"NAME"),
            Some(b"COMMENT"),
            Some(b"MolfileName"),
            Some(b"MolfileLine2"),
            Some(b"MolfileComment"),
            Some(b"MolfileIntRegNo"),
            Some(b"x"),
        ];
        let label = labels[t[7].parse::<usize>().unwrap()];
        let mut ios = stream(&unhex(t[8]));
        let mut perr = String::new();
        let mut name = Vec::new();
        let mut id: u64 = 777;
        let mut value = Vec::new();
        let (m, err) = read_molfile(
            &mut ios,
            None,
            bcoord,
            polymers,
            npzz,
            ReadMolfileExtras {
                name: use_name.then_some((&mut name, 64)),
                id: use_id.then_some(&mut id),
                sdf_label: label,
                sdf_value: use_val.then_some(&mut value),
            },
            &mut perr,
            nowarn,
        );
        let mut out = Vec::new();
        if let Some(m) = m {
            dump_hdr(&m.hdr, &mut out);
            dump_ctab(&m.ctab, &mut out);
            dump_v3(&m.ctab, &mut out);
            dump_all(&m.ctab, &mut out);
        }
        out.push(format!("N {} {id} {}", hex(&name), hex(&value)));
        out.push(format!("R {err} {}", hex(perr.as_bytes())));
        assert_eq!(out, want, "case {ci}: {inp}");
    }
}

fn dump_inp_atoms(
    at: &[rdoxide::external::inchi::inpdef::InpAtom],
    n: usize,
    out: &mut Vec<String>,
) {
    for (i, a) in at.iter().enumerate().take(n.min(100)) {
        let name = &a.elname[..a
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(a.elname.len())];
        let mut l = format!(
            "A {i} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
            hex(name),
            a.el_number,
            a.orig_at_number,
            a.valence,
            a.chem_bonds_valence,
            a.num_h,
            a.num_iso_h[0],
            a.num_iso_h[1],
            a.num_iso_h[2],
            a.iso_atw_diff,
            a.charge,
            a.radical,
            pd(a.x),
            pd(a.y)
        );
        l += &format!(" {}", pd(a.z));
        for j in 0..a.valence.max(0) as usize {
            l += &format!(" {}:{}:{}", a.neighbor[j], a.bond_type[j], a.bond_stereo[j]);
        }
        out.push(l);
    }
}

#[test]
fn molfile_to_inp_atoms_matches_c() {
    use rdoxide::external::inchi::mol2atom::*;
    use rdoxide::external::inchi::strutil::OrigAtomData;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mol2atom.json")).unwrap();
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let inp = c["in"].as_str().unwrap();
        let t: Vec<&str> = inp.split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let p = |k: usize| t[k].parse::<i32>().unwrap();
        let (use_od, merge, bcoord, noaddh) = (p(0) == 1, p(1) == 1, p(2) == 1, p(3) == 1);
        let (polymers, npzz, large, label_kind) = (p(4), p(5), p(6) == 1, p(7) as usize);
        let (use_val, use_id, nowarn) = (p(8) == 1, p(9) == 1, p(10) == 1);
        let labels: [Option<&[u8]>; 8] = [
            None,
            Some(b""),
            Some(b"CAS"),
            Some(b"x"),
            Some(b"MolfileName"),
            Some(b"MolfileLine2"),
            Some(b"MolfileComment"),
            Some(b"MolfileIntRegNo"),
        ];
        let mut ios = stream(&unhex(t[11]));
        let mut perr = String::new();
        let mut value = Vec::new();
        let mut id: u64 = 777;
        let mut molnum: i64 = -5;
        let mut flags: u64 = if large { FLAG_SET_INP_LARGE_MOLS } else { 0 };
        let mut od = use_od.then(OrigAtomData::default);
        let mut err = 0;
        let na = create_orig_inp_data_from_molfile(
            &mut ios,
            &mut od,
            merge,
            bcoord,
            noaddh,
            polymers,
            npzz,
            MolfileAtomsArgs {
                sdf_label: labels[label_kind],
                sdf_value: use_val.then_some(&mut value),
                id: use_id.then_some(&mut id),
                inp_atom_flags: Some(&mut flags),
            },
            &mut molnum,
            &mut err,
            &mut perr,
            nowarn,
        );
        let mut out = vec![format!(
            "R {na} {err} {molnum} {flags} {id} {} {}",
            hex(&value),
            hex(perr.as_bytes())
        )];
        if let (true, 0, Some(od)) = (use_od, err, od.as_ref()) {
            out.push(format!(
                "O {} {} {}",
                od.num_inp_atoms, od.num_inp_bonds, od.num_dimensions
            ));
            dump_inp_atoms(&od.at, od.num_inp_atoms, &mut out);
            if let Some(co) = &od.sz_coord {
                for i in 0..od.num_inp_atoms.min(100) {
                    let n = co[i].iter().position(|&x| x == 0).unwrap_or(32);
                    out.push(format!("K {i} {}", hex(&co[i][..n])));
                }
            }
            if let Some(pl) = &od.polymer {
                out.push(format!(
                    "P {} {} {} {} {} {} {}",
                    pl.n,
                    pl.n_pzz,
                    pl.really_do_frame_shift,
                    pl.frame_shift_scheme,
                    pl.treat,
                    pl.is_in_reconn,
                    pl.edit_repeats
                ));
                for u in &pl.units {
                    let mut l = format!(
                        "U {} {} {} {} {} {} {} ",
                        u.id, u.type_, u.subtype, u.conn, u.label, u.na, u.nb
                    );
                    for x in u.xbr1.iter().chain(u.xbr2.iter()) {
                        l += &format!("{} ", pd(*x));
                    }
                    l += &hex(&u.smt);
                    for x in &u.alist {
                        l += &format!(" {x}");
                    }
                    l += " |";
                    for x in &u.blist {
                        l += &format!(" {x}");
                    }
                    out.push(l);
                }
            }
            if let Some(v) = &od.v3000 {
                let mut l = format!(
                    "V {} {} {} {} {} {} {} {} {} {} {} {}",
                    v.n_non_star_atoms,
                    v.n_star_atoms,
                    v.n_sgroups,
                    v.n_3d_constraints,
                    v.n_collections,
                    v.n_non_haptic_bonds,
                    v.n_haptic_bonds,
                    v.n_steabs,
                    v.n_sterel,
                    v.n_sterac,
                    i32::from(!v.atom_index_orig.is_empty()),
                    i32::from(!v.atom_index_fin.is_empty())
                );
                if !v.atom_index_orig.is_empty() {
                    let raw: Vec<u8> = v
                        .atom_index_orig
                        .iter()
                        .flat_map(|x| x.to_le_bytes())
                        .collect();
                    for k in 0..16 {
                        l += &format!(" {}", raw.get(k).copied().unwrap_or(0));
                    }
                }
                out.push(l);
                fn lists(ls: &[Vec<i32>], n: i32, len_at: usize, add: i32) -> String {
                    let mut s = format!(" {n}");
                    for x in ls.iter().take(n.max(0) as usize) {
                        s += " [";
                        for k in 0..((x[len_at] + add).max(0) as usize).min(40) {
                            s += &format!("{},", x[k]);
                        }
                        s += "]";
                    }
                    s
                }
                out.push(format!(
                    "W{}",
                    lists(&v.lists_haptic_bonds, v.n_haptic_bonds, 2, 3)
                ));
                out.push(format!(
                    "S{}{}{}",
                    lists(&v.lists_steabs, v.n_steabs, 1, 2),
                    lists(&v.lists_sterel, v.n_sterel, 1, 2),
                    lists(&v.lists_sterac, v.n_sterac, 1, 2)
                ));
            }
        }
        assert_eq!(out, want, "case {ci}: {inp}");
    }
}

#[test]
fn molfile_text_to_inchi_matches_c() {
    use rdoxide::external::inchi::inchi_dll_b::make_inchi_from_molfile_text;
    let data: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/data/inchi_mol2inchi.json")).unwrap();
    let h = |s: &Option<String>| match s {
        None => "-".to_string(),
        Some(t) if t.is_empty() => "=".to_string(),
        Some(t) => hex(t.as_bytes()),
    };
    for (ci, c) in data.as_array().unwrap().iter().enumerate() {
        let t: Vec<&str> = c["in"].as_str().unwrap().split(' ').collect();
        let want: Vec<&str> = c["out"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let opts = (t[0] == "1").then(|| if t[1] == "=" { Vec::new() } else { unhex(t[1]) });
        let mol = unhex(t[2]);
        let (ret, out) = make_inchi_from_molfile_text(&mol, opts.as_deref());
        let got = [
            format!("R {ret}"),
            format!("I {}", h(&out.sz_inchi)),
            format!("X {}", h(&out.sz_aux_info)),
            format!("M {}", h(&out.sz_message)),
            format!("L {}", h(&out.sz_log)),
        ];
        assert_eq!(got.as_slice(), want.as_slice(), "case {ci}: {}", c["in"]);
    }
}
