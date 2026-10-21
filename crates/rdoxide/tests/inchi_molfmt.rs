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
