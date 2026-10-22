//! End-to-end check of the InChI -> structure restoration
//! (`AllInchiToStructure` and everything below it): a standard InChI is
//! read, every component is restored and the InChI of the restored
//! structure must have the same formula and charge.

use rdoxide::external::inchi::ichi_io::{
    INCHI_IOS_TYPE_STRING, InchiIostream, IosFile, inchi_ios_init,
};
use rdoxide::external::inchi::ichiparm::{InputParms, read_command_line_parms};
use rdoxide::external::inchi::ichiread::*;
use rdoxide::external::inchi::ichirvr2::StrFromInChI;
use rdoxide::external::inchi::ichirvr7::{
    add_proton_and_iso_h_balance_to_mob_h_struct, all_inchi_to_structure,
    compare_all_disconnected_orig_inchi_to_rev_inchi, compare_all_orig_inchi_to_rev_inchi,
};
use rdoxide::external::inchi::ichirvrs::set_up_srm;
use rdoxide::external::inchi::runichi3::StructData;

fn parse(s: &str) -> Option<InpInChI> {
    let mut ios = InchiIostream::default();
    inchi_ios_init(&mut ios, INCHI_IOS_TYPE_STRING, IosFile::None);
    let bytes = s.as_bytes().to_vec();
    ios.s.n_used_length = bytes.len() as i32;
    ios.s.p_str = Some(bytes);
    ios.s.n_ptr = 0;
    let mut line = SegmLine::default();
    let mut pstr = None;
    let (mut state, mut nerr) = (-77, -9);
    let mut d = InpInChI::default();
    let (mut stdf, mut hso, mut bits) = (-5, -5, 0u8);
    let r = inchi_line2data(
        &mut ios, &mut line, &mut pstr, &mut state, &mut nerr, &mut d, true, 0, &mut stdf,
        &mut hso, &mut bits,
    );
    ((r == RI_ERR_EOL || r == RI_ERR_EOF) && nerr == 0).then_some(d)
}

fn sig(c: &rdoxide::external::inchi::ichi::InChI) -> String {
    let na = c.n_number_of_atoms as usize;
    let st = |x: &Option<rdoxide::external::inchi::ichi::InChIStereo>| {
        x.as_ref()
            .filter(|x| x.n_number_of_stereo_centers + x.n_number_of_stereo_bonds != 0)
            .map(|x| {
                let (nc, nb) = (
                    x.n_number_of_stereo_centers as usize,
                    x.n_number_of_stereo_bonds as usize,
                );
                format!(
                    "{:?}{:?}{:?}{:?}{:?}",
                    &x.n_number[..nc],
                    &x.t_parity[..nc],
                    &x.n_bond_atom1[..nb],
                    &x.n_bond_atom2[..nb],
                    &x.b_parity[..nb]
                )
            })
    };
    format!(
        "{:?} q{} {:?} {:?} {:?} {:?} {:?} {:?}",
        c.sz_hill_formula,
        c.n_total_charge,
        &c.n_atom[..na],
        &c.n_conn_table[..c.len_conn_table as usize],
        &c.n_num_h[..na],
        // a lone 0 (no t-groups) is how a recalculated structure with added
        // protons stores an empty t-group list
        c.n_tautomer[..c.len_tautomer as usize]
            .iter()
            .copied()
            .filter(|_| c.len_tautomer > 1)
            .collect::<Vec<_>>(),
        st(&c.stereo),
        st(&c.stereo_isotopic)
    )
}

fn default_parms() -> InputParms {
    let mut ip = InputParms::default();
    let mut log = InchiIostream::default();
    inchi_ios_init(&mut log, INCHI_IOS_TYPE_STRING, IosFile::None);
    let r = read_command_line_parms(&["inchi"], &mut ip, true, &mut log);
    assert_eq!(r, 0);
    ip
}

#[test]
fn restore_simple_inchi() {
    let cases = [
        "InChI=1S/CH4/h1H4",
        "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
        "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)",
        "InChI=1S/H2O/h1H2",
        "InChI=1S/C2H4/c1-2/h1-2H2",
        "InChI=1S/CH5N/c1-2/h2H2,1H3",
        "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1D",
        "InChI=1S/C4H8O2/c1-3-4(5)6-2/h3H2,1-2H3",
        "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
        "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2",
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)",
        "InChI=1S/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)",
        "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)",
        "InChI=1S/ClH/h1H",
        "InChI=1S/Na.ClH/h;1H/q+1;/p-1",
        "InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3",
        "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+",
        "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-",
        "InChI=1S/C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3/t4-/m1/s1",
        "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3",
        "InChI=1S/C6H5NO2/c8-7(9)6-4-2-1-3-5-6/h1-5H",
        "InChI=1S/C6H7N/c7-6-4-2-1-3-5-6/h1-5H,7H2",
        "InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
        "InChI=1S/H3N/h1H3/p+1",
        "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1",
        "InChI=1S/CH3O/c1-2/h1H3/q-1",
        "InChI=1S/O4S/c1-5(2,3)4/q-2",
        "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
        "InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)",
    ];
    let ip = default_parms();
    let srm = set_up_srm();
    let mut n_restored = 0;
    let mut bad: Vec<String> = Vec::new();
    for s in cases {
        let Some(mut d) = parse(s) else {
            bad.push(format!("UNREADABLE {s}"));
            continue;
        };
        let mut ps: [[Vec<StrFromInChI>; 2]; 2] = Default::default();
        let ne = all_inchi_to_structure(&ip, 1, &srm, false, &mut ps, &mut d, 1, None);
        assert!(ne >= 0, "{s}: fatal {ne}");
        // no Fixed-H layer: the proton balance is added to the Mobile-H structures
        let rb = add_proton_and_iso_h_balance_to_mob_h_struct(
            &ip,
            &mut StructData::default(),
            1,
            false,
            &mut ps,
            &d,
        );
        assert!(rb >= 0, "{s}: balance {rb}");
        let rc = compare_all_orig_inchi_to_rev_inchi(&ps, &mut d, false);
        assert!(rc >= 0, "{s}: compare {rc}");
        if d.compare_inchi_flags[0] != [0, 0] {
            bad.push(format!("{s}: compare flags {:?}", d.compare_inchi_flags[0]));
        }
        // no reconnected layer: nothing to compare
        assert_eq!(
            compare_all_disconnected_orig_inchi_to_rev_inchi(&ps, &mut d, false),
            0,
            "{s}"
        );
        assert_eq!(d.compare_inchi_flags[1], [0, 0], "{s}");
        let comp = &ps[0][1][0];
        eprintln!(
            "{s}: errors {ne}, n_error {}, atoms {}, at2 {}, rev {}",
            comp.n_error,
            comp.num_atoms,
            comp.at2.len(),
            comp.rev_inchi_component_exists(0, 1, 0)
        );
        if ne == 0 && comp.rev_inchi_component_exists(0, 1, 0) {
            let rev = comp.rev_inchi.p_inchi[0][0][1].as_ref().unwrap();
            let orig = &d.p_inp_inchi[0][1][0];
            if sig(rev) != sig(orig) {
                bad.push(format!("{s}\n   got  {}\n   want {}", sig(rev), sig(orig)));
                continue;
            }
            n_restored += 1;
        }
    }
    eprintln!("restored: {n_restored}/{}", cases.len());
    for b in &bad {
        eprintln!("MISMATCH {b}");
    }
    assert!(bad.is_empty(), "restored InChI differs");
    assert_eq!(n_restored, cases.len());
}

#[test]
fn compare_reports_differences() {
    use rdoxide::external::inchi::ichirvr7::diff;
    let ip = default_parms();
    let mut d = parse("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").unwrap();
    let mut ps: [[Vec<StrFromInChI>; 2]; 2] = Default::default();
    assert_eq!(
        all_inchi_to_structure(&ip, 1, &set_up_srm(), false, &mut ps, &mut d, 1, None),
        0
    );
    // the original InChI claims one more H on the first atom
    d.p_inp_inchi[0][1][0].n_num_h[0] += 1;
    assert!(compare_all_orig_inchi_to_rev_inchi(&ps, &mut d, false) >= 0);
    assert_ne!(d.compare_inchi_flags[0][1] & diff::POSITION_H, 0);
}
