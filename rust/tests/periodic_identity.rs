//! Periodic bond images and the symmetry site map (spec § 5.1, § 5.7).
//!
//! A bond of a periodic structure may join atom i to an image of atom j in a
//! neighbouring cell (`pairs[].image`), and a `structure.symmetry` payload may
//! say which site, operation and lattice translation generated each atom
//! (`sites`, `atom_sites`). These tests write archives with the Rust writer,
//! check that the valid ones validate and that each rule the validator
//! enforces is reported for an archive that breaks it. The cases mirror
//! `python/tests/test_periodic_identity.py`.

use qvf::{validate_bytes, ValidateOptions, Writer};
use serde_json::{json, Map, Value};

// --- fixtures -----------------------------------------------------------------------

const A_GRAPHENE: f64 = 2.46;
const A_CUBIC: f64 = 3.35;

fn graphene_lattice() -> [[f64; 3]; 3] {
    [[A_GRAPHENE, 0.0, 0.0], [-A_GRAPHENE / 2.0, A_GRAPHENE * 3f64.sqrt() / 2.0, 0.0], [0.0, 0.0, 15.0]]
}

fn cubic_lattice() -> [[f64; 3]; 3] {
    [[A_CUBIC, 0.0, 0.0], [0.0, A_CUBIC, 0.0], [0.0, 0.0, A_CUBIC]]
}

/// Each carbon bonds to the other through three images (1.42 Å each).
fn graphene_bonds() -> Vec<Value> {
    vec![
        json!({"i": 0, "j": 1, "image": [0, 0, 0]}),
        json!({"i": 0, "j": 1, "image": [-1, 0, 0]}),
        json!({"i": 0, "j": 1, "image": [0, 1, 0]}),
    ]
}

/// One atom bonded to its own six neighbours: three bonds, each written once.
fn cubic_bonds() -> Vec<Value> {
    vec![
        json!({"i": 0, "j": 0, "image": [1, 0, 0]}),
        json!({"i": 0, "j": 0, "image": [0, 1, 0]}),
        json!({"i": 0, "j": 0, "image": [0, 0, 1]}),
    ]
}

fn cart(frac: [f64; 3], lattice: &[[f64; 3]; 3]) -> [f64; 3] {
    let mut x = [0.0; 3];
    for (r, xr) in x.iter_mut().enumerate() {
        *xr = (0..3).map(|c| frac[c] * lattice[c][r]).sum();
    }
    x
}

fn atom(symbol: &str, z: u32, position: [f64; 3]) -> Value {
    json!({"symbol": symbol, "atomic_number": z, "position": position})
}

fn graphene_atoms() -> Vec<Value> {
    let l = graphene_lattice();
    vec![atom("C", 6, cart([1.0 / 3.0, 2.0 / 3.0, 0.5], &l)), atom("C", 6, cart([2.0 / 3.0, 1.0 / 3.0, 0.5], &l))]
}

fn polonium() -> Vec<Value> {
    vec![atom("Po", 84, [0.0, 0.0, 0.0])]
}

fn structure_payload(atoms: &[Value], pbc: Option<[bool; 3]>, lattice: Option<[[f64; 3]; 3]>) -> Value {
    let mut p = json!({"atoms": atoms});
    if let Some(pbc) = pbc {
        p["pbc"] = json!(pbc);
    }
    if let Some(l) = lattice {
        p["lattice_vectors"] = json!(l);
    }
    p
}

#[derive(Clone, Copy, PartialEq)]
enum Layout {
    /// The bonds are the structure section's own `bonds` member.
    Inline,
    /// A standalone `bonds` section with `structure_ref`.
    Standalone,
    /// A `bond_orders` section (the archive's only structure).
    BondOrders,
}

fn writer() -> Writer<Vec<u8>> {
    Writer::new(Vec::new(), qvf::Source::new("qvf-rust-test", "0.1.0")).unwrap()
}

fn bonds_archive(
    atoms: &[Value],
    pairs: &[Value],
    pbc: Option<[bool; 3]>,
    lattice: Option<[[f64; 3]; 3]>,
    layout: Layout,
) -> Vec<u8> {
    let mut w = writer();
    w.add_section("structure", "structure", Map::new()).unwrap();
    w.add_json("structure", "structure", &structure_payload(atoms, pbc, lattice)).unwrap();
    match layout {
        Layout::Inline => {
            w.add_json("structure", "bonds", &json!({"pairs": pairs})).unwrap();
        }
        Layout::Standalone => {
            let mut fields = Map::new();
            fields.insert("structure_ref".into(), json!("structure"));
            w.add_section("bonds", "bonds", fields).unwrap();
            w.add_json("bonds", "bonds", &json!({"pairs": pairs})).unwrap();
        }
        Layout::BondOrders => {
            let pairs: Vec<Value> = pairs
                .iter()
                .map(|p| {
                    let mut p = p.clone();
                    p["order"] = json!(1.0);
                    p
                })
                .collect();
            w.add_section("bond_orders", "bond_orders", Map::new()).unwrap();
            w.add_json("bond_orders", "bond_orders", &json!({"method": "mayer", "pairs": pairs})).unwrap();
        }
    }
    w.finish().unwrap()
}

fn graphene(pairs: &[Value], layout: Layout) -> Vec<u8> {
    bonds_archive(&graphene_atoms(), pairs, Some([true, true, false]), Some(graphene_lattice()), layout)
}

/// `(code, location)` of every error.
fn errors(bytes: Vec<u8>) -> Vec<(String, String)> {
    validate_bytes(bytes, &ValidateOptions::default())
        .findings
        .iter()
        .filter(|f| f.code.as_str().starts_with('E'))
        .map(|f| (f.code.as_str().to_string(), f.location.clone()))
        .collect()
}

fn codes(bytes: Vec<u8>) -> Vec<String> {
    let mut c: Vec<String> = errors(bytes).into_iter().map(|e| e.0).collect();
    c.sort();
    c.dedup();
    c
}

fn assert_clean(bytes: Vec<u8>) {
    let report = validate_bytes(bytes, &ValidateOptions::default());
    let all: Vec<String> = report.findings.iter().map(|f| format!("{} {} {}", f.code, f.location, f.message)).collect();
    assert!(report.ok && all.is_empty(), "{all:#?}");
}

fn malformed(at: &str) -> Vec<(String, String)> {
    vec![("E_MEMBER_MALFORMED".to_string(), at.to_string())]
}

// --- periodic bond images: archives that validate -------------------------------------

#[test]
fn graphene_cell_bonds_through_three_images() {
    assert_clean(graphene(&graphene_bonds(), Layout::Inline));
    assert_clean(graphene(&graphene_bonds(), Layout::Standalone));
}

#[test]
fn cubic_cell_bonds_to_its_own_images() {
    assert_clean(bonds_archive(&polonium(), &cubic_bonds(), Some([true; 3]), Some(cubic_lattice()), Layout::Inline));
}

#[test]
fn molecule_bonds_without_images_keep_their_meaning() {
    let water = [atom("O", 8, [0.0, 0.0, 0.117]), atom("H", 1, [0.757, 0.0, -0.470]), atom("H", 1, [-0.757, 0.0, -0.470])];
    // A repeated pair without images stays valid, as it was before images.
    let pairs = [json!({"i": 0, "j": 1}), json!({"i": 0, "j": 2}), json!({"i": 1, "j": 0})];
    assert_clean(bonds_archive(&water, &pairs, None, None, Layout::Inline));
}

#[test]
fn bond_orders_carry_images() {
    assert_clean(graphene(&graphene_bonds(), Layout::BondOrders));
}

#[test]
fn image_along_a_leading_aperiodic_axis_is_zero() {
    // The periodic axes need not lead (spec § 5.1).
    let lattice = [[20.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]];
    let atoms = [atom("H", 1, [0.0, 0.0, 0.0])];
    let pairs = [json!({"i": 0, "j": 0, "image": [0, 0, 1]})];
    assert_clean(bonds_archive(&atoms, &pairs, Some([false, true, true]), Some(lattice), Layout::Inline));
}

// --- periodic bond images: archives that must fail -----------------------------------

#[test]
fn self_bond_with_a_zero_image_is_refused() {
    let pairs = [json!({"i": 1, "j": 1, "image": [0, 0, 0]})];
    assert_eq!(errors(graphene(&pairs, Layout::Inline)), malformed("/sections/0/members/bonds"));
}

#[test]
fn image_that_is_not_three_integers_is_the_schema_finding() {
    for image in [json!([0.5, 0, 0]), json!([1, 0])] {
        let pairs = [json!({"i": 0, "j": 1, "image": image})];
        assert_eq!(codes(graphene(&pairs, Layout::Inline)), ["E_PAYLOAD_SCHEMA"]);
    }
}

#[test]
fn image_on_a_molecule_is_refused() {
    let atoms = [atom("H", 1, [0.0, 0.0, 0.0]), atom("H", 1, [0.0, 0.0, 0.74])];
    let pairs = [json!({"i": 0, "j": 1, "image": [0, 0, 0]})];
    assert_eq!(errors(bonds_archive(&atoms, &pairs, None, None, Layout::Standalone)), malformed("/sections/1/members/bonds"));
}

#[test]
fn image_along_an_aperiodic_axis_is_refused() {
    let pairs = [json!({"i": 0, "j": 1, "image": [0, 0, 1]})];
    assert_eq!(errors(graphene(&pairs, Layout::Inline)), malformed("/sections/0/members/bonds"));
}

#[test]
fn one_bond_written_twice_is_refused() {
    for second in [
        json!({"i": 0, "j": 1, "image": [-1, 0, 0]}), // the same record twice
        json!({"i": 1, "j": 0, "image": [1, 0, 0]}),  // its reversal
    ] {
        let pairs = [graphene_bonds()[1].clone(), second];
        assert_eq!(errors(graphene(&pairs, Layout::Inline)), malformed("/sections/0/members/bonds"));
    }
}

#[test]
fn self_bond_and_its_reversal_are_one_bond() {
    let pairs = [json!({"i": 0, "j": 0, "image": [1, 0, 0]}), json!({"i": 0, "j": 0, "image": [-1, 0, 0]})];
    let bytes = bonds_archive(&polonium(), &pairs, Some([true; 3]), Some(cubic_lattice()), Layout::Inline);
    assert_eq!(errors(bytes), malformed("/sections/0/members/bonds"));
}

#[test]
fn bond_orders_follow_the_same_rules() {
    let mut reversed = graphene_bonds()[0].clone();
    reversed["i"] = json!(1);
    reversed["j"] = json!(0);
    let pairs = [graphene_bonds()[0].clone(), reversed];
    assert_eq!(errors(graphene(&pairs, Layout::BondOrders)), malformed("/sections/1/members/bond_orders"));
}

#[test]
fn bonds_of_an_unresolved_structure_are_not_checked_for_images() {
    // With two structures and no structure_ref, a bond_orders section belongs
    // to neither, so its images are not checked.
    let mut w = writer();
    for id in ["a", "b"] {
        w.add_section(id, "structure", Map::new()).unwrap();
        w.add_json(id, "structure", &structure_payload(&graphene_atoms(), None, None)).unwrap();
    }
    w.add_section("bond_orders", "bond_orders", Map::new()).unwrap();
    let pairs = json!([{"i": 0, "j": 1, "image": [1, 0, 0], "order": 1.0}]);
    w.add_json("bond_orders", "bond_orders", &json!({"method": "mayer", "pairs": pairs})).unwrap();
    assert_clean(w.finish().unwrap());
}

// --- the site map -------------------------------------------------------------------

type Op = ([[f64; 3]; 3], [f64; 3]);

/// An operation such as `x+1/2,-y,z` as (R, t) on fractional coordinates.
fn parse_op(text: &str) -> Op {
    let mut rot = [[0.0; 3]; 3];
    let mut shift = [0.0; 3];
    for (r, part) in text.split(',').enumerate() {
        let chars: Vec<char> = part.chars().collect();
        let mut sign = 1.0;
        let mut i = 0;
        while i < chars.len() {
            let ch = chars[i];
            match ch {
                '+' => sign = 1.0,
                '-' => sign = -1.0,
                'x' | 'y' | 'z' => {
                    rot[r]["xyz".find(ch).unwrap()] = sign;
                    sign = 1.0;
                }
                _ if ch.is_ascii_digit() => {
                    let mut j = i;
                    while j < chars.len() && (chars[j].is_ascii_digit() || chars[j] == '/') {
                        j += 1;
                    }
                    let frac: String = chars[i..j].iter().collect();
                    let (num, den) = frac.split_once('/').unwrap();
                    shift[r] += sign * num.parse::<f64>().unwrap() / den.parse::<f64>().unwrap();
                    sign = 1.0;
                    i = j - 1;
                }
                _ => {}
            }
            i += 1;
        }
    }
    (rot, shift)
}

/// Expand `sites` (label, symbol, Z, fractional position) with the operations
/// into the reference cell, removing copies that coincide, and record each
/// atom's site, operation and translation.
fn expand(
    op_texts: &[&str],
    sites: &[(&str, &str, u32, [f64; 3])],
    lattice: &[[f64; 3]; 3],
) -> (Vec<Value>, Vec<Value>, Vec<Value>) {
    let ops: Vec<Op> = op_texts.iter().map(|t| parse_op(t)).collect();
    let (mut atoms, mut atom_sites) = (Vec::new(), Vec::new());
    for (s, (_, symbol, z, pos)) in sites.iter().enumerate() {
        let mut seen: Vec<[f64; 3]> = Vec::new();
        for (o, (rot, shift)) in ops.iter().enumerate() {
            let raw: Vec<f64> = (0..3).map(|r| (0..3).map(|c| rot[r][c] * pos[c]).sum::<f64>() + shift[r]).collect();
            let n: Vec<i64> = raw.iter().map(|x| -((x + 1e-9).floor() as i64)).collect();
            let frac = [raw[0] + n[0] as f64, raw[1] + n[1] as f64, raw[2] + n[2] as f64];
            if seen.iter().any(|f| (0..3).map(|r| (frac[r] - f[r]).abs()).fold(0.0, f64::max) < 1e-6) {
                continue;
            }
            seen.push(frac);
            atoms.push(atom(symbol, *z, cart(frac, lattice)));
            atom_sites.push(json!({"site": s, "operation": o, "translation": n}));
        }
    }
    let operations = ops.iter().map(|(r, t)| json!({"rotation": r, "translation": t})).collect();
    (atoms, operations, atom_sites)
}

/// P 42/m n m (No. 136), the sixteen operations of the standard setting.
const P42MNM: [&str; 16] = [
    "x,y,z",
    "-x,-y,z",
    "-y+1/2,x+1/2,z+1/2",
    "y+1/2,-x+1/2,z+1/2",
    "-x+1/2,y+1/2,-z+1/2",
    "x+1/2,-y+1/2,-z+1/2",
    "y,x,-z",
    "-y,-x,-z",
    "-x,-y,-z",
    "x,y,-z",
    "y+1/2,-x+1/2,-z+1/2",
    "-y+1/2,x+1/2,-z+1/2",
    "x+1/2,-y+1/2,z+1/2",
    "-x+1/2,y+1/2,z+1/2",
    "-y,-x,z",
    "y,x,z",
];
const RUTILE_LATTICE: [[f64; 3]; 3] = [[4.594, 0.0, 0.0], [0.0, 4.594, 0.0], [0.0, 0.0, 2.959]];
/// Ti on 2a, a special position; O on 4f, deposited outside the reference
/// cell so that its images carry nonzero translations.
const RUTILE_SITES: [(&str, &str, u32, [f64; 3]); 2] =
    [("Ti1", "Ti", 22, [0.0, 0.0, 0.0]), ("O1", "O", 8, [1.3049, 0.3049, 0.0])];

/// P 1 21/n 1, a non-default setting of P 21/c (No. 14), in a skew cell.
const P21N: [&str; 4] = ["x,y,z", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x+1/2,-y+1/2,z+1/2"];
const P21N_SITES: [(&str, &str, u32, [f64; 3]); 1] = [("C1", "C", 6, [0.11, 0.22, 0.33])];

fn p21n_lattice() -> [[f64; 3]; 3] {
    let beta = 97f64.to_radians();
    [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [7.0 * beta.cos(), 0.0, 7.0 * beta.sin()]]
}

fn write_symmetry(structure: &Value, data: &Value) -> Vec<u8> {
    let mut w = writer();
    w.add_section("structure", "structure", Map::new()).unwrap();
    w.add_json("structure", "structure", structure).unwrap();
    w.add_section("symmetry", "structure.symmetry", Map::new()).unwrap();
    w.add_json("symmetry", "data", data).unwrap();
    w.finish().unwrap()
}

fn symmetry_archive(
    op_texts: &[&str],
    sites: &[(&str, &str, u32, [f64; 3])],
    lattice: [[f64; 3]; 3],
    pbc: [bool; 3],
    edit: impl FnOnce(&mut Value),
    atoms_edit: impl FnOnce(&mut Vec<Value>),
) -> Vec<u8> {
    let (mut atoms, operations, atom_sites) = expand(op_texts, sites, &lattice);
    let site_list: Vec<Value> = sites.iter().map(|(label, _, _, pos)| json!({"label": label, "position": pos})).collect();
    let mut data = json!({"source": "deposited", "operations": operations, "sites": site_list, "atom_sites": atom_sites});
    edit(&mut data);
    atoms_edit(&mut atoms);
    write_symmetry(&structure_payload(&atoms, Some(pbc), Some(lattice)), &data)
}

fn rutile(edit: impl FnOnce(&mut Value)) -> Vec<u8> {
    symmetry_archive(&P42MNM, &RUTILE_SITES, RUTILE_LATTICE, [true; 3], edit, |_| {})
}

// --- the site map: archives that validate ---------------------------------------------

#[test]
fn rutile_site_map_on_a_special_position() {
    let (atoms, _, atom_sites) = expand(&P42MNM, &RUTILE_SITES, &RUTILE_LATTICE);
    assert_eq!(atoms.len(), 6); // Ti on 2a, O on 4f
    let sites: Vec<i64> = atom_sites.iter().map(|e| e["site"].as_i64().unwrap()).collect();
    assert_eq!(sites, [0, 0, 1, 1, 1, 1]);
    assert!(atom_sites.iter().any(|e| e["translation"] != json!([0, 0, 0])));
    assert_clean(rutile(|_| {}));
}

#[test]
fn non_default_setting_in_a_skew_cell() {
    assert_clean(symmetry_archive(&P21N, &P21N_SITES, p21n_lattice(), [true; 3], |_| {}, |_| {}));
}

#[test]
fn molecule_site_map_in_cartesian_coordinates() {
    let mut data = json!({
        "point_group": "C2v", "source": "detected", "symprec": 1e-3,
        "operations": [{"rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},
                       {"rotation": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]}],
        "sites": [{"label": "O1", "position": [0.0, 0.0, 0.117]},
                  {"label": "H1", "position": [0.757, 0.0, -0.470]}],
        "atom_sites": [{"site": 0, "operation": 0}, {"site": 1, "operation": 0},
                       {"site": 1, "operation": 1, "translation": [0, 0, 0]}]
    });
    let water = structure_payload(
        &[atom("O", 8, [0.0, 0.0, 0.117]), atom("H", 1, [0.757, 0.0, -0.470]), atom("H", 1, [-0.757, 0.0, -0.470])],
        None,
        None,
    );
    assert_clean(write_symmetry(&water, &data));
    // A molecule has no periodic axis, so any lattice translation is refused.
    data["atom_sites"][2]["translation"] = json!([1, 0, 0]);
    assert_eq!(errors(write_symmetry(&water, &data)), malformed("/sections/1/members/data"));
}

#[test]
fn symmetry_without_a_site_map_is_unchanged() {
    assert_clean(rutile(|d| {
        let o = d.as_object_mut().unwrap();
        o.remove("atom_sites");
        o.remove("sites");
    }));
}

#[test]
fn site_tolerance_widens_the_position_check() {
    let shift = |atoms: &mut Vec<Value>| {
        let x = atoms[3]["position"][0].as_f64().unwrap();
        atoms[3]["position"][0] = json!(x + 0.05);
    };
    let tight = symmetry_archive(&P42MNM, &RUTILE_SITES, RUTILE_LATTICE, [true; 3], |_| {}, shift);
    assert_eq!(errors(tight), malformed("/sections/1/members/data"));
    let wide = symmetry_archive(&P42MNM, &RUTILE_SITES, RUTILE_LATTICE, [true; 3], |d| d["site_tolerance"] = json!(0.1), shift);
    assert_clean(wide);
}

// --- the site map: archives that must fail -------------------------------------------

#[test]
fn site_map_rules() {
    type Edit = Box<dyn FnOnce(&mut Value)>;
    let entry = |key: &'static str, value: Value| -> Edit { Box::new(move |d: &mut Value| d["atom_sites"][2][key] = value) };
    let remove = |key: &'static str| -> Edit {
        Box::new(move |d: &mut Value| {
            d.as_object_mut().unwrap().remove(key);
        })
    };
    let cases: Vec<(&str, Edit, &str)> = vec![
        ("operation out of range", entry("operation", json!(16)), "E_MEMBER_MALFORMED"),
        ("site out of range", entry("site", json!(2)), "E_MEMBER_MALFORMED"),
        ("a wrong position", entry("operation", json!(3)), "E_MEMBER_MALFORMED"),
        ("fractional translation", entry("translation", json!([0, 0, 0.5])), "E_PAYLOAD_SCHEMA"),
        ("translation moves the atom", entry("translation", json!([1, 0, 0])), "E_MEMBER_MALFORMED"),
        (
            "one entry short",
            Box::new(|d: &mut Value| {
                d["atom_sites"].as_array_mut().unwrap().pop();
            }),
            "E_MEMBER_MALFORMED",
        ),
        ("no sites", remove("sites"), "E_MEMBER_MALFORMED"),
        ("no operations", remove("operations"), "E_MEMBER_MALFORMED"),
        ("unknown source", Box::new(|d: &mut Value| d["source"] = json!("guessed")), "E_PAYLOAD_SCHEMA"),
        (
            "site without a position",
            Box::new(|d: &mut Value| {
                d["sites"][0].as_object_mut().unwrap().remove("position");
            }),
            "E_PAYLOAD_SCHEMA",
        ),
    ];
    for (what, edit, code) in cases {
        assert_eq!(codes(rutile(edit)), [code], "{what}");
    }
}

#[test]
fn translation_along_an_aperiodic_axis_is_refused() {
    let bytes = symmetry_archive(
        &P42MNM,
        &RUTILE_SITES,
        RUTILE_LATTICE,
        [true, true, false],
        |d| d["atom_sites"][0]["translation"] = json!([0, 0, 1]),
        |atoms| {
            let z = atoms[0]["position"][2].as_f64().unwrap();
            atoms[0]["position"][2] = json!(z + RUTILE_LATTICE[2][2]);
        },
    );
    assert_eq!(errors(bytes), malformed("/sections/1/members/data"));
}

#[test]
fn an_edited_structure_disagrees_with_its_site_map() {
    // A writer that moves an atom must rewrite or drop the map (spec § 5.7).
    let bytes = symmetry_archive(
        &P42MNM,
        &RUTILE_SITES,
        RUTILE_LATTICE,
        [true; 3],
        |_| {},
        |atoms| {
            atoms[0]["position"] = json!([1.0, 1.0, 1.0]);
        },
    );
    assert_eq!(errors(bytes), malformed("/sections/1/members/data"));
}

#[test]
fn checks_run_at_the_references_level_only() {
    // As in the Python reference, which runs them in its references stage.
    let pairs = [json!({"i": 1, "j": 1, "image": [0, 0, 0]})];
    let at = |level| {
        let options = ValidateOptions { level, ..ValidateOptions::default() };
        validate_bytes(graphene(&pairs, Layout::Inline), &options).has(qvf::Code::E_MEMBER_MALFORMED)
    };
    assert!(!at(qvf::Level::Payloads));
    assert!(!at(qvf::Level::Semantics));
    assert!(at(qvf::Level::References));
    assert!(at(qvf::Level::All));
    let bytes = rutile(|d| d["atom_sites"][2]["operation"] = json!(16));
    let options = ValidateOptions { level: qvf::Level::Semantics, ..ValidateOptions::default() };
    assert!(!validate_bytes(bytes.clone(), &options).has(qvf::Code::E_MEMBER_MALFORMED));
    assert!(validate_bytes(bytes, &ValidateOptions::default()).has(qvf::Code::E_MEMBER_MALFORMED));
}
