//! Rules of an `energy_decomposition` payload that JSON Schema cannot express (spec § 5.16).
//! Mirrors `python/qvf/_eda.py` and `js/src/eda.ts`.

use std::collections::{HashMap, HashSet};

use serde_json::Value;

use crate::error::Code;
use crate::json::{as_obj, int_of, num, num_of, show, str_of, Obj};

const SECTION: &str = "spec § 5.16";

/// Default sum-rule tolerance, hartree (§ 5.16).
pub const DEFAULT_SUM_TOLERANCE_EH: f64 = 1e-6;

/// Hartree per unit of an energy unit a decomposition may use (§ 4.3), legacy spellings included.
pub fn hartree_per_energy_unit(units: &str) -> Option<f64> {
    Some(match units {
        "hartree" | "Eh" | "Ha" | "Hartree" => 1.0,
        "ev" | "eV" => 1.0 / 27.211386245988,
        "kcal_per_mol" => 1.0 / 627.5094740631,
        "kj_per_mol" => 1.0 / 2625.4996394799,
        _ => return None,
    })
}

/// An `{value, units}` energy in hartree, or `None` when it is not one.
pub fn energy_eh(energy: Option<&Value>) -> Option<f64> {
    let e = as_obj(energy)?;
    let value = num_of(e.get("value")).filter(|v| v.is_finite())?;
    Some(value * hartree_per_energy_unit(str_of(e.get("units"))?)?)
}

fn list(v: Option<&Value>) -> &[Value] {
    v.and_then(Value::as_array).map(Vec::as_slice).unwrap_or(&[])
}

fn idx(v: Option<&Value>) -> Option<i64> {
    int_of(v).map(|x| x as i64)
}

/// `(index, spin)` of every orbital pair, or `None` without a pair table.
pub fn pair_keys(payload: Option<&Value>) -> Option<HashSet<(i64, String)>> {
    let pairs = as_obj(as_obj(payload)?.get("orbital_pairs"))?;
    let mut keys = HashSet::new();
    for p in list(pairs.get("pairs")) {
        if let Some(p) = p.as_object() {
            if let Some(i) = idx(p.get("index")) {
                keys.insert((i, str_of(p.get("spin")).unwrap_or("restricted").to_string()));
            }
        }
    }
    Some(keys)
}

fn energies(p: &Obj) -> Vec<(String, &Value)> {
    let mut found = Vec::new();
    for (k, f) in list(p.get("fragments")).iter().enumerate() {
        if let Some(e) = f.get("energy") {
            found.push((format!("fragments[{k}].energy"), e));
        }
    }
    if let Some(e) = as_obj(p.get("energies")).and_then(|e| e.get("system")) {
        found.push(("energies.system".into(), e));
    }
    for (k, t) in list(p.get("terms")).iter().enumerate() {
        if let Some(e) = t.get("energy") {
            found.push((format!("terms[{k}].energy"), e));
        }
    }
    for (k, r) in list(p.get("sum_rules")).iter().enumerate() {
        if let Some(e) = r.get("tolerance") {
            found.push((format!("sum_rules[{k}].tolerance"), e));
        }
    }
    if let Some(pairs) = as_obj(p.get("orbital_pairs")) {
        if let Some(e) = pairs.get("total") {
            found.push(("orbital_pairs.total".into(), e));
        }
        for (k, q) in list(pairs.get("pairs")).iter().enumerate() {
            if let Some(e) = q.get("energy") {
                found.push((format!("orbital_pairs.pairs[{k}].energy"), e));
            }
        }
    }
    found
}

fn join_ints(v: &[i64]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
}

/// `(code, message)` for every rule of § 5.16 an `energy_decomposition` payload breaks; `n_atoms`
/// is the atom count of the section's structure when it resolves.
pub fn payload_errors(payload: Option<&Value>, n_atoms: Option<usize>) -> Vec<(Code, String)> {
    let Some(p) = as_obj(payload) else { return vec![] };
    let mut errors: Vec<(Code, String)> = Vec::new();
    let bad =
        |msg: String, errors: &mut Vec<(Code, String)>| errors.push((Code::E_MEMBER_MALFORMED, format!("{msg} ({SECTION})")));

    // --- fragments ---------------------------------------------------------------------
    let fragments = list(p.get("fragments"));
    let n_frag = fragments.len() as i64;
    let atoms_of: Vec<Option<HashSet<i64>>> = fragments
        .iter()
        .map(|f| f.get("atoms").and_then(Value::as_array).map(|a| a.iter().filter_map(|x| idx(Some(x))).collect()))
        .collect();
    let top: Vec<usize> = fragments
        .iter()
        .enumerate()
        .filter(|(_, f)| f.as_object().is_some_and(|o| !o.contains_key("parent")))
        .map(|(k, _)| k)
        .collect();
    if p.get("fragments").is_some_and(Value::is_array) && top.len() < 2 {
        bad(
            format!("fragments has {} interacting fragment(s) (no parent); a decomposition needs at least two", top.len()),
            &mut errors,
        );
    }
    let mut owner: HashMap<i64, usize> = HashMap::new();
    for &k in &top {
        let mut atoms: Vec<i64> = atoms_of[k].iter().flatten().copied().collect();
        atoms.sort();
        for a in atoms {
            if let Some(&prev) = owner.get(&a) {
                bad(
                    format!("fragments[{k}] and fragments[{prev}] both contain atom {a}; interacting fragments are disjoint"),
                    &mut errors,
                );
            } else {
                owner.insert(a, k);
            }
        }
    }
    for (k, f) in fragments.iter().enumerate() {
        let Some(f) = f.as_object() else { continue };
        if let Some(n) = n_atoms {
            let mut atoms: Vec<i64> = atoms_of[k].iter().flatten().copied().collect();
            atoms.sort();
            for a in atoms {
                if a < 0 || a >= n as i64 {
                    bad(format!("fragments[{k}].atoms names atom {a}, but the structure has {n} atoms"), &mut errors);
                }
            }
        }
        if f.contains_key("parent") {
            match idx(f.get("parent")) {
                Some(pi) if pi >= 0 && pi < n_frag && pi as usize != k => {
                    let pu = pi as usize;
                    if !top.contains(&pu) {
                        bad(format!("fragments[{k}].parent {pi} is itself a sub-fragment; a sub-fragment subdivides an interacting fragment"), &mut errors);
                    } else if let (Some(mine), Some(theirs)) = (&atoms_of[k], &atoms_of[pu]) {
                        let mut extra: Vec<i64> = mine.difference(theirs).copied().collect();
                        extra.sort();
                        if !extra.is_empty() {
                            bad(
                                format!(
                                    "fragments[{k}] has atoms [{}] that its parent fragments[{pi}] does not contain",
                                    join_ints(&extra)
                                ),
                                &mut errors,
                            );
                        }
                    }
                }
                _ => bad(
                    format!("fragments[{k}].parent {} is not the index of another fragment", show(f.get("parent"))),
                    &mut errors,
                ),
            }
        }
    }

    // --- terms -------------------------------------------------------------------------
    let terms = list(p.get("terms"));
    let n_terms = terms.len() as i64;
    let mut keys: HashMap<String, usize> = HashMap::new();
    let mut interactions = 0usize;
    for (k, t) in terms.iter().enumerate() {
        let Some(t) = t.as_object() else { continue };
        let directed = t.get("directed").and_then(Value::as_bool) == Some(true);
        let mut scope = String::from("null");
        if let Some(frs) = t.get("fragments").and_then(Value::as_array) {
            let mut ints: Vec<i64> = frs.iter().filter_map(|x| idx(Some(x))).collect();
            for &f in &ints {
                if f < 0 || f >= n_frag {
                    bad(format!("terms[{k}].fragments names fragment {f}, but there are {n_frag} fragments"), &mut errors);
                }
            }
            if directed && frs.len() != 2 {
                bad(
                    format!("terms[{k}] is directed but names {} fragments; a directed term names [donor, acceptor]", frs.len()),
                    &mut errors,
                );
            }
            if !directed {
                ints.sort();
            }
            scope = format!("[{}]", join_ints(&ints));
        }
        if t.contains_key("parent") {
            match idx(t.get("parent")) {
                Some(pi) if pi >= 0 && pi < n_terms && pi as usize != k => {}
                _ => bad(format!("terms[{k}].parent {} is not the index of another term", show(t.get("parent"))), &mut errors),
            }
        }
        let key = format!(
            "{}|{scope}|{}|{}",
            t.get("name").map(Value::to_string).unwrap_or_default(),
            t.get("spin").map(Value::to_string).unwrap_or_default(),
            t.get("parent").map(Value::to_string).unwrap_or_default()
        );
        if let Some(&prev) = keys.get(&key) {
            bad(format!("terms[{k}] repeats terms[{prev}] (same name, fragments, spin and parent)"), &mut errors);
        } else {
            keys.insert(key, k);
        }
        if str_of(t.get("name")) == Some("interaction")
            && !t.contains_key("fragments")
            && !t.contains_key("spin")
            && !t.contains_key("parent")
        {
            interactions += 1;
        }
    }
    for (k, t) in terms.iter().enumerate() {
        let mut seen: HashSet<i64> = HashSet::from([k as i64]);
        let mut node = idx(t.get("parent"));
        while let Some(n) = node.filter(|n| *n >= 0 && *n < n_terms) {
            if !seen.insert(n) {
                bad(format!("terms[{k}] is its own ancestor through parent"), &mut errors);
                break;
            }
            node = idx(terms[n as usize].get("parent"));
        }
    }
    if p.get("terms").is_some_and(Value::is_array) && interactions != 1 {
        bad(
            format!("terms has {interactions} interaction terms (name 'interaction' without fragments, spin or parent); exactly one is required"),
            &mut errors,
        );
    }

    // --- units -------------------------------------------------------------------------
    for (where_, e) in energies(p) {
        if let Some(u) = str_of(as_obj(Some(e)).and_then(|o| o.get("units"))) {
            if hartree_per_energy_unit(u).is_none() {
                errors.push((
                    Code::E_UNIT_CONTRADICTED,
                    format!(
                        "{where_} is in '{u}', which is not an energy unit (hartree, ev, kcal_per_mol, kj_per_mol) ({SECTION})"
                    ),
                ));
            }
        }
    }

    // --- sum rules ---------------------------------------------------------------------
    for (k, r) in list(p.get("sum_rules")).iter().enumerate() {
        let Some(r) = r.as_object() else { continue };
        let total = idx(r.get("total"));
        let parts_v = r.get("parts").and_then(Value::as_array);
        let parts: Vec<i64> = parts_v.map(|a| a.iter().filter_map(|x| idx(Some(x))).collect()).unwrap_or_default();
        let all: Vec<i64> = total.into_iter().chain(parts.iter().copied()).collect();
        let out_of_range: Vec<i64> = all.iter().copied().filter(|i| *i < 0 || *i >= n_terms).collect();
        if !out_of_range.is_empty() {
            bad(
                format!("sum_rules[{k}] names term(s) [{}], but there are {n_terms} terms", join_ints(&out_of_range)),
                &mut errors,
            );
            continue;
        }
        let (Some(total), Some(_)) = (total, parts_v) else { continue };
        if parts.contains(&total) {
            bad(format!("sum_rules[{k}] lists its total terms[{total}] among its parts"), &mut errors);
            continue;
        }
        let values: Vec<Option<f64>> = all.iter().map(|&i| energy_eh(terms[i as usize].get("energy"))).collect();
        let tol = if r.contains_key("tolerance") { energy_eh(r.get("tolerance")) } else { Some(DEFAULT_SUM_TOLERANCE_EH) };
        // A negative tolerance is reported by the schema check.
        let (Some(tol), true) = (tol.filter(|t| *t >= 0.0), values.iter().all(Option::is_some)) else { continue };
        let v: Vec<f64> = values.into_iter().flatten().collect();
        let sum: f64 = v[1..].iter().sum();
        let residual = v[0] - sum;
        if residual.abs() > tol {
            errors.push((
                Code::E_SUM_RULE,
                format!(
                    "sum_rules[{k}]: terms[{total}] is {} hartree but its parts [{}] add up to {} hartree (residual {residual:.3e}, tolerance {tol:e}) ({SECTION})",
                    v[0],
                    join_ints(&parts),
                    sum
                ),
            ));
        }
    }

    // --- orbital pairs -----------------------------------------------------------------
    let mut n_k: Option<usize> = None;
    if let Some(kp) = as_obj(p.get("periodic")).and_then(|x| as_obj(x.get("kpoints"))) {
        if let Some(kps) = kp.get("kpoints").and_then(Value::as_array) {
            n_k = Some(kps.len());
            if let Some(w) = kp.get("weights").and_then(Value::as_array) {
                if w.len() != kps.len() {
                    bad(format!("periodic.kpoints has {} k-points but {} weights", kps.len(), w.len()), &mut errors);
                }
            }
        }
    }
    if let Some(pairs) = as_obj(p.get("orbital_pairs")) {
        if str_of(pairs.get("type")) == Some("nocv") && !pairs.contains_key("eigenvalue_convention") {
            bad("orbital_pairs of type 'nocv' must state eigenvalue_convention".into(), &mut errors);
        }
        let mut seen: HashMap<String, usize> = HashMap::new();
        for (k, q) in list(pairs.get("pairs")).iter().enumerate() {
            let Some(q) = q.as_object() else { continue };
            let spin = str_of(q.get("spin")).unwrap_or("restricted");
            let key = format!("{}/{spin}", q.get("index").map(Value::to_string).unwrap_or_default());
            if let Some(&prev) = seen.get(&key) {
                bad(
                    format!("orbital_pairs.pairs[{k}] repeats index {} of spin '{spin}' (pairs[{prev}])", show(q.get("index"))),
                    &mut errors,
                );
            } else {
                seen.insert(key, k);
            }
            for field in ["donor", "acceptor"] {
                if let Some(f) = idx(q.get(field)) {
                    if f < 0 || f >= n_frag {
                        bad(
                            format!("orbital_pairs.pairs[{k}].{field} names fragment {f}, but there are {n_frag} fragments"),
                            &mut errors,
                        );
                    }
                }
            }
            for field in ["kpoint_eigenvalues", "kpoint_energies"] {
                if !q.contains_key(field) {
                    continue;
                }
                match (n_k, q.get(field).and_then(Value::as_array)) {
                    (None, _) => bad(format!("orbital_pairs.pairs[{k}].{field} needs periodic.kpoints"), &mut errors),
                    (Some(n), Some(a)) if a.len() != n => bad(
                        format!("orbital_pairs.pairs[{k}].{field} has {} entries for {} k-points", a.len(), num(n as f64)),
                        &mut errors,
                    ),
                    _ => {}
                }
            }
        }
    }
    errors
}
