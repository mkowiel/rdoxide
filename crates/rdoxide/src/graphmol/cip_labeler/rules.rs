//! The CIP sequence rules and the sort that applies them (upstream
//! `CIPLabeler/{Sort.{h,cpp},Priority.h,rules/*}`).
//!
//! Upstream rules are objects, each owning a `Sort` over the rules up to and
//! including itself. Here a rule is a small value ([`RuleRef`]) that names
//! its kind and which sorter it owns, and the rules are applied by the
//! functions of this module with the digraph passed along.

use super::budget;
use super::descriptor::Descriptor;
use super::digraph::{Digraph, EdgeId, NodeId};
use crate::{RDKitError, Result};
use std::cmp::Ordering;

/// The pseudo-asymmetric part of a sort result (`Priority`)
#[derive(Debug, Clone, Copy)]
pub struct Priority {
    unique: bool,
    pseudo_asym: bool,
}

impl Priority {
    pub fn is_unique(&self) -> bool {
        self.unique
    }

    pub fn is_pseudo_asymmetric(&self) -> bool {
        self.pseudo_asym
    }
}

fn three_way<T: PartialOrd>(x: T, y: T) -> i32 {
    if x < y {
        -1
    } else if x == y {
        0
    } else {
        1
    }
}

// --------------------------------------------------------------- pair lists

/// `PairList`: the descriptors met in a hierarchical digraph, relative to a
/// reference.
#[derive(Debug, Clone, Default)]
struct PairList {
    descriptors: Vec<Descriptor>,
}

impl PairList {
    /// `PairList::ref`
    fn reference(d: Descriptor) -> Descriptor {
        match d {
            Descriptor::R | Descriptor::M | Descriptor::SeqCis => Descriptor::R,
            Descriptor::S | Descriptor::P | Descriptor::SeqTrans => Descriptor::S,
            _ => Descriptor::None,
        }
    }

    fn new(reference: Descriptor) -> Self {
        let mut list = PairList::default();
        list.add(reference);
        list
    }

    fn ref_descriptor(&self) -> Result<Descriptor> {
        match self.descriptors.first() {
            Some(&d) => Ok(Self::reference(d)),
            None => Err(RDKitError::Runtime(
                "Cannot get a reference from an empty PairList".into(),
            )),
        }
    }

    fn add(&mut self, d: Descriptor) -> bool {
        match d {
            Descriptor::R
            | Descriptor::S
            | Descriptor::M
            | Descriptor::P
            | Descriptor::SeqTrans
            | Descriptor::SeqCis => {
                self.descriptors.push(Self::reference(d));
                true
            }
            _ => false,
        }
    }

    fn compare_to(&self, that: &PairList) -> Result<i32> {
        if self.descriptors.len() != that.descriptors.len() {
            return Err(RDKitError::Runtime(
                "Descriptor lists should be the same length!".into(),
            ));
        }
        if self.descriptors.is_empty() {
            return Ok(0);
        }
        let this_ref = self.descriptors[0];
        let that_ref = that.descriptors[0];
        for i in 1..self.descriptors.len() {
            let this_like = this_ref == self.descriptors[i];
            let that_like = that_ref == that.descriptors[i];
            if this_like && !that_like {
                return Ok(1);
            }
            if !this_like && that_like {
                return Ok(-1);
            }
        }
        Ok(0)
    }
}

// -------------------------------------------------------------------- rules

/// The kind of a sequence rule. Rules 4b and 5 carry the reference descriptor
/// they compare against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RuleKind {
    Rule1a,
    Rule1b,
    Rule2,
    Rule3,
    Rule4a,
    Rule4b(Descriptor),
    Rule4c,
    Rule5New(Descriptor),
    Rule6,
}

/// A rule together with the sorter it owns: the rules of its set up to and
/// including itself (`Some(position)`), or itself alone (`None`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RuleRef {
    kind: RuleKind,
    prefix: Option<usize>,
}

/// `Rules`: an ordered set of rules.
#[derive(Debug)]
pub struct Rules {
    refs: Vec<RuleRef>,
}

impl Rules {
    pub fn new(kinds: &[RuleKind]) -> Self {
        Rules {
            refs: kinds
                .iter()
                .enumerate()
                .map(|(i, &kind)| RuleRef {
                    kind,
                    prefix: Some(i),
                })
                .collect(),
        }
    }

    /// The rules up to rule 2, which only look at the constitution
    pub fn constitutional() -> Self {
        Self::new(&[RuleKind::Rule1a, RuleKind::Rule1b, RuleKind::Rule2])
    }

    /// All the rules
    pub fn all() -> Self {
        Self::new(&[
            RuleKind::Rule1a,
            RuleKind::Rule1b,
            RuleKind::Rule2,
            RuleKind::Rule3,
            RuleKind::Rule4a,
            RuleKind::Rule4b(Descriptor::None),
            RuleKind::Rule4c,
            RuleKind::Rule5New(Descriptor::None),
            RuleKind::Rule6,
        ])
    }

    pub fn num_sub_rules(&self) -> usize {
        self.refs.len()
    }

    /// `Rules::sort`: orders the edges of a node with all the rules.
    pub fn sort(
        &self,
        d: &mut Digraph<'_>,
        node: NodeId,
        edges: &mut [EdgeId],
    ) -> Result<Priority> {
        prioritize(self, d, &self.refs, node, edges, true)
    }

    /// `Rules::getSorter()->getGroups`
    pub fn groups(&self, d: &mut Digraph<'_>, sorted: &[EdgeId]) -> Result<Vec<Vec<EdgeId>>> {
        get_groups(self, d, &self.refs, sorted)
    }

    /// The rules of the sorter that belongs to a rule
    fn sorter_of<'s>(&'s self, rule: &'s RuleRef) -> &'s [RuleRef] {
        match rule.prefix {
            Some(k) => &self.refs[..=k],
            None => std::slice::from_ref(rule),
        }
    }
}

// --------------------------------------------------------------------- sort

/// `Sort::prioritize`: insertion sort of the edges, highest priority first.
fn prioritize(
    set: &Rules,
    d: &mut Digraph<'_>,
    rules: &[RuleRef],
    node: NodeId,
    edges: &mut [EdgeId],
    deep: bool,
) -> Result<Priority> {
    let mut unique = true;
    let mut num_pseudo_asym = 0;
    for i in 0..edges.len() {
        for j in (1..=i).rev() {
            let cmp = compare_substituents(set, d, rules, node, edges[j - 1], edges[j], deep)?;
            if !(-1..=1).contains(&cmp) {
                num_pseudo_asym += 1;
            }
            if cmp < 0 {
                edges.swap(j, j - 1);
            } else {
                if cmp == 0 {
                    unique = false;
                }
                break;
            }
        }
    }
    Ok(Priority {
        unique,
        pseudo_asym: num_pseudo_asym == 1,
    })
}

fn compare_substituents(
    set: &Rules,
    d: &mut Digraph<'_>,
    rules: &[RuleRef],
    node: NodeId,
    a: EdgeId,
    b: EdgeId,
    deep: bool,
) -> Result<i32> {
    let a_beg = d.edge(a).is_beg(node);
    let b_beg = d.edge(b).is_beg(node);
    if !a_beg && b_beg {
        return Ok(1);
    } else if a_beg && !b_beg {
        return Ok(-1);
    }
    for rule in rules {
        let cmp = if deep {
            recursive_compare(set, d, rule, a, b)?
        } else {
            compare(set, d, rule, a, b)?
        };
        if cmp != 0 {
            return Ok(cmp);
        }
    }
    Ok(0)
}

/// `Sort::getGroups`: splits sorted edges in groups of equivalent edges.
fn get_groups(
    set: &Rules,
    d: &mut Digraph<'_>,
    rules: &[RuleRef],
    sorted: &[EdgeId],
) -> Result<Vec<Vec<EdgeId>>> {
    let mut groups: Vec<Vec<EdgeId>> = Vec::new();
    let mut prev: Option<EdgeId> = None;
    for &edge in sorted {
        let new_group = match prev {
            None => true,
            Some(p) => {
                let beg = d.edge(p).beg;
                compare_substituents(set, d, rules, beg, p, edge, true)? != 0
            }
        };
        if new_group {
            groups.push(Vec::new());
        }
        prev = Some(edge);
        groups.last_mut().unwrap().push(edge);
    }
    Ok(groups)
}

// ------------------------------------------------------ recursive comparison

/// `SequenceRule::recursiveCompare`
fn recursive_compare(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    a: EdgeId,
    b: EdgeId,
) -> Result<i32> {
    if !budget::consume() {
        return Err(budget::max_iterations_exceeded());
    }
    let cmp = compare(set, d, rule, a, b)?;
    if cmp != 0 {
        return Ok(cmp);
    }
    let sorter = set.sorter_of(rule);
    let mut a_queue = vec![a];
    let mut b_queue = vec![b];
    let mut pos = 0;
    while pos < a_queue.len() && pos < b_queue.len() {
        let a = a_queue[pos];
        let b = b_queue[pos];
        pos += 1;
        let a_node = d.edge(a).end;
        let b_node = d.edge(b).end;
        let mut as_ = d.edges(a_node)?;
        let mut bs = d.edges(b_node)?;
        prioritize(set, d, sorter, a_node, &mut as_, false)?;
        prioritize(set, d, sorter, b_node, &mut bs, false)?;
        let size_diff = three_way(as_.len() as i32, bs.len() as i32);
        for (&a_edge, &b_edge) in as_.iter().zip(&bs) {
            if are_up_edges(d, a_node, b_node, a_edge, b_edge)? {
                continue;
            }
            let cmp = compare(set, d, rule, a_edge, b_edge)?;
            if cmp != 0 {
                return Ok(cmp);
            }
        }
        if size_diff != 0 {
            return Ok(size_diff);
        }
        prioritize(set, d, sorter, a_node, &mut as_, true)?;
        prioritize(set, d, sorter, b_node, &mut bs, true)?;
        for (&a_edge, &b_edge) in as_.iter().zip(&bs) {
            if are_up_edges(d, a_node, b_node, a_edge, b_edge)? {
                continue;
            }
            let cmp = compare(set, d, rule, a_edge, b_edge)?;
            if cmp != 0 {
                return Ok(cmp);
            }
            a_queue.push(a_edge);
            b_queue.push(b_edge);
        }
    }
    Ok(0)
}

fn are_up_edges(
    d: &Digraph<'_>,
    a_node: NodeId,
    b_node: NodeId,
    a_edge: EdgeId,
    b_edge: EdgeId,
) -> Result<bool> {
    if d.edge(a_edge).is_end(a_node) {
        if !d.edge(b_edge).is_end(b_node) {
            return Err(RDKitError::Runtime("Something unexpected!".into()));
        }
        return Ok(true);
    }
    Ok(false)
}

/// The label of the bond of an edge (`SequenceRule::getBondLabel`)
fn bond_label(d: &Digraph<'_>, edge: EdgeId) -> Descriptor {
    let e = d.edge(edge);
    if e.bond.is_none() {
        return Descriptor::None;
    }
    e.aux
}

/// `SequenceRule::compare` of the rule
fn compare(set: &Rules, d: &mut Digraph<'_>, rule: &RuleRef, a: EdgeId, b: EdgeId) -> Result<i32> {
    match rule.kind {
        RuleKind::Rule1a => Ok(three_way(d.end_of(a).atomic_num, d.end_of(b).atomic_num)),
        RuleKind::Rule1b => Ok(compare_rule_1b(d, a, b)),
        RuleKind::Rule2 => Ok(compare_rule_2(d, a, b)),
        RuleKind::Rule3 => Ok(three_way(
            ord_rule_3(d.end_of(a).aux),
            ord_rule_3(d.end_of(b).aux),
        )),
        RuleKind::Rule4a => Ok(compare_ordinals(d, a, b, ord_rule_4a)),
        RuleKind::Rule4c => Ok(compare_ordinals(d, a, b, ord_rule_4c)),
        RuleKind::Rule6 => Ok(compare_rule_6(d, a, b)),
        RuleKind::Rule4b(reference) => compare_rule_4b(set, d, rule, reference, a, b),
        RuleKind::Rule5New(reference) => compare_rule_5(set, d, rule, reference, a, b),
    }
}

/// Compares the bond labels and then the labels of the end nodes.
fn compare_ordinals(d: &Digraph<'_>, a: EdgeId, b: EdgeId, ord: impl Fn(Descriptor) -> i32) -> i32 {
    let cmp = three_way(ord(bond_label(d, a)), ord(bond_label(d, b)));
    if cmp != 0 {
        return cmp;
    }
    three_way(ord(d.end_of(a).aux), ord(d.end_of(b).aux))
}

/// Rule 1b, with the 2013 IUPAC recommendations off
fn compare_rule_1b(d: &Digraph<'_>, a: EdgeId, b: EdgeId) -> i32 {
    use super::digraph::RING_DUPLICATE;
    let (a_end, b_end) = (d.end_of(a), d.end_of(b));
    let (a_ring, b_ring) = (a_end.is_set(RING_DUPLICATE), b_end.is_set(RING_DUPLICATE));
    if a_ring && b_ring {
        -three_way(a_end.dist, b_end.dist)
    } else if a_ring {
        1
    } else if b_ring {
        -1
    } else {
        0
    }
}

fn compare_rule_2(d: &Digraph<'_>, a: EdgeId, b: EdgeId) -> i32 {
    let mol = d.mol();
    let (a_end, b_end) = (d.end_of(a), d.end_of(b));
    let a_num = a_end.atomic_number(mol);
    let b_num = b_end.atomic_number(mol);
    if a_num == 0 && b_num == 0 {
        return 0;
    } else if a_num == 0 || b_num == 0 {
        return three_way(a_num, b_num);
    }
    if a_end.mass_num(mol) == 0 && b_end.mass_num(mol) == 0 {
        return 0;
    }
    three_way(a_end.atomic_mass, b_end.atomic_mass)
}

fn ord_rule_3(lab: Descriptor) -> i32 {
    match lab {
        Descriptor::E => 1,
        Descriptor::Z => 2,
        _ => 0,
    }
}

fn ord_rule_4a(lab: Descriptor) -> i32 {
    match lab {
        Descriptor::Unknown | Descriptor::Ns | Descriptor::None => 0,
        Descriptor::PseudoR
        | Descriptor::PseudoS
        | Descriptor::PseudoM
        | Descriptor::PseudoP
        | Descriptor::E
        | Descriptor::Z => 1,
        // R, S, M, P, seqTrans and seqCis (anything else is no atom or bond
        // descriptor and never gets here)
        _ => 2,
    }
}

fn ord_rule_4c(lab: Descriptor) -> i32 {
    match lab {
        Descriptor::PseudoM | Descriptor::PseudoR => 2,
        Descriptor::PseudoP | Descriptor::PseudoS => 1,
        _ => 0,
    }
}

fn compare_rule_6(d: &Digraph<'_>, a: EdgeId, b: EdgeId) -> i32 {
    let Some(reference) = d.rule6_ref() else {
        return 0;
    };
    let a_atom = d.end_of(a).atom;
    let b_atom = d.end_of(b).atom;
    if Some(reference) == a_atom && Some(reference) != b_atom {
        2 // a is the reference (it has priority and rule 6 follows rule 5)
    } else if Some(reference) != a_atom && Some(reference) == b_atom {
        -2
    } else {
        0
    }
}

// ---------------------------------------------------------- rules 4b and 5

/// `getRefSorter`: the sorter of a rule with the rule itself replaced by
/// `replacement`.
fn ref_sorter(set: &Rules, rule: &RuleRef, replacement: RuleRef) -> Result<Vec<RuleRef>> {
    let rules = set.sorter_of(rule);
    if !rules.contains(rule) {
        return Err(RDKitError::invariant("rule not in rule set"));
    }
    let mut new_rules: Vec<RuleRef> = rules.iter().filter(|r| *r != rule).copied().collect();
    new_rules.push(replacement);
    Ok(new_rules)
}

/// The "like"/"unlike" comparison of rules 4b and 5 for nodes that are not
/// below the root: the edge to the node that is like the reference wins.
fn compare_to_reference(d: &Digraph<'_>, reference: Descriptor, a: EdgeId, b: EdgeId) -> i32 {
    if reference == Descriptor::None {
        return 0;
    }
    let a_desc = d.end_of(a).aux;
    let b_desc = d.end_of(b).aux;
    if a_desc != Descriptor::None
        && b_desc != Descriptor::None
        && a_desc != Descriptor::Ns
        && b_desc != Descriptor::Ns
    {
        let like_a = PairList::reference(reference) == PairList::reference(a_desc);
        let like_b = PairList::reference(reference) == PairList::reference(b_desc);
        if like_a && !like_b {
            return 1;
        }
        if like_b && !like_a {
            return -1;
        }
    }
    0
}

/// Breadth first collection of the descriptors below a node, in the order
/// given by the rules with the reference of the list.
fn fill_pairs(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    beg: NodeId,
    list: &mut PairList,
    replacement: impl Fn(Descriptor) -> RuleKind,
) -> Result<()> {
    let replacement_rule = RuleRef {
        kind: replacement(list.ref_descriptor()?),
        prefix: None,
    };
    let sorter = ref_sorter(set, rule, replacement_rule)?;
    let mut queue = vec![beg];
    let mut i = 0;
    while i < queue.len() {
        let node = queue[i];
        i += 1;
        list.add(d.node(node).aux);
        let mut edges = d.edges(node)?;
        prioritize(set, d, &sorter, node, &mut edges, true)?;
        for edge in edges {
            if d.edge(edge).is_beg(node) && !d.end_of(edge).is_terminal() {
                queue.push(d.edge(edge).end);
            }
        }
    }
    Ok(())
}

fn compare_rule_5(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    reference: Descriptor,
    a: EdgeId,
    b: EdgeId,
) -> Result<i32> {
    let (a_beg, a_end) = (d.edge(a).beg, d.edge(a).end);
    let (b_beg, b_end) = (d.edge(b).beg, d.edge(b).end);
    if d.current_root() != a_beg || d.current_root() != b_beg {
        return Ok(compare_to_reference(d, reference, a, b));
    }
    let mut list_ra = PairList::new(Descriptor::R);
    let mut list_rb = PairList::new(Descriptor::R);
    let mut list_sa = PairList::new(Descriptor::S);
    let mut list_sb = PairList::new(Descriptor::S);
    let make = RuleKind::Rule5New;
    fill_pairs(set, d, rule, a_end, &mut list_ra, make)?;
    fill_pairs(set, d, rule, a_end, &mut list_sa, make)?;
    fill_pairs(set, d, rule, b_end, &mut list_rb, make)?;
    fill_pairs(set, d, rule, b_end, &mut list_sb, make)?;
    let cmp_r = list_ra.compare_to(&list_rb)?;
    let cmp_s = list_sa.compare_to(&list_sb)?;
    Ok(match cmp_r.cmp(&0) {
        Ordering::Less => {
            if cmp_s < 0 {
                -1
            } else {
                -2
            }
        }
        Ordering::Greater => {
            if cmp_s > 0 {
                1
            } else {
                2
            }
        }
        Ordering::Equal => 0,
    })
}

fn compare_rule_4b(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    reference: Descriptor,
    a: EdgeId,
    b: EdgeId,
) -> Result<i32> {
    let (a_beg, a_end) = (d.edge(a).beg, d.edge(a).end);
    let (b_beg, b_end) = (d.edge(b).beg, d.edge(b).end);
    if d.current_root() != a_beg || d.current_root() != b_beg {
        return Ok(compare_to_reference(d, reference, a, b));
    }
    let mut list1 = new_pair_lists(&reference_descriptors(set, d, rule, a_end)?);
    let mut list2 = new_pair_lists(&reference_descriptors(set, d, rule, b_end)?);
    if list1.is_empty() != list2.is_empty() {
        return Err(RDKitError::Runtime(
            "Substituents should be topologically equivalent!".into(),
        ));
    }
    if list1.len() == 1 {
        let ref_a = list1[0].ref_descriptor()?;
        let ref_b = match list2.first() {
            Some(list) => list.ref_descriptor()?,
            None => {
                return Err(RDKitError::Runtime(
                    "Substituents should be topologically equivalent!".into(),
                ));
            }
        };
        return compare_pairs(set, d, rule, a_end, b_end, ref_a, ref_b);
    } else if list1.len() > 1 {
        for list in &mut list1 {
            fill_pairs(set, d, rule, a_end, list, RuleKind::Rule4b)?;
        }
        for list in &mut list2 {
            fill_pairs(set, d, rule, b_end, list, RuleKind::Rule4b)?;
        }
        sort_descending(&mut list1)?;
        sort_descending(&mut list2)?;
        for (l1, l2) in list1.iter().zip(&list2) {
            let cmp = l1.compare_to(l2)?;
            if cmp != 0 {
                return Ok(cmp);
            }
        }
    }
    Ok(0)
}

/// `std::sort(list.rbegin(), list.rend())`, i.e. the greatest list first
fn sort_descending(lists: &mut [PairList]) -> Result<()> {
    let mut error = None;
    lists.sort_by(|x, y| match y.compare_to(x) {
        Ok(-1) => Ordering::Less,
        Ok(1) => Ordering::Greater,
        Ok(_) => Ordering::Equal,
        Err(e) => {
            error.get_or_insert(e);
            Ordering::Equal
        }
    });
    error.map_or(Ok(()), Err)
}

fn new_pair_lists(descriptors: &[Descriptor]) -> Vec<PairList> {
    descriptors.iter().map(|&d| PairList::new(d)).collect()
}

/// `Rule4b::getReferenceDescriptors`: looking level by level into the digraph
/// for the first descriptors, the most common one is the reference (both, if
/// they are equally common).
fn reference_descriptors(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    node: NodeId,
) -> Result<Vec<Descriptor>> {
    let mut prev: Vec<Vec<NodeId>> = vec![vec![node]];
    while !prev.is_empty() {
        for nodes in &prev {
            let mut result = Vec::new();
            if reference_of_level(d, nodes, &mut result) {
                return Ok(result);
            }
        }
        prev = next_level(set, d, rule, &prev)?;
    }
    Ok(Vec::new())
}

fn reference_of_level(d: &Digraph<'_>, nodes: &[NodeId], result: &mut Vec<Descriptor>) -> bool {
    let (mut right, mut left) = (0, 0);
    for &node in nodes {
        match d.node(node).aux {
            Descriptor::R | Descriptor::M | Descriptor::SeqCis => right += 1,
            Descriptor::S | Descriptor::P | Descriptor::SeqTrans => left += 1,
            _ => {}
        }
    }
    if right + left == 0 {
        return false;
    }
    match right.cmp(&left) {
        Ordering::Greater => result.push(Descriptor::R),
        Ordering::Less => result.push(Descriptor::S),
        Ordering::Equal => {
            result.push(Descriptor::R);
            result.push(Descriptor::S);
        }
    }
    true
}

fn next_level(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    prev_level: &[Vec<NodeId>],
) -> Result<Vec<Vec<NodeId>>> {
    let sorter = set.sorter_of(rule);
    let mut next = Vec::with_capacity(4 * prev_level.len());
    for prev in prev_level {
        let mut tmp: Vec<Vec<Vec<EdgeId>>> = Vec::new();
        for &node in prev {
            let mut edges = d.non_terminal_out_edges(node)?;
            prioritize(set, d, sorter, node, &mut edges, true)?;
            tmp.push(get_groups(set, d, sorter, &edges)?);
        }
        let size = tmp.first().map_or(0, Vec::len);
        if tmp.iter().any(|groups| groups.len() != size) {
            return Err(RDKitError::Runtime("Something unexpected!".into()));
        }
        for i in 0..size {
            let eq: Vec<NodeId> = tmp
                .iter()
                .flat_map(|groups| groups[i].iter().map(|&e| d.edge(e).end))
                .collect();
            if !eq.is_empty() {
                next.push(eq);
            }
        }
    }
    Ok(next)
}

/// `Rule4b::comparePairs`
fn compare_pairs(
    set: &Rules,
    d: &mut Digraph<'_>,
    rule: &RuleRef,
    a: NodeId,
    b: NodeId,
    ref_a: Descriptor,
    ref_b: Descriptor,
) -> Result<i32> {
    let sorter_a = ref_sorter(
        set,
        rule,
        RuleRef {
            kind: RuleKind::Rule4b(ref_a),
            prefix: None,
        },
    )?;
    let sorter_b = ref_sorter(
        set,
        rule,
        RuleRef {
            kind: RuleKind::Rule4b(ref_b),
            prefix: None,
        },
    )?;
    let mut a_queue = vec![a];
    let mut b_queue = vec![b];
    let mut pos = 0;
    while pos < a_queue.len() && pos < b_queue.len() {
        let a_node = a_queue[pos];
        let b_node = b_queue[pos];
        pos += 1;
        let des_a = PairList::reference(d.node(a_node).aux);
        let des_b = PairList::reference(d.node(b_node).aux);
        if des_a == ref_a && des_b != ref_b {
            return Ok(1);
        } else if des_a != ref_a && des_b == ref_b {
            return Ok(-1);
        }
        for (node, sorter, queue) in [
            (a_node, &sorter_a, &mut a_queue),
            (b_node, &sorter_b, &mut b_queue),
        ] {
            let mut edges = d.edges(node)?;
            prioritize(set, d, sorter, node, &mut edges, true)?;
            for edge in edges {
                if d.edge(edge).is_beg(node) && !d.end_of(edge).is_terminal() {
                    queue.push(d.edge(edge).end);
                }
            }
        }
    }
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_lists_ignore_what_is_not_a_reference() {
        let mut list = PairList::default();
        for d in [Descriptor::None, Descriptor::Unknown, Descriptor::Ns] {
            assert!(!list.add(d));
        }
        assert!(!list.add(Descriptor::PseudoR));
        assert!(!list.add(Descriptor::PseudoS));
        assert!(list.add(Descriptor::R));
        assert!(list.add(Descriptor::S));
    }

    #[test]
    fn pair_lists_compare_like_and_unlike() {
        let mut list1 = PairList::default();
        let mut list2 = PairList::default();
        for d in [Descriptor::R, Descriptor::M, Descriptor::R, Descriptor::S] {
            list1.add(d);
        }
        for d in [Descriptor::R, Descriptor::P, Descriptor::S, Descriptor::M] {
            list2.add(d);
        }
        // R:llu against R:uul
        assert_eq!(list1.compare_to(&list2).unwrap(), 1);
        assert_eq!(list2.compare_to(&list1).unwrap(), -1);
        assert!(list1.compare_to(&PairList::new(Descriptor::R)).is_err());
    }

    #[test]
    fn rule_sets_own_prefix_sorters() {
        let rules = Rules::all();
        assert_eq!(rules.num_sub_rules(), 9);
        assert_eq!(Rules::constitutional().num_sub_rules(), 3);
        let rule4b = rules.refs[5];
        assert_eq!(rules.sorter_of(&rule4b).len(), 6);
        let alone = RuleRef {
            kind: RuleKind::Rule4b(Descriptor::R),
            prefix: None,
        };
        assert_eq!(rules.sorter_of(&alone), &[alone]);
    }
}
