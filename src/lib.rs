use itertools::Itertools;
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Display,
    io::{sink, Write},
};

use thiserror::Error;

#[derive(Debug)]
pub struct SubLevelType(pub u8);

static WHY: &str = "Should be able to express as u8 (why would you need sublevel with L=50, lol?)";

impl SubLevelType {
    pub fn max_electrons(&self) -> u8 {
        self.0
            .checked_mul(2)
            .and_then(|r| r.checked_add(1))
            .and_then(|r| r.checked_mul(2))
            .expect(WHY)
    }

    pub fn mls(&self) -> impl IntoIterator<Item = i8> {
        let l: i8 = self.0.try_into().expect(WHY);
        -l..=l
    }
}

impl Display for SubLevelType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("(L={})", self.0);
        f.write_str(match self.0 {
            0 => "s",
            1 => "p",
            2 => "d",
            3 => "f",
            4 => "g",
            5 => "h",
            6 => "i",
            _ => &s,
        })
    }
}

pub struct SubLevel {
    tp: SubLevelType,
    electrons: u8,
}

impl Display for SubLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{}^{{{}}}", self.tp, self.electrons))
    }
}

#[derive(Debug, Error)]
pub enum LevelError {
    #[error("There could be at most {} electrons on the {0} sublevel", .0.max_electrons())]
    ToMuch(SubLevelType),
}

impl SubLevel {
    pub fn new(t: SubLevelType, electrons: u8) -> Result<Self, LevelError> {
        if electrons <= t.max_electrons() {
            Ok(Self { tp: t, electrons })
        } else {
            Err(LevelError::ToMuch(t))
        }
    }
}

#[derive(Debug, Hash, PartialEq, Eq)]
struct TermMomentum(usize);

impl Display for TermMomentum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("(L={})", self.0);
        f.write_str(match self.0 {
            0 => "S",
            1 => "P",
            2 => "D",
            3 => "F",
            4 => "G",
            5 => "H",
            6 => "I",
            _ => &s,
        })
    }
}

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct TermType {
    momentum: TermMomentum,
    multiplet: usize,
}

impl Display for TermType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("^{{{}}}{}", self.multiplet, self.momentum))
    }
}

pub fn ee_terms(l: SubLevel) -> Result<Vec<TermType>, std::io::Error> {
    ee_terms_log(l, sink)
}

static SEPARATOR: &[u8] = " ----- \n".as_bytes();
static SPINS: [i8; 2] = [-1, 1]; // SPINS ARE DOUBLED IN THE CODE!!!!!

pub fn ee_terms_log<W: Write>(
    l: SubLevel,
    log: impl Fn() -> W,
) -> Result<Vec<TermType>, std::io::Error> {
    writeln!(log(), "Sublevel: {l}")?;
    log().write(SEPARATOR)?;

    let single_states =
        l.tp.mls()
            .into_iter()
            .cartesian_product(SPINS)
            .collect_vec();
    let single_states_num = single_states.len();
    writeln!(log(), "Single electron states ({single_states_num} total)")?;
    single_states
        .iter()
        .enumerate()
        .map(|(i, (ml, ms))| writeln!(log(), "{i}: ({ml}, {ms}/2)"))
        .try_collect()?;
    log().write(SEPARATOR)?;

    let level_states = (0..single_states_num)
        .combinations(l.electrons as usize)
        .map(|state| {
            let mut repr = (String::new(), 0, 0);
            state.into_iter().for_each(|next| {
                repr.0.push_str(&format!("{} ", next + 1));
                repr.1 += single_states[next].0;
                repr.2 += single_states[next].1;
            });
            repr
        })
        .collect_vec();
    writeln!(log(), "Level states")?;
    writeln!(log(), "({} total)", level_states.len())?;
    level_states
        .iter()
        .map(|(name, ml, ms)| {
            writeln!(
                log(),
                "{}: ({}, {})",
                name,
                ml,
                if ms & 1 == 0 {
                    (ms / 2).to_string()
                } else {
                    format!("{}/2", ms)
                }
            )
        })
        .try_collect()?;
    log().write(SEPARATOR)?;

    // Here's a fancy approach with itertool's groups, but it ends up with some states lost for some reason :idk:
    /*
    let mut sorted_states: BTreeMap<i8, BTreeMap<i8, Vec<String>>> = level_states
        .into_iter()
        .group_by(|(_, ml, _)| *ml)
        .into_iter()
        .map(|(ml, l_states)| {
            (
                ml,
                l_states
                    .into_iter()
                    .group_by(|(_, _, ms)| *ms)
                    .into_iter()
                    .map(|(ms, ls_states)| {
                        (
                            ms,
                            ls_states.into_iter().map(|(name, _, _)| name).collect_vec(),
                        )
                    })
                    .collect(),
            )
        })
        .collect();
     */
    let mut sorted_states: BTreeMap<i8, BTreeMap<i8, Vec<String>>> = BTreeMap::new();
    level_states.into_iter().for_each(|(name, ml, ms)| {
        sorted_states
            .entry(ml)
            .or_default()
            .entry(ms)
            .or_default()
            .push(name);
    });

    writeln!(log(), "Terms:")?;
    let mut term_states: HashMap<TermType, Vec<String>> = HashMap::new();
    while let Some((&l, _)) = sorted_states.last_key_value() {
        let l_states = sorted_states
            .get(&l)
            .expect("Should be entry for this key, enforced above");
        let (&s, _) = l_states
            .last_key_value()
            .expect("There's at least one state, so entry should exist");
        let term = TermType {
            momentum: TermMomentum(l.try_into().expect("Max momentum must be nonnegative!")),
            multiplet: (s + 1).try_into().expect("Max spin must be nonnegative!"), // SPIN IS DOUBLED!
        };
        writeln!(log(), "{term}")?;
        let this_term_states = term_states.entry(term).or_default();

        for l in -l..=l {
            let l_states = sorted_states
                .get_mut(&l)
                .expect("Should be entry for this key, enforced above");
            for s in (-s..=s).step_by(2) {
                let sl_states = l_states
                    .get_mut(&s)
                    .expect("Should be entry states with this spin!");
                let this_state = sl_states
                    .pop()
                    .expect("Should be at least one state, will be enforced now");
                writeln!(log(), "- {this_state}")?;
                this_term_states.push(this_state);
                if sl_states.is_empty() {
                    l_states.remove_entry(&s);
                }
            }

            if l_states.is_empty() {
                sorted_states.remove_entry(&l);
            }
        }
    }

    Ok(term_states.into_keys().collect())
}

#[cfg(test)]
mod tests {
    use crate::ee_terms_log;

    #[test]
    fn it_works() {
        ee_terms_log(
            crate::SubLevel {
                tp: crate::SubLevelType(1),
                electrons: 3,
            },
            std::io::stdout,
        )
        .expect("Should be ok");
    }
}
