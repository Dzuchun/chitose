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
        f.write_str(match self.0 {
            0 => "s",
            1 => "p",
            2 => "d",
            3 => "f",
            4 => "g",
            5 => "h",
            6 => "i",
            _ => "(ORBITALS GO BRRRR)",
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
        f.write_str(match self.0 {
            0 => "S",
            1 => "P",
            2 => "D",
            3 => "F",
            4 => "G",
            5 => "H",
            6 => "I",
            _ => "(MOMENTUM GO BRRRR)",
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
    ee_terms_log(l, sink())
}

static SEPARATOR: &[u8] = " ----- \n".as_bytes();
static SPINS: [i8; 2] = [-1, 1]; // SPINS ARE DOUBLED IN THE CODE!!!!!

pub fn ee_terms_log(l: SubLevel, mut log: impl Write) -> Result<Vec<TermType>, std::io::Error> {
    writeln!(log, "Sublevel: {l}")?;
    log.write(SEPARATOR)?;
    writeln!(log, "Single electron states (ml, ms)")?;
    let mut single_states: Vec<(i8, i8)> = Vec::new();
    let mut sn = 1;
    for ml in l.tp.mls() {
        for ms in SPINS {
            writeln!(log, "{}: ({}, {}/2)", sn, ml, ms)?;
            single_states.push((ml, ms));
            sn += 1;
        }
    }
    let single_states_num = sn - 1;

    log.write(SEPARATOR)?;
    writeln!(log, "Level states")?;
    fn backtrack(
        ind: usize,
        iterator: &mut Vec<usize>,
        min: usize,
        max: usize,
        level_states: &mut Vec<(String, i8, i8)>,
        single_states: &[(i8, i8)],
    ) {
        if ind == iterator.len() {
            let state = iterator.iter().fold((String::new(), 0, 0), |acc, &next| {
                (
                    format!("{}{} ", acc.0, next + 1),
                    acc.1 + single_states[next].0,
                    acc.2 + single_states[next].1,
                )
            });
            level_states.push(state);
            return;
        }
        for val in min..max {
            iterator[ind] = val;
            backtrack(
                ind + 1,
                iterator,
                val + 1,
                max + 1,
                level_states,
                single_states,
            );
        }
    }

    let mut level_states = Vec::new();
    backtrack(
        0,
        &mut (0..l.electrons as usize).collect(),
        0,
        single_states_num + 1 - (l.electrons as usize),
        &mut level_states,
        single_states.as_slice(),
    );
    let mut sorted_states: BTreeMap<i8, BTreeMap<i8, Vec<String>>> = BTreeMap::new();
    writeln!(log, "({} total)", level_states.len())?;
    for (name, ml, ms) in level_states {
        writeln!(
            log,
            "{}: ({}, {})",
            name,
            ml,
            if ms & 1 == 0 {
                (ms / 2).to_string()
            } else {
                format!("{}/2", ms)
            }
        )?;

        sorted_states
            .entry(ml)
            .or_default()
            .entry(ms)
            .or_default()
            .push(name);
    }

    log.write(SEPARATOR)?;
    writeln!(log, "Terms:")?;
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
        writeln!(log, "{term}")?;
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
                writeln!(log, "- {this_state}")?;
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
mod tests {}
