use chitose::{ee_terms, ee_terms_log, SubLevel, SubLevelType};
use clap::Parser;

#[derive(Debug, Parser)]
struct Config {
    /// Type of a sublevel (0 for s, 1 for p, etc)
    #[arg(short = 'l')]
    orbital: u8,
    /// Number of electrons
    #[arg(short = 'n')]
    electrons: u8,
    /// If set, prints all of the states
    #[arg(short, default_value_t = false)]
    verbose: bool,
}

pub fn main() {
    let config = Config::parse();
    let level_type = SubLevelType(config.orbital);
    let level = SubLevel::new(level_type, config.electrons).unwrap();
    let terms = if config.verbose {
        ee_terms_log(level, std::io::stdout)
    } else {
        ee_terms(level)
    }
    .unwrap();
    println!("\nFound terms:");
    for term in terms {
        println!("{}", term);
    }
}
