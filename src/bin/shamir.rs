//! The CLI app for splitting a secret and reassemble them
use clap::{Parser, Subcommand};

/// Splits secret or reassemble them
#[derive(Debug, Parser)]
#[command(version = "0.1")]
pub struct Args {
    #[command(subcommand)]
    command: ShamirSubcommands,
}

#[derive(Debug, Subcommand)]
pub enum ShamirSubcommands {
    Split {
        /// The number of shares needed to assemble the secret
        #[arg(short = 't')]
        threshold: usize,

        /// The number of shares to split the secret into. This number must be no less than the
        /// threshold
        #[arg(short = 'n')]
        redundancy: usize,

        /// The directory to write the shares to
        #[arg(short = 'o')]
        outdir: String,

        /// The secret. If a file is supplied, the content of the file will be the secret. If no
        /// file is supplied, use stdin.
        input: Option<String>,
    },
    Assemble {
        /// The paths to the shares
        input: Vec<String>,
    },
}

fn main() {
    let _args = Args::parse();
}
