//! The CLI app for splitting a secret and reassemble them
use std::{
    error::Error,
    fs::File,
    io::{Read, Write},
    path::PathBuf,
};

use clap::{Parser, Subcommand};
use rand::rngs::OsRng;
use shamirsecretsharing::secretsharing::{SecretShare, SecretSharing256};
use zip::write::FileOptions;

/// owner can read, write, execute; group can read and execute; everyone can read and execute
const DEFAULT_ARCHIVE_UNIX_PERM: u32 = 0o755;

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

        /// Secret shares will be compressed into a single zip archive and written to the specified
        /// file path; default to "shamir.zip"
        #[arg(short = 'o')]
        outpath: Option<String>,

        /// The secret. If a file is supplied, the content of the file will be the secret. If no
        /// file is supplied, use stdin.
        input: Option<String>,
    },
    Assemble {
        /// The paths to the shares
        inputs: Vec<PathBuf>,
    },
}

fn open_file_or_stdin(path: Option<String>) -> Result<Box<dyn Read>, Box<dyn Error>> {
    let reader: Box<dyn Read> = match path {
        Some(path_str) => {
            let file = std::fs::File::open(path_str)?;
            Box::new(file)
        }
        None => Box::new(std::io::stdin()),
    };

    Ok(reader)
}

fn handle_split(
    threshold: usize,
    redundancy: usize,
    outpath: Option<String>,
    input: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut input = open_file_or_stdin(input)?;
    let outfile = File::create(&outpath.unwrap_or("./shamir.zip".to_string()))?;
    let mut archive = zip::ZipWriter::new(outfile);
    let archive_options: FileOptions<()> = FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        .unix_permissions(DEFAULT_ARCHIVE_UNIX_PERM);

    let mut secret = String::new();
    input.read_to_string(&mut secret)?;
    let mut secretsharing = SecretSharing256::init_with_rng(&mut OsRng, threshold);
    secretsharing.encrypt(&secret.as_bytes())?;
    secretsharing.safe_split(redundancy);
    let shares = secretsharing.stringify_shards()?;
    // println!("count of shares {}", shares.len());

    let _ = shares
        .iter()
        .enumerate()
        .map(|(i, share)| {
            let filename = format!("{i}.txt");
            let share_str = share.to_string()?;
            // TODO: add option for verbosity?
            // println!("{share_str}");
            archive.start_file(&filename, archive_options)?;
            archive.write_all(&share_str.as_bytes())?;

            Ok::<(), Box<dyn Error>>(())
        })
        .collect::<Result<Vec<()>, _>>()?;
    archive.finish()?;

    Ok(())
}

fn handle_assemble(inputs: Vec<PathBuf>) -> Result<(), Box<dyn std::error::Error>> {
    let shares = inputs
        .iter()
        .map(|path| {
            let mut file = File::open(path)?;
            let mut file_str = String::new();
            file.read_to_string(&mut file_str)?;
            let share = SecretShare::from_string(&file_str)?;
            Ok(share)
        })
        .collect::<Result<Vec<SecretShare>, Box<dyn Error>>>()?;
    let decryption = SecretSharing256::decrypt_from_secret_shares(&shares)?;

    print!("{}", String::from_utf8(decryption)?);

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    match args.command {
        ShamirSubcommands::Split {
            threshold,
            redundancy,
            outpath,
            input,
        } => handle_split(threshold, redundancy, outpath, input),
        ShamirSubcommands::Assemble { inputs } => handle_assemble(inputs),
    }?;

    Ok(())
}
