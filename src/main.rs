use mentori::normalize;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::path::PathBuf;
use structopt::clap::crate_description;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(about = crate_description!())]
struct Options {
    /// Silence all output
    #[structopt(short = "q", long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[structopt(short = "v", long = "verbose", parse(from_occurrences))]
    verbose: usize,

    /// Path to file to process.
    #[structopt(parse(from_os_str))]
    pub input: Option<PathBuf>,
}

pub fn main() {
    let options = Options::from_args();

    stderrlog::new()
        .module(module_path!())
        .quiet(options.quiet)
        .verbosity(options.verbose)
        .timestamp(stderrlog::Timestamp::Millisecond)
        .init()
        .unwrap();

    let mut vcf = if let Some(path) = options.input {
        bcf::Reader::from_path(path).unwrap()
    } else {
        bcf::Reader::from_stdin().unwrap()
    };

    let header = bcf::Header::from_template(vcf.header());
    let mut output = bcf::Writer::from_stdout(&header, true, bcf::Format::Vcf).unwrap();

    for r in vcf.records() {
        let mut record = r.unwrap();

        if normalize(&mut record).is_ok() {
            output.write(&record).unwrap_or(());
        }
    }
}
