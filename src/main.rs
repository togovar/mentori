use mentori::normalize;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::path::PathBuf;
use std::str;
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

    /// Check if the VCF has duplicated record
    #[structopt(short = "c", long = "check")]
    check: bool,

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

    if options.check {
        let mut buf = Vec::new();
        let mut pos = 0;

        for r in vcf.records() {
            let record = r.unwrap();

            if pos != record.pos() {
                buf.clear();
            }

            let alleles = record.alleles();
            if let Some((reference, alternate)) = alleles.split_first() {
                let keys = (
                    record.rid().unwrap(),
                    String::from_utf8(reference.to_vec()).unwrap(),
                    String::from_utf8(alternate[0].to_vec()).unwrap(),
                );

                if buf.contains(&keys) {
                    let name = record.header().rid2name(record.rid().unwrap()).unwrap();

                    println!(
                        "VCF has duplicated record: {}\t{}\t{}\t{}\t{}",
                        str::from_utf8(name).unwrap(),
                        record.pos() + 1,
                        str::from_utf8(record.id().as_slice()).unwrap(),
                        str::from_utf8(reference).unwrap(),
                        str::from_utf8(alternate[0]).unwrap()
                    )
                }

                buf.push(keys);
            }

            pos = record.pos();
        }
    } else {
        let header = bcf::Header::from_template(vcf.header());
        let mut output = bcf::Writer::from_stdout(&header, true, bcf::Format::Vcf).unwrap();

        for r in vcf.records() {
            let mut record = r.unwrap();

            if normalize(&mut record).is_ok() {
                output.write(&record).unwrap_or(());
            }
        }
    }
}
