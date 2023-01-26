use mentori::normalize;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::collections::HashMap;
use std::hash::Hash;
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

type Cache = HashMap<(u32, i64, String, String), Vec<(u32, i64, String, String, String)>>;

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
        let mut pos = 0;
        let mut map: Cache = HashMap::new();

        for r in vcf.records() {
            let record = r.unwrap();

            if record.allele_count() > 2 {
                continue;
            }

            let alleles = record.alleles();

            if let Some((reference, alternate)) = alleles.split_first() {
                if pos != record.pos() {
                    print_duplicates(record.header(), &mut map);
                    map.clear();
                }

                let key = (
                    record.rid().unwrap(),
                    record.pos(),
                    String::from_utf8(reference.to_vec()).unwrap(),
                    String::from_utf8(alternate[0].to_vec()).unwrap(),
                );

                if !map.contains_key(&key) {
                    map.insert(key.to_owned(), Vec::new());
                }

                let vec = map.get_mut(&key).unwrap();
                vec.push((
                    record.rid().unwrap(),
                    record.pos(),
                    String::from_utf8(record.id()).unwrap(),
                    String::from_utf8(reference.to_vec()).unwrap(),
                    String::from_utf8(alternate[0].to_vec()).unwrap(),
                ));

                pos = record.pos();
            }
        }

        print_duplicates(vcf.header(), &mut map);
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

fn print_duplicates(header: &bcf::header::HeaderView, map: &mut Cache) {
    let dups = map
        .iter()
        .filter(|(_k, v)| v.len() > 1)
        .collect::<HashMap<_, _>>();

    for dup in dups.values() {
        for r in dup.iter() {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                str::from_utf8(header.rid2name(r.0).unwrap()).unwrap(),
                r.1 + 1,
                r.2,
                r.3,
                r.4
            )
        }
    }
}
