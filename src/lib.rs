pub mod errors;

use errors::Error;
use errors::Result;
use rust_htslib::bcf;
use std::str;

pub fn normalize(record: &mut bcf::Record) -> Result<()> {
    let alleles = record.alleles();

    if let Some(_) = record.rid() {
        if record.allele_count() == 2 {
            if let Some((reference, alternate)) = alleles.split_first() {
                let (shift, reference, alternate) =
                    remove_shared_bases(str::from_utf8(reference)?, str::from_utf8(alternate[0])?);

                let reference = reference.to_owned();
                let alternate = alternate.to_owned();

                let alleles = [reference.as_bytes(), alternate.as_bytes()];
                record.set_pos(record.pos() + shift);
                record.set_alleles(&alleles)?;
            }
        }
    } else {
        Err(Error::ReferenceIdNotFoundError)?;
    }

    Ok(())
}

fn remove_shared_bases<'a>(reference: &'a str, alternate: &'a str) -> (i64, &'a str, &'a str) {
    let mut r_itr = reference.char_indices().rev();
    let mut a_itr = alternate.char_indices().rev();

    let mut r_to: usize = reference.len();
    let mut a_to: usize = alternate.len();

    while let (Some((i1, c1)), Some((i2, c2))) = (r_itr.next(), a_itr.next()) {
        if i1 == 0 || i2 == 0 || c1 != c2 {
            break;
        } else {
            r_to = i1;
            a_to = i2;
        }
    }

    let mut r_itr = reference.char_indices();
    let mut a_itr = alternate.char_indices();

    let mut r_from: usize = 0;
    let mut a_from: usize = 0;

    while let (Some((i1, c1)), Some((i2, c2))) = (r_itr.next(), a_itr.next()) {
        if i1 + 1 == r_to || i2 + 1 == a_to || c1 != c2 {
            break;
        } else {
            r_from = i1 + 1;
            a_from = i2 + 1;
        }
    }

    (
        r_from as i64,
        &reference[r_from..r_to],
        &alternate[a_from..a_to],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remove_shared_bases_1() {
        let (shift, reference, alternate) = remove_shared_bases("T", "A");

        assert_eq!(shift, 0);
        assert_eq!(reference, "T");
        assert_eq!(alternate, "A");
    }

    #[test]
    fn test_remove_shared_bases_2() {
        let (shift, reference, alternate) = remove_shared_bases("TA", "AA");

        assert_eq!(shift, 0);
        assert_eq!(reference, "T");
        assert_eq!(alternate, "A");
    }

    #[test]
    fn test_remove_shared_bases_3() {
        let (shift, reference, alternate) = remove_shared_bases("ACCCTAAC", "A");

        assert_eq!(shift, 0);
        assert_eq!(reference, "ACCCTAAC");
        assert_eq!(alternate, "A");
    }

    #[test]
    fn test_remove_shared_bases_4() {
        let (shift, reference, alternate) = remove_shared_bases("ACCCTAAC", "CCCCTAAC");

        assert_eq!(shift, 0);
        assert_eq!(reference, "A");
        assert_eq!(alternate, "C");
    }

    #[test]
    fn test_remove_shared_bases_5() {
        let (shift, reference, alternate) = remove_shared_bases("ACCCTAAC", "ACCCCTAAC");

        assert_eq!(shift, 0);
        assert_eq!(reference, "A");
        assert_eq!(alternate, "AC");
    }

    #[test]
    fn test_remove_shared_bases_6() {
        let (shift, reference, alternate) = remove_shared_bases("CTA", "CT");

        assert_eq!(shift, 1);
        assert_eq!(reference, "TA");
        assert_eq!(alternate, "T");
    }
}
