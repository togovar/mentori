use std::result;
use thiserror::Error;

pub type Result<T, E = Error> = result::Result<T, E>;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    RustHtslibError(#[from] rust_htslib::errors::Error),

    #[error(transparent)]
    FromUtf8Error(#[from] std::string::FromUtf8Error),

    #[error(transparent)]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error("Reference ID not found")]
    ReferenceIdNotFoundError,
}
