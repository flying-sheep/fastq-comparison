use std::io;
use std::io::prelude::*;
use std::ascii::AsciiExt;
use std::fs;
use std::fmt;
use std::path::Path;
use std::convert::AsRef;

use super::Record as RecordTrait;


/// A FastQ reader.
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    sep_line: String,
}


impl Reader<fs::File> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}


impl<R: io::Read> Reader<R> {
    /// Read from a given `io::Read`.
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            sep_line: String::new(),
        }
    }

    /// Read into a given record.
    /// Returns an error if the record in incomplete or syntax is violated.
    /// The content of the record can be checked via the record object.
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        try!(self.reader.read_line(&mut record.header));

        if !record.header.is_empty() {
            if !record.header.starts_with('@') {
                return Err(io::Error::new(io::ErrorKind::Other, "Expected @ at record start."));
            }
            try!(self.reader.read_line(&mut record.seq));
            try!(self.reader.read_line(&mut self.sep_line));
            try!(self.reader.read_line(&mut record.qual));
            if record.qual.is_empty() {
                return Err(io::Error::new(io::ErrorKind::Other,
                                          "Incomplete record. Each FastQ record has to consist \
                                           of 4 lines: header, sequence, separator and qualities."));
            }
        }

        Ok(())
    }

    /// Return an iterator over the records of this FastQ file.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}


/// A FastQ record.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    header: String,
    seq: String,
    qual: String,
}


impl super::Record for Record {
    /// Create a new, empty FastQ record.
    fn new() -> Self {
        Record {
            header: String::new(),
            seq: String::new(),
            qual: String::new(),
        }
    }

    /// Check if record is empty.
    fn is_empty(&self) -> bool {
        self.header.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    /// Check validity of FastQ record.
    fn check(&self) -> Result<(), &str> {
        if self.id().is_none() {
            return Err("Expecting id for FastQ record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }
        if !self.qual.is_ascii() {
            return Err("Non-ascii character found in qualities.");
        }
        if self.seq().len() != self.qual().len() {
            return Err("Unequal length of sequence an qualities.");
        }

        Ok(())
    }

    /// Return the id of the record.
    fn id(&self) -> Option<&str> {
        self.header[1..].trim_right().splitn(2, ' ').next()
    }

    /// Return descriptions if present.
    fn desc(&self) -> Option<&str> {
        self.header[1..].trim_right().splitn(2, ' ').skip(1).next()
    }

    /// Return the sequence of the record.
    fn seq(&self) -> &[u8] {
        self.seq.trim_right().as_bytes()
    }

    /// Return the base qualities of the record.
    fn qual(&self) -> &[u8] {
        self.qual.trim_right().as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.header.clear();
        self.seq.clear();
        self.qual.clear();
    }
}


impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "@{}\n{}\n+\n{}", self.header, self.seq, self.qual)
    }
}


/// An iterator over the records of a FastQ file.
pub struct Records<R: io::Read> {
    reader: Reader<R>,
}


impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}
