use std::ascii::AsciiExt;
use std::error::Error;
use std::io::{self,Write,BufRead};

pub struct Record {
	id: String,
	desc: Option<String>,
	seq: String,
	qual: String,
}

impl Record {
	pub fn from_strings(id: String, desc: Option<String>, seq: String, qual: String) -> Record {
		Record { id: id, desc: desc, seq: seq, qual: qual }
	}
}

impl super::Record for Record {
	fn new() -> Record {
		Record { id: String::new(), desc: None, seq: String::new(), qual: String::new() }
	}
	
	fn id(&self) -> Option<&str> { Some(self.id.as_ref()) }
	fn desc(&self) -> Option<&str> { self.desc.as_ref().map(|s| s.as_str()) }
	fn seq(&self) -> &[u8] { self.seq.as_bytes() }
	fn qual(&self) -> &[u8] { self.qual.as_bytes() }
	
	fn check(&self) -> Result<(), &'static str> {
		if !self.id.is_empty() {
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
	
	fn is_empty(&self) -> bool {
		self.id.is_empty() && self.desc.is_none() && self.seq.is_empty() && self.qual.is_empty()
	}
	
	fn clear(&mut self) {
		self.id.clear();
		self.desc = None;
		self.seq.clear();
		self.qual.clear();
	}
}

quick_error!(
	#[derive(Debug)]
	pub enum ParseError {
		NoAt(byte: u8) {
			description("No @ at FASTQ start")
			display("Encountered {:?} instead of @", byte)
		}
		NoPlus(prev: String, byte: u8) {
			description("No + after FASTQ sequence")
			display("Encountered {:?} instead of + after {:?}", byte, prev)
		}
		Incomplete(prev: String) {
			description("Incomplete FASTQ record")
			display("Premature EOF after {:?}", prev)
		}
		LengthMismatch(a: String, b: String) {
			description("FASTQ with differing length")
			display("The lengths do not match:\nseq:  {:?}\nqual: {:?}", a, b)
		}
		Io(err: io::Error) {
			from()
			cause(err)
			description(err.description())
		}
	}
);

quick_error!(
	#[derive(Debug,Clone)]
	pub enum FakeError {
		Inner(desc: String, disp: String) {
			description(desc)
			display("{}", disp)
		}
	}
);

impl<'e> From<&'e io::Error> for ParseError {
	fn from(e: &io::Error) -> ParseError {
		ParseError::Io(io::Error::new(e.kind(), FakeError::Inner(e.description().to_owned(), format!("{}", e))))
	}
}

impl Clone for ParseError {
	fn clone(&self) -> ParseError {
		use self::ParseError::*;
		match *self {
			NoAt(b)	=> NoAt(b),
			NoPlus(ref p, b)	=> NoPlus(p.clone(), b),
			Incomplete(ref p)	=> Incomplete(p.clone()),
			LengthMismatch(ref s, ref q)	=> LengthMismatch(s.clone(), q.clone()),
			Io(ref e)	=> e.into(),
		}
	}
}

macro_rules! try_some(($expr:expr) => {
	match $expr {
		Ok(v)  => v,
		Err(e) => return Some(Err(e.into())),
	}
});

pub struct FastqReader<R>(pub R);

impl<R: BufRead> Iterator for FastqReader<R> {
	type Item = Result<Record, ParseError>;
	
	fn next(&mut self) -> Option<Result<Record, ParseError>> {
		let &mut FastqReader(ref mut it) = self;
		
		let mut at = [0];
		if try_some!(it.read(&mut at)) == 0 { return None };
		if at[0] != b'@' { return Some(Err(ParseError::NoAt(at[0]))) };
		
		let mut header = try_some!(read_line_without_nl(it, || "@<nothing>".to_owned()));
		
		let desc = header.split_whitespace().next_back().map(|desc| desc.to_owned());
		if let Some(ref desc) = desc {
			let l = header.len();
			header.truncate(l - desc.len() - 1);
		}
		
		let seq = try_some!(read_line_without_nl(it, || format!("@{}\n<nothing>\n+\n<nothing>", header)));
	
		let mut qual_head = String::new();
		if let Err(e) = it.read_line(&mut qual_head) { return Some(Err(e.into())) };
		if !qual_head.starts_with('+') { return Some(Err(ParseError::NoPlus(format!("@{}\n{}", header, seq), qual_head.bytes().next().unwrap()))) }
		
		let qual = try_some!(read_line_without_nl(it, || format!("@{}\n{}\n+\n<nothing>", header, seq)));
		
		Some(if seq.len() == qual.len() {
			Ok(Record::from_strings(header, desc, seq, qual))
		} else {
			Err(ParseError::LengthMismatch(seq, qual))
		})
	}
}

#[inline]
fn read_line_without_nl<R, F>(r: &mut R, f: F) -> Result<String, ParseError> where R: BufRead, F: Fn() -> String {
	let mut string = String::new();
	try!(r.read_line(&mut string));
	if string.len() <= 1 { return Err(ParseError::Incomplete(f())) }
	assert!(string.ends_with('\n'));
	string.pop();
	Ok(string)
}
