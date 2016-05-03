#[macro_use] extern crate quick_error;

pub mod fancy_parser;
pub mod unfancy_parser;

trait Record {
	/// Create a new, empty FastQ record.
	fn new() -> Self;
	
	/// Check if record is empty.
	fn is_empty(&self) -> bool;
	
	/// Check validity of FastQ record.
	fn check(&self) -> Result<(), &str>;
	
	/// Return the id of the record.
	fn id(&self) -> Option<&str>;
	/// Return descriptions if present.
	fn desc(&self) -> Option<&str>;
	/// Return the sequence of the record.
	fn seq(&self) -> &[u8];
	/// Return the base qualities of the record.
	fn qual(&self) -> &[u8];
	
	/// Clear the record.
	fn clear(&mut self);
}
