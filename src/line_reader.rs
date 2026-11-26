//! Buffered line reading without a fresh allocation per line.

use std::fs::File;
use std::io::{BufRead, BufReader};

/// Line reader equivalent to `BufRead::lines` (same terminator handling and
/// UTF-8 errors) that reuses one buffer instead of allocating per line.
pub struct LineReader {
    reader: BufReader<File>,
    buffer: Vec<u8>,
}

impl LineReader {
    pub fn open(path: &str) -> std::io::Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(1 << 20, File::open(path)?),
            buffer: Vec::new(),
        })
    }

    /// Open `path` positioned at byte `offset`.
    pub fn open_at(path: &str, offset: u64) -> std::io::Result<Self> {
        use std::io::{Seek, SeekFrom};
        let mut file = File::open(path)?;
        if offset != 0 {
            file.seek(SeekFrom::Start(offset))?;
        }
        Ok(Self {
            reader: BufReader::with_capacity(1 << 20, file),
            buffer: Vec::new(),
        })
    }

    pub fn next_line(&mut self) -> Option<std::io::Result<&str>> {
        self.buffer.clear();
        match self.reader.read_until(b'\n', &mut self.buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(error) => return Some(Err(error)),
        }
        let mut line = &self.buffer[..];
        if let Some(rest) = line.strip_suffix(b"\n") {
            line = rest.strip_suffix(b"\r").unwrap_or(rest);
        }
        Some(std::str::from_utf8(line).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "stream did not contain valid UTF-8",
            )
        }))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_reader_matches_bufread_lines() -> std::io::Result<()> {
        let dir = tempfile::tempdir()?;
        for content in [
            &b""[..],
            b"\n",
            b"a\r\nb\nc",
            b"x\r\r\n\ry\r",
            b"h,e\n\n,,id,ann\r\nlast,line,id2,ann2",
            b"ok\n\xff\xfe\nafter\n",
            b" \t padded \t\r\n##FASTA\n>c1 desc\nACGT\n",
        ] {
            let path = dir.path().join("lines.txt");
            std::fs::write(&path, content)?;
            let path = path.to_str().unwrap();
            let expected: Vec<_> = BufReader::new(File::open(path)?)
                .lines()
                .map(|line| line.map_err(|error| error.kind()))
                .collect();
            let mut lines = LineReader::open(path)?;
            let mut actual = Vec::new();
            while let Some(line) = lines.next_line() {
                actual.push(line.map(str::to_owned).map_err(|error| error.kind()));
            }
            assert_eq!(actual, expected);
        }
        Ok(())
    }
}
