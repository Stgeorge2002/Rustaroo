//! Buffered line reading without a fresh allocation per line.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Take};

/// Line reader equivalent to `BufRead::lines` (same terminator handling and
/// UTF-8 errors) that reuses one buffer instead of allocating per line.
pub struct LineReader {
    reader: BufReader<Take<File>>,
    buffer: Vec<u8>,
    /// Byte offset just past the last line returned.
    position: u64,
}

impl LineReader {
    pub fn open(path: &str) -> std::io::Result<Self> {
        Self::open_range(path, 0, u64::MAX)
    }

    /// Open `path` positioned at byte `offset`.
    pub fn open_at(path: &str, offset: u64) -> std::io::Result<Self> {
        Self::open_range(path, offset, u64::MAX)
    }

    /// Read only the lines in bytes `start..end` of `path`. With both ends
    /// on line boundaries (see `line_chunks`), the lines of consecutive
    /// ranges are exactly the lines of the whole range.
    pub fn open_range(path: &str, start: u64, end: u64) -> std::io::Result<Self> {
        let mut file = File::open(path)?;
        if start != 0 {
            file.seek(SeekFrom::Start(start))?;
        }
        Ok(Self {
            reader: BufReader::with_capacity(1 << 20, file.take(end.saturating_sub(start))),
            buffer: Vec::new(),
            position: start,
        })
    }

    /// Byte offset just past the last line returned (its terminator included).
    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn next_line(&mut self) -> Option<std::io::Result<&str>> {
        self.buffer.clear();
        match self.reader.read_until(b'\n', &mut self.buffer) {
            Ok(0) => return None,
            Ok(read) => self.position += read as u64,
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

/// Split bytes `start..` of `path` into at most `parts` ranges that each
/// begin at a line start (the byte after a `\n`, or `start` itself), for
/// reading with `LineReader::open_range` in parallel.
pub fn line_chunks(path: &str, start: u64, parts: usize) -> std::io::Result<Vec<(u64, u64)>> {
    let mut file = File::open(path)?;
    let len = file.metadata()?.len();
    if start >= len || parts <= 1 {
        return Ok(vec![(start, len.max(start))]);
    }
    let mut bounds = vec![start];
    let mut probe = [0u8; 1 << 16];
    for part in 1..parts {
        let target = start + ((len - start) as u128 * part as u128 / parts as u128) as u64;
        // The first line start at or after `target`: just past the first
        // newline at or after `target - 1`.
        let mut at = target.max(start + 1) - 1;
        if at < *bounds.last().unwrap() {
            at = *bounds.last().unwrap();
        }
        file.seek(SeekFrom::Start(at))?;
        let boundary = loop {
            let read = file.read(&mut probe)?;
            if read == 0 {
                break len;
            }
            if let Some(index) = memchr::memchr(b'\n', &probe[..read]) {
                break at + index as u64 + 1;
            }
            at += read as u64;
        };
        if boundary > *bounds.last().unwrap() && boundary < len {
            bounds.push(boundary);
        }
    }
    bounds.push(len);
    Ok(bounds.windows(2).map(|pair| (pair[0], pair[1])).collect())
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
            // Any chunking gives the same lines, in order.
            for parts in [1, 2, 3, 5, 64] {
                for start in [0u64, 1] {
                    let start = start.min(content.len() as u64);
                    let mut whole = LineReader::open_at(path, start)?;
                    let mut expected = Vec::new();
                    while let Some(line) = whole.next_line() {
                        expected.push(line.map(str::to_owned).map_err(|error| error.kind()));
                    }
                    let mut chunked = Vec::new();
                    let chunks = line_chunks(path, start, parts)?;
                    assert_eq!(chunks.first().unwrap().0, start);
                    assert_eq!(
                        chunks.last().unwrap().1,
                        content.len().max(start as usize) as u64
                    );
                    for (from, to) in chunks {
                        let mut lines = LineReader::open_range(path, from, to)?;
                        while let Some(line) = lines.next_line() {
                            chunked.push(line.map(str::to_owned).map_err(|error| error.kind()));
                        }
                        assert_eq!(lines.position(), to);
                    }
                    assert_eq!(chunked, expected, "{content:?} {parts} {start}");
                }
            }
        }
        Ok(())
    }
}
